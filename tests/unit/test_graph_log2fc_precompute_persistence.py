"""Regression tests for adataGraph.py's pre-pool log2FC precompute/write-back step. This runs
BEFORE the graph-type multiprocessing.Pool is created ("to prevent saving issues later if
multiprocessing is used"), and is meant to persist any freshly-computed combo back to the h5ad so
a later `graph` run (e.g. regenerating figures) hits the cache instead of recomputing.

Bug: the change-detection snapshotted the "before" state with a shallow `dict.copy()`, then
compared it to the (mutated in place) "after" state. Since adataLog2FC.main() mutates the SAME
nested dict objects the shallow copy shares references to, adding a new readtype/cutoff under an
already-existing config_name/compare path (the common case -- 'default'/'group' by default) never
registers as a difference, so the write-back essentially never fires even when real computation
just happened. Confirmed against the real hg38 server run: only 2 of 5 default --diffrts
readtypes were ever actually persisted, with zero write-back log messages anywhere."""
import logging
from types import SimpleNamespace
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from trnagraph.modules.adataGraph import anndataGrapher
from trnagraph.modules.toolsSchemas import VariantTag

READCOUNT_CUTOFF = 80


def _make_adata(n_per_group=4, precomputed_readtypes=('total_unique',)):
    """Builds obs with nreads_total_unique_{raw,norm} and nreads_total_{raw,norm} (the two
    readtypes the volcano overview always needs), and pre-seeds uns['log2FC'] with whichever of
    them are listed in precomputed_readtypes -- reproducing "some combos already cached under an
    existing config_name/compare path, others not" exactly as in the real object."""
    rng = np.random.default_rng(0)
    groups = ['A', 'B']
    samples = [f'{g}_rep{i}' for g in groups for i in range(n_per_group)]
    sample_group = {s: s.rsplit('_rep', 1)[0] for s in samples}
    trna_means = {'trnaHigh': {'A': 800, 'B': 50}}
    trna_means.update({f'trnaFlat{i}': {'A': 400, 'B': 400} for i in range(6)})

    rows = []
    for trna, group_means in trna_means.items():
        for sample in samples:
            group = sample_group[sample]
            mean = group_means[group]
            raw = rng.negative_binomial(n=10, p=10 / (10 + mean))
            rows.append({
                'trna': trna, 'sample': sample, 'group': group,
                'nreads_total_unique_raw': raw, 'nreads_total_unique_norm': raw,
                'nreads_total_raw': raw, 'nreads_total_norm': raw,
            })
    obs = pd.DataFrame(rows)
    obs.index = [f"{r.trna}_{r['sample']}" for _, r in obs.iterrows()]
    adata = ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs)

    from trnagraph.modules.toolsTG import adataLog2FC
    log2fc_dict = {}
    for rt in precomputed_readtypes:
        adataLog2FC(adata, 'group', f'nreads_{rt}_norm', readcount_cutoff=READCOUNT_CUTOFF, config_name='default', overwrite=True).main()
    return adata


def _make_grapher(adata, graphtypes=('volcano',), diffrts=('total_unique',), regen_uns=False, threads=4):
    grapher = anndataGrapher.__new__(anndataGrapher)
    grapher.logger = logging.getLogger('trnagraph.modules.adataGraph')
    grapher.adata = adata
    grapher.adata_original = adata
    grapher.config_name = 'default'
    grapher.variant_spec = VariantTag(raw='norm:full', norm='norm', tag='full')
    grapher.args = SimpleNamespace(
        anndata='fake.h5ad', graphtypes=list(graphtypes), heatgrp='group', volgrp='group', diffrts=list(diffrts),
        heatcutoff=READCOUNT_CUTOFF, volcutoff=READCOUNT_CUTOFF, regen_uns=regen_uns, threads=threads,
    )
    return grapher


def test_persists_when_a_new_readtype_is_computed_under_an_existing_config_path():
    """The exact bug scenario: 'default'/'group' already exists in uns['log2FC'] (from
    total_unique being precomputed), and the overview loop needs to ALSO compute 'total' under
    that same existing path -- this must be detected as a real change and written back."""
    adata = _make_adata(precomputed_readtypes=('total_unique',))  # 'total' is NOT precomputed
    grapher = _make_grapher(adata)

    with patch.object(adata, 'write') as mock_write:
        grapher._precompute_and_persist_log2fc()

    mock_write.assert_called_once()
    assert 'nreads_total_norm' in adata.uns['log2FC']['default']['group']


def test_does_not_write_when_everything_is_already_cached():
    adata = _make_adata(precomputed_readtypes=('total_unique', 'total'))  # both already cached
    grapher = _make_grapher(adata)

    with patch.object(adata, 'write') as mock_write:
        grapher._precompute_and_persist_log2fc()

    mock_write.assert_not_called()


def test_precompute_passes_threads_as_n_cpus():
    """This precompute runs before the graph-type pool exists, so it's safe to use real
    parallelism here -- it should use the command's own --threads value, not the safe-by-default
    n_cpus=1 that graph-time on-demand calls (inside the pool) must use."""
    adata = _make_adata(precomputed_readtypes=())
    grapher = _make_grapher(adata, threads=7)

    import trnagraph.modules.toolsTG as toolsTG_module
    with patch.object(adata, 'write'), \
         patch('trnagraph.modules.toolsTG.DeseqDataSet', wraps=toolsTG_module.DeseqDataSet) as mock_dds:
        grapher._precompute_and_persist_log2fc()

    assert mock_dds.call_args_list
    for call in mock_dds.call_args_list:
        assert call.kwargs.get('n_cpus') == 7
