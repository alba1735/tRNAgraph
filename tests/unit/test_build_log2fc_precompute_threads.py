"""Regression test for threading --threads into the build-time log2FC precompute
(adataBuild.py's _adata_build_() and merge_variant_into_adata(), both of which precompute
'group'/nreads_total_{unique_,}norm at the standard cutoffs). This precompute runs once, outside
any multiprocessing.Pool, so unlike the graph-time on-demand calls inside plotsVolcano.py/
plotsHeatmap.py (which must stay at adataLog2FC's safe default of n_cpus=1 to avoid a nested-pool
deadlock), it's safe to use real parallelism here -- matching the same command's own --threads
budget rather than silently defaulting to something else."""
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd

import trnagraph.modules.toolsTG as toolsTG
from trnagraph.modules.adataBuild import _precompute_default_log2fc


def _make_adata(n_per_group=4):
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
    return ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs)


def test_precompute_default_log2fc_passes_threads_as_n_cpus():
    adata = _make_adata()

    with patch('trnagraph.modules.toolsTG.DeseqDataSet', wraps=toolsTG.DeseqDataSet) as mock_dds:
        _precompute_default_log2fc(adata, threads=5)

    assert mock_dds.call_args_list
    for call in mock_dds.call_args_list:
        assert call.kwargs.get('n_cpus') == 5


def test_precompute_default_log2fc_covers_both_overview_readtypes_at_standard_cutoffs():
    adata = _make_adata()

    _precompute_default_log2fc(adata, threads=1)

    cached = adata.uns['log2FC']['default']['group']
    assert set(cached.keys()) >= {'nreads_total_unique_norm', 'nreads_total_norm'}
    for readtype in ('nreads_total_unique_norm', 'nreads_total_norm'):
        assert set(cached[readtype].keys()) == {'20', '40', '80', '100', '200'}
