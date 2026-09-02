"""Each Venn stores its own membership; two plans must not share one slot.

uns['multivariate'] is keyed [config][analysis][entry]. Piece 2 specified that third level as
the VARIANT TAG, which is right for a per-variant analysis -- but a Venn spans variants, so
every Venn in a run would key on the same 'full' and the second would silently overwrite the
first. The stored membership then describes only whichever diagram happened to be drawn last,
while the figures and tables on disk show both, and nothing announces the disagreement.

Found by running -g venn end to end on the bacterial demo object: two figures and two tables
were written, and the object came back holding one entry.
"""
import anndata as ad
import numpy as np
import pandas as pd

from trnagraph.modules import plotsVenn
from trnagraph.modules.toolsSchemas import MultivariateConfig


def _adata(tmp_path):
    rows = []
    for trna in ['tRNA-Ala-AGC-1', 'tRNA-Gly-GCC-1']:
        for sample in ['s1', 's2']:
            row = {'trna': trna, 'sample': sample, 'group': 'g'}
            for rt in ('total', 'fiveprime', 'threeprime'):
                row[f'nreads_{rt}_norm'] = 100.0
            rows.append(row)
    obs = pd.DataFrame(rows)
    obs.index = [f'o{i}' for i in range(len(obs))]
    adata = ad.AnnData(X=np.zeros((len(obs), 1), dtype='float32'), obs=obs)
    adata.uns['size_splits'] = {'u60': {}, 'o60': {}}
    for tag in ('u60', 'o60'):
        adata.obsm[f'size_split_{tag}'] = pd.DataFrame(
            {c: obs[c].values for c in obs.columns if c.startswith('nreads_')}, index=obs.index)
    return adata


def test_both_plans_keep_their_own_membership(tmp_path):
    adata = _adata(tmp_path)

    plotsVenn.visualizer(adata, MultivariateConfig(presence_cutoff=1),
                         output=f'{tmp_path}/', config_name='default')

    entries = adata.uns['multivariate']['default']['venn']
    assert set(entries) == {'fragment_vs_full_length', 'fiveprime_vs_threeprime'}, \
        'one entry per diagram drawn, keyed by the diagram'


def test_each_entry_records_the_plan_that_produced_it(tmp_path):
    adata = _adata(tmp_path)

    plotsVenn.visualizer(adata, MultivariateConfig(presence_cutoff=1),
                         output=f'{tmp_path}/', config_name='default')

    entries = adata.uns['multivariate']['default']['venn']
    for name, entry in entries.items():
        assert entry['provenance']['plan'] == name
        assert set(entry['sets']), 'membership is stored, not just provenance'
