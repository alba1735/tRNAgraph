"""A declared set's own variant wins; --variant fills the gaps and is otherwise ignored.

`graph` resolves exactly one --variant up front and every other plot type is built from it. The
Venn is the deliberate exception: its whole purpose can be to contrast u60 with o60, so a
diagram whose sets name their own variants cannot be forced onto one.

Rather than make --variant an error alongside such a declaration, it is IGNORED for the sets
that name a variant and used as the default for those that do not -- so a config crossing
timepoint with read type, and saying nothing about variants, still follows the flag. The run
says which happened, because a silently disregarded flag is worse than a refused one.
"""
import anndata as ad
import numpy as np
import pandas as pd

from trnagraph.modules import plotsVenn
from trnagraph.modules.toolsSchemas import MultivariateConfig


def _adata():
    rows = []
    for level in ('D0', 'D70'):
        for replicate in range(2):
            rows.append({'trna': 'tRNA-Ala-AGC-1', 'sample': f'{level}_{replicate}',
                         'timepoint': level, 'nreads_total_norm': 100.0})
    obs = pd.DataFrame(rows)
    obs.index = [f'o{i}' for i in range(len(obs))]
    adata = ad.AnnData(X=np.zeros((len(obs), 1), dtype='float32'), obs=obs)
    adata.uns['size_splits'] = {'u60': {}, 'o60': {}}
    for tag in ('u60', 'o60'):
        adata.obsm[f'size_split_{tag}'] = pd.DataFrame(
            {'nreads_total_norm': obs['nreads_total_norm'].values}, index=obs.index)
    return adata


def _block(sets):
    return MultivariateConfig.model_validate(
        {'grouping': 'timepoint', 'venn': [{'name': 'v', 'sets': sets}]})


def test_a_declared_variant_wins_over_the_flag():
    block = _block([{'level': 'D0', 'variant': 'norm:u60'},
                    {'level': 'D0', 'variant': 'norm:o60'}])

    plans = plotsVenn.declared_venn_plans(_adata(), block, variant_tag='full')

    assert [s.tag for s in plans[0].sets] == ['u60', 'o60']


def test_a_set_without_a_variant_takes_the_flags():
    block = _block([{'level': 'D0'}, {'level': 'D70'}])

    plans = plotsVenn.declared_venn_plans(_adata(), block, variant_tag='u60')

    assert [s.tag for s in plans[0].sets] == ['u60', 'u60']


def test_the_run_reports_that_the_flag_was_overridden(caplog):
    block = _block([{'level': 'D0', 'variant': 'norm:u60'},
                    {'level': 'D0', 'variant': 'norm:o60'}])

    with caplog.at_level('INFO'):
        plotsVenn.declared_venn_plans(_adata(), block, variant_tag='full')

    assert any('variant' in record.message.lower() for record in caplog.records), \
        'a disregarded flag has to be stated'


def test_no_report_when_nothing_was_overridden(caplog):
    block = _block([{'level': 'D0'}, {'level': 'D70'}])

    with caplog.at_level('INFO'):
        plotsVenn.declared_venn_plans(_adata(), block, variant_tag='u60')

    assert not [r for r in caplog.records if 'ignor' in r.message.lower()]
