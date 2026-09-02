"""The declared order picks the DESeq2 reference level, not just the legend order.

This is the claim that justifies storing order on the object instead of passing it per call.
adataLog2FC derives its contrasts from `pivot_table(columns=compare).columns` and takes
`itertools.combinations` over them, so the column order IS the contrast order, and pair[0] is
the baseline every fold change is measured against (log2FoldChange = log2(pair[1] / pair[0])).

With a plain string column that order is whatever pandas produces for the unique values; with
an ordered Categorical it is the declared one. The fixture uses D0/D7/D14/D35, where chronology
and alphabet genuinely disagree -- 'D14' sorts before 'D7'.

Counts are drawn from a negative binomial rather than set to constants: with zero within-group
variance PyDESeq2's dispersion estimates collapse and shrinkage returns nonsense, which is an
artifact of degenerate synthetic data rather than anything the pipeline does. Same reasoning as
test_compare_plots.py's fixture.
"""
import anndata as ad
import numpy as np
import pandas as pd

from trnagraph.modules import toolsTG


TIMEPOINTS = ['D0', 'D7', 'D14', 'D35']
TRNAS = [f'tRNA-Ala-AGC-{i}' for i in range(1, 7)]
REPLICATES = 3


def _make_adata(ordered):
    rng = np.random.default_rng(0)
    rows = []
    for t_index, timepoint in enumerate(TIMEPOINTS):
        for replicate in range(REPLICATES):
            for i, trna in enumerate(TRNAS):
                # Half the features move with time, half stay flat: a dataset where everything
                # moves together is read as a library-size difference and normalized away.
                mean = 100.0 * (2 ** t_index) if i < len(TRNAS) // 2 else 400.0
                value = float(rng.negative_binomial(n=20, p=20 / (20 + mean)))
                rows.append({'trna': trna, 'sample': f'{timepoint}_{replicate}',
                             'timepoint': timepoint,
                             'nreads_total_norm': value, 'nreads_total_raw': value})
    obs = pd.DataFrame(rows)
    obs.index = [f'obs{i}' for i in range(len(obs))]
    if ordered:
        toolsTG.apply_category_order(obs, {'timepoint': TIMEPOINTS})
    var = pd.DataFrame({'coverage': ['uniquecoverage']}, index=['v0'])
    return ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs, var=var)


def _pairs(ordered):
    adata = _make_adata(ordered)
    _, pairs = toolsTG.adataLog2FC(
        adata, 'timepoint', 'nreads_total_norm', readcount_cutoff=0, shrink='none'
    ).log2fc_df()
    return pairs


def test_contrasts_follow_the_declared_order():
    assert [p[0] for p in _pairs(ordered=True)][:3] == ['D0', 'D0', 'D0']
    assert _pairs(ordered=True)[0] == ('D0', 'D7'), 'the earliest timepoint is the baseline'


def test_without_a_declared_order_the_contrasts_are_alphabetical():
    """The behaviour being corrected: 'D14' sorts before 'D7', so the baseline of the first
    contrast is a timepoint in the middle of the course."""
    assert _pairs(ordered=False)[0] == ('D0', 'D14')
