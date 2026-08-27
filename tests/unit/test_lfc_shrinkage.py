"""apeGLM shrinkage of the log2 fold changes.

tRNAgraph stored raw maximum-likelihood LFCs while tRAX shrinks its own via DESeq2's
`betaPrior=TRUE` (`analyzecounts.R:104`) -- an unrecorded divergence on the headline DE
statistic. PyDESeq2 offers only apeGLM (Zhu, Ibrahim & Love 2019,
doi:10.1093/bioinformatics/bty895), which shrinks a design-matrix *coefficient* rather than
an arbitrary contrast, so a pair can only be shrunk from a fit whose reference level is that
pair's own baseline. Pairs are therefore grouped by baseline and each distinct baseline costs
one fit.

Measured on the real hg38 object: p-values are bit-identical for pairs sharing the default
reference, and differ by at most 9e-3 with zero significance flips for the pair that needed
its own refit -- that difference comes from refitting, not from shrinkage.
"""
import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")

from trnagraph.modules.toolsTG import adataLog2FC


def _adata(groups=('A', 'A', 'A', 'B', 'B', 'B', 'C', 'C', 'C'), n_trna=25, seed=0):
    """Nine samples over three groups, with a few features given a real effect."""
    rng = np.random.default_rng(seed)
    trnas = [f'tRNA-Test-{i}' for i in range(n_trna)]
    rows = []
    for sample_index, group in enumerate(groups):
        sample = f's{sample_index}'
        for trna_index, trna in enumerate(trnas):
            base = 200 if trna_index % 5 else 8
            lift = 12 if (group == 'C' and trna_index % 5 == 0) else 1
            raw = float(rng.poisson(base * lift))
            rows.append({'trna': trna, 'sample': sample, 'group': group,
                         'nreads_total_unique_raw': raw,
                         'nreads_total_unique_norm': raw})
    obs = pd.DataFrame(rows)
    obs.index = [f'{r.sample}_{r.trna}' for r in obs.itertuples()]
    return anndata.AnnData(X=np.zeros((len(obs), 1)), obs=obs)


def _run(shrink, **kwargs):
    adata = kwargs.pop('adata', None) or _adata()
    return adataLog2FC(adata, compare='group', readtype='nreads_total_unique_norm',
                       readcount_cutoff=0, lfc_shrink=shrink, **kwargs).log2fc_df()


def test_shrinkage_is_on_by_default():
    """tRAX shrinks; matching it is the default, not an opt-in."""
    adata = _adata()
    assert adataLog2FC(adata, compare='group', readtype='nreads_total_unique_norm').lfc_shrink


def test_shrinkage_pulls_estimates_toward_zero():
    mle, _ = _run(False)
    shrunk, _ = _run(True)

    columns = [c for c in mle.columns if c.startswith('log2_')]
    assert shrunk[columns].abs().median().median() <= mle[columns].abs().median().median()


def test_every_pair_gets_a_fold_change_and_a_pvalue():
    """The bug this caught: the pvalue assignment sat outside the per-pair loop, so only one
    column per reference group was populated and the rest were silently all-NaN."""
    for shrink in (False, True):
        df, pairs = _run(shrink)
        assert len(pairs) == 3
        for a, b in pairs:
            assert df[f'log2_{a}-{b}'].notna().any(), f'log2_{a}-{b} all NaN (shrink={shrink})'
            assert df[f'pval_{a}-{b}'].notna().any(), f'pval_{a}-{b} all NaN (shrink={shrink})'


def test_shrinkage_does_not_flip_significance_calls():
    mle, pairs = _run(False)
    shrunk, _ = _run(True)

    for a, b in pairs:
        column = f'pval_{a}-{b}'
        joined = pd.concat([mle[column].rename('m'), shrunk[column].rename('s')], axis=1).dropna()
        if joined.empty:
            continue
        assert (((joined.m < 0.05) != (joined.s < 0.05)).sum()) == 0


def test_pairs_and_shape_are_unchanged_by_shrinkage():
    mle, mle_pairs = _run(False)
    shrunk, shrunk_pairs = _run(True)

    assert mle_pairs == shrunk_pairs
    assert list(mle.columns) == list(shrunk.columns)
    assert list(mle.index) == list(shrunk.index)


# --- caching -------------------------------------------------------------------------

def test_cache_records_the_shrinkage_state():
    adata = _adata()
    calc = adataLog2FC(adata, compare='group', readtype='nreads_total_unique_norm',
                       readcount_cutoff=0, lfc_shrink=True)
    calc.main()

    entry = adata.uns['log2FC']['default']['group']['nreads_total_unique_norm']['0']
    assert entry['lfc_shrink'] is True


def test_switching_shrinkage_invalidates_the_cache():
    adata = _adata()
    adataLog2FC(adata, compare='group', readtype='nreads_total_unique_norm',
                readcount_cutoff=0, lfc_shrink=False).main()

    second = adataLog2FC(adata, compare='group', readtype='nreads_total_unique_norm',
                         readcount_cutoff=0, lfc_shrink=True)
    second.main()

    assert second.computed_fresh, 'a shrunk run must not reuse an unshrunken cached frame'


def test_a_cache_written_before_shrinkage_existed_is_treated_as_unshrunken():
    """Objects built before apeGLM landed carry no marker; they must not be served as shrunk."""
    adata = _adata()
    calc = adataLog2FC(adata, compare='group', readtype='nreads_total_unique_norm',
                       readcount_cutoff=0, lfc_shrink=True)
    calc.main()
    entry = adata.uns['log2FC']['default']['group']['nreads_total_unique_norm']['0']
    del entry['lfc_shrink']

    again = adataLog2FC(adata, compare='group', readtype='nreads_total_unique_norm',
                        readcount_cutoff=0, lfc_shrink=True)
    again.main()

    assert again.computed_fresh


def test_matching_shrinkage_state_hits_the_cache():
    adata = _adata()
    adataLog2FC(adata, compare='group', readtype='nreads_total_unique_norm',
                readcount_cutoff=0, lfc_shrink=True).main()

    second = adataLog2FC(adata, compare='group', readtype='nreads_total_unique_norm',
                         readcount_cutoff=0, lfc_shrink=True)
    second.main()

    assert not second.computed_fresh
