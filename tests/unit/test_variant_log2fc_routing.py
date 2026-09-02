"""A split variant's fold changes are cached under that split, never over the full variant's.

The full variant's fits live in adata.uns['log2FC']; a split's live in
adata.uns['size_splits'][tag]['log2FC']. build_variant_view()'s docstring is explicit that the
resolved copy must never be written back to the original path, because for a split that would
overwrite the real full/default data with the split's overlaid values. The cross-variant
accessor writes results itself, so it has to honour that split by hand -- and it is the one
piece of this design that can silently corrupt an object if it gets it wrong.

Counts are drawn from a negative binomial, with half the features flat: constant values collapse
PyDESeq2's dispersion estimates, and a dataset where everything moves together is read as a
library-size difference and normalized away. Same reasoning as test_compare_plots.py's fixture.
"""
import anndata as ad
import numpy as np
import pandas as pd
import pytest

from trnagraph.modules import toolsTG


GROUPS = ['ctrl', 'treat']
TRNAS = [f'tRNA-Ala-AGC-{i}' for i in range(1, 7)]
REPLICATES = 3


def _values(rng, scale, moving, group_index):
    mean = scale * (4 ** group_index) if moving else scale
    return float(rng.negative_binomial(n=20, p=20 / (20 + mean)))


def _adata():
    """One object carrying a full variant and two split tags, each with different counts."""
    rng = np.random.default_rng(0)
    rows, u60, o60 = [], [], []
    for g, group in enumerate(GROUPS):
        for replicate in range(REPLICATES):
            for i, trna in enumerate(TRNAS):
                moving = i < len(TRNAS) // 2
                full = _values(rng, 400.0, moving, g)
                rows.append({'trna': trna, 'sample': f'{group}_{replicate}', 'group': group,
                             'nreads_total_norm': full, 'nreads_total_raw': full})
                under = _values(rng, 100.0, moving, g)
                over = _values(rng, 250.0, moving, g)
                u60.append({'nreads_total_norm': under, 'nreads_total_raw': under})
                o60.append({'nreads_total_norm': over, 'nreads_total_raw': over})
    obs = pd.DataFrame(rows)
    obs.index = [f'o{i}' for i in range(len(obs))]
    adata = ad.AnnData(X=np.zeros((len(obs), 1), dtype='float32'), obs=obs,
                       var=pd.DataFrame({'coverage': ['uniquecoverage']}, index=['v0']))
    adata.obsm['size_split_u60'] = pd.DataFrame(u60, index=obs.index)
    adata.obsm['size_split_o60'] = pd.DataFrame(o60, index=obs.index)
    adata.uns['size_splits'] = {'u60': {}, 'o60': {}}
    return adata


def _compute(adata, tag):
    return toolsTG.variant_log2fc(adata, tag, 'group', 'nreads_total_norm',
                                  readcount_cutoff=0, shrink='none')


def test_a_split_fit_is_cached_under_that_split():
    adata = _adata()

    _compute(adata, 'u60')

    assert 'log2FC' in adata.uns['size_splits']['u60']
    assert 'default' in adata.uns['size_splits']['u60']['log2FC']


def test_a_split_fit_never_lands_in_the_full_variants_cache():
    """The corruption case: a split's values written where the real data belongs."""
    adata = _adata()

    _compute(adata, 'u60')

    assert not adata.uns.get('log2FC'), "the full variant's cache must be untouched"


def test_two_tags_do_not_disturb_each_other():
    adata = _adata()

    _compute(adata, 'u60')
    _compute(adata, 'o60')

    assert 'log2FC' in adata.uns['size_splits']['u60']
    assert 'log2FC' in adata.uns['size_splits']['o60']


def test_the_full_variant_still_uses_the_unsuffixed_location():
    adata = _adata()

    _compute(adata, 'full')

    assert 'default' in adata.uns['log2FC']
    assert 'log2FC' not in adata.uns['size_splits'].get('u60', {})


def test_each_tag_fits_its_own_counts():
    """Different underlying counts must give different fold changes -- otherwise the accessor
    is reading one variant's numbers for every tag, which no routing assertion would catch."""
    adata = _adata()

    under, _ = _compute(adata, 'u60')
    over, _ = _compute(adata, 'o60')

    column = under.columns[0]
    assert not np.allclose(under[column].values, over[column].values, equal_nan=True)


def test_a_cached_fit_is_returned_rather_than_recomputed():
    """Seeded with a sentinel rather than mocked: what matters is that the accessor READS the
    tag's cache, which is observable by putting a recognisable frame there and asking for it."""
    adata = _adata()
    sentinel = pd.DataFrame({'log2_ctrl-treat': [42.0], 'pval_ctrl-treat': [0.5]},
                            index=['tRNA-Ala-AGC-1'])
    adata.uns['size_splits']['u60']['log2FC'] = {
        'default': {'group': {'nreads_total_norm': {'0': {'df': sentinel, 'shrink': 'none'}},
                              'pairs': [('ctrl', 'treat')]}}}

    df, pairs = _compute(adata, 'u60')

    assert df.equals(sentinel), 'the cached frame was returned'
    assert pairs == [('ctrl', 'treat')]


def test_the_accessor_never_resolves_a_full_variant_view(monkeypatch):
    """build_variant_view() copies the whole object -- roughly 460MB on the human dataset.
    Calling it per tag is the failure mode this accessor exists to avoid, so the design is
    pinned rather than left to reviewer memory."""
    def explode(*args, **kwargs):
        raise AssertionError('build_variant_view must not be called by the cross-variant path')
    monkeypatch.setattr(toolsTG, 'build_variant_view', explode)

    _compute(_adata(), 'u60')
