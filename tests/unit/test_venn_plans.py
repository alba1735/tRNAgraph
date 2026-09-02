"""Which of the two always-drawn Venns a given object can actually support.

Both are worth drawing whenever the data allows: fragment vs full-length answers "is this tRNA
present as a fragment, as a mature molecule, or both", and 5' vs 3' answers the same question
for the two halves. Neither is declared in config -- they are drawn automatically whenever
`-g venn` runs (the plot family itself stays config-gated).

Their prerequisites differ, and neither absence is an error. Fragment vs full-length needs a
read-length split, which a plain build does not have; 5' vs 3' needs both end-specific read
types. A build lacking one must still draw the other, and must SAY which it skipped and why --
silently drawing one Venn where the user expected two is how a missing figure gets mistaken for
a biological result.
"""
import anndata as ad
import numpy as np
import pandas as pd

from trnagraph.modules import plotsVenn


def _adata(split_tags=(), readtypes=('total', 'fiveprime', 'threeprime')):
    rows = []
    for trna in ['tRNA-Ala-AGC-1', 'tRNA-Gly-GCC-1']:
        for sample in ['s1', 's2']:
            row = {'trna': trna, 'sample': sample, 'group': 'g'}
            for rt in readtypes:
                row[f'nreads_{rt}_norm'] = 100.0
                row[f'nreads_{rt}_raw'] = 100.0
            rows.append(row)
    obs = pd.DataFrame(rows)
    obs.index = [f'o{i}' for i in range(len(obs))]
    adata = ad.AnnData(X=np.zeros((len(obs), 1), dtype='float32'), obs=obs)
    if split_tags:
        adata.uns['size_splits'] = {t: {} for t in split_tags}
        for t in split_tags:
            adata.obsm[f'size_split_{t}'] = pd.DataFrame(
                {c: obs[c].values for c in obs.columns if c.startswith('nreads_')},
                index=obs.index)
    return adata


def test_a_split_build_supports_both():
    plans, skipped = plotsVenn.simple_venn_plans(_adata(split_tags=('u60', 'o60')))

    assert [p.name for p in plans] == ['fragment_vs_full_length', 'fiveprime_vs_threeprime']
    assert skipped == []


def test_an_unsplit_build_still_draws_the_readtype_venn():
    plans, skipped = plotsVenn.simple_venn_plans(_adata())

    assert [p.name for p in plans] == ['fiveprime_vs_threeprime']
    assert len(skipped) == 1
    assert 'fragment_vs_full_length' in skipped[0]
    assert 'readlengthsplit' in skipped[0] or 'split' in skipped[0], 'says what is missing'


def test_missing_end_specific_readtypes_skip_that_venn_by_name():
    plans, skipped = plotsVenn.simple_venn_plans(
        _adata(split_tags=('u60', 'o60'), readtypes=('total',)))

    assert [p.name for p in plans] == ['fragment_vs_full_length']
    assert any('fiveprime' in message for message in skipped)


def test_a_plan_names_the_two_sets_it_will_draw():
    """Each set is a (variant tag, readtype) pair -- what distinguishes the two circles."""
    plans, _ = plotsVenn.simple_venn_plans(_adata(split_tags=('u60', 'o60')))
    fragment = plans[0]

    assert [s.tag for s in fragment.sets] == ['u60', 'o60']
    assert {s.readtype for s in fragment.sets} == {'nreads_total_norm'}


def test_neither_venn_is_available_on_a_bare_object():
    plans, skipped = plotsVenn.simple_venn_plans(_adata(readtypes=('total',)))

    assert plans == []
    assert len(skipped) == 2


def _adata_with_both_bases():
    """An object carrying unique AND all-reads columns, so a basis choice is observable."""
    import anndata as ad
    rows = []
    for trna in ['tRNA-Ala-AGC-1']:
        for sample in ['s1', 's2']:
            row = {'trna': trna, 'sample': sample, 'group': 'g'}
            for rt in ('total', 'fiveprime', 'threeprime'):
                row[f'nreads_{rt}_norm'] = 100.0
                row[f'nreads_{rt}_unique_norm'] = 50.0
            rows.append(row)
    obs = pd.DataFrame(rows)
    obs.index = [f'o{i}' for i in range(len(obs))]
    return ad.AnnData(X=np.zeros((len(obs), 1), dtype='float32'), obs=obs)


def test_the_default_basis_is_unique_counts():
    """Unique (transcript-specific) counts are the project-wide default; --allreads opts out.
    The earlier fixtures carry no unique columns, so resolve_readtype falls back and the choice
    is invisible -- this one makes it observable."""
    plans, _ = plotsVenn.simple_venn_plans(_adata_with_both_bases())

    assert [s.readtype for s in plans[0].sets] == [
        'nreads_fiveprime_unique_norm', 'nreads_threeprime_unique_norm']


def test_allreads_selects_the_all_reads_columns():
    from trnagraph.modules import toolsTG

    plans, _ = plotsVenn.simple_venn_plans(_adata_with_both_bases(),
                                           read_basis=toolsTG.read_basis(True))

    assert [s.readtype for s in plans[0].sets] == [
        'nreads_fiveprime_norm', 'nreads_threeprime_norm']
