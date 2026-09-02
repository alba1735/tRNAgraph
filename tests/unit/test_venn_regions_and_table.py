"""Venn membership as a table, because the figure alone is not the result.

A Venn's scientific content is WHICH tRNAs sit in each region. That is unreadable off the
picture -- a reader can see that 34 features are unique to fragments and has no way to act on
it. The table is the artefact that goes into a paper; the figure is how it is presented.

Regions are EXCLUSIVE: a feature belongs to exactly one, named by the full combination of sets
containing it. This is the same membership an UpSet plot encodes, which is what lets the two
render the same numbers when a complex Venn gains an UpSet companion.

The table carries a provenance header naming the object and every parameter that produced it,
so a stale file left in results/ is self-identifying rather than merely wrong -- results/ is a
convenience for readers and papers, never an input, and the object stays authoritative.
"""
import pandas as pd

from trnagraph.modules import plotsVenn


SETS = {
    'A': ['t1', 't2', 't3'],
    'B': ['t2', 't3', 't4'],
}


def test_regions_are_exclusive_and_named_by_their_combination():
    regions = plotsVenn.exclusive_regions(SETS)

    assert regions['A'] == ['t1']
    assert regions['B'] == ['t4']
    assert regions['A & B'] == ['t2', 't3']


def test_every_feature_appears_in_exactly_one_region():
    regions = plotsVenn.exclusive_regions(SETS)
    placed = [f for members in regions.values() for f in members]

    assert sorted(placed) == ['t1', 't2', 't3', 't4']
    assert len(placed) == len(set(placed)), 'no feature is counted twice'


def test_an_empty_region_is_omitted_rather_than_written_as_a_blank_row():
    regions = plotsVenn.exclusive_regions({'A': ['t1'], 'B': []})

    assert 'B' not in regions
    assert regions == {'A': ['t1']}


def test_three_sets_produce_the_seven_populated_regions():
    regions = plotsVenn.exclusive_regions(
        {'A': ['x', 'ab', 'ac', 'abc'], 'B': ['b', 'ab', 'bc', 'abc'], 'C': ['c', 'ac', 'bc', 'abc']})

    assert regions['A & B & C'] == ['abc']
    assert regions['A & B'] == ['ab']
    assert regions['A'] == ['x']
    assert len(regions) == 7


def test_the_table_carries_a_provenance_header(tmp_path):
    path = tmp_path / 'venn.tsv'

    plotsVenn.write_membership_table(path, SETS, provenance={
        'object': 'demo.h5ad', 'config': 'default', 'readtype': 'nreads_total_norm',
        'cutoff': 20})

    text = path.read_text()
    header = [line for line in text.splitlines() if line.startswith('#')]
    assert any('demo.h5ad' in line for line in header)
    assert any('cutoff' in line and '20' in line for line in header)


def test_the_table_is_readable_as_a_dataframe(tmp_path):
    """It exists to be opened by a person or pasted into a paper, so it has to parse."""
    path = tmp_path / 'venn.tsv'
    plotsVenn.write_membership_table(path, SETS, provenance={'object': 'demo.h5ad'})

    df = pd.read_csv(path, sep='\t', comment='#')

    assert list(df.columns) == ['region', 'n', 'features']
    assert set(df['region']) == {'A', 'B', 'A & B'}
    assert df.loc[df['region'] == 'A & B', 'n'].item() == 2


def _adata_with_build_dir(build_dir):
    import anndata as ad
    import numpy as np
    obs = pd.DataFrame({'trna': ['t1'], 'sample': ['s1']}, index=['o0'])
    adata = ad.AnnData(X=np.zeros((1, 1), dtype='float32'), obs=obs)
    if build_dir is not None:
        adata.uns['trnagraphruninfo'] = {'trnagraph_directory': str(build_dir),
                                         'expname': 'demo'}
    return adata


def test_the_table_lands_beside_the_other_result_files(tmp_path):
    """results/multivariate/ is a SIBLING of results/u60 and results/o60, not inside either:
    a Venn spanning two variants belongs to neither, and filing it under one would say
    something false about what produced it."""
    (tmp_path / 'results').mkdir()

    directory, message = plotsVenn.resolve_results_dir(_adata_with_build_dir(tmp_path))

    assert directory == str(tmp_path / 'results' / 'multivariate')
    assert message is None


def test_a_build_directory_that_no_longer_exists_warns_and_skips(tmp_path):
    """Not hypothetical: the demo object records a scratch path from the run that built it.
    There is no fallback-directory convention in tRNAgraph, and inventing one would scatter
    output. The object still holds the membership, so only the convenience copy is lost."""
    missing = tmp_path / 'gone'

    directory, message = plotsVenn.resolve_results_dir(_adata_with_build_dir(missing))

    assert directory is None
    assert str(missing) in message, 'names the path it looked for'


def test_an_object_without_build_provenance_warns_and_skips():
    directory, message = plotsVenn.resolve_results_dir(_adata_with_build_dir(None))

    assert directory is None
    assert 'trnagraphruninfo' in message
