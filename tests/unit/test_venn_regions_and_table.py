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

    assert list(df.columns) == ['rank', 'region', 'n', 'features']
    assert set(df['region']) == {'A', 'B', 'A & B'}
    assert df.loc[df['region'] == 'A & B', 'n'].item() == 2


def test_the_largest_region_comes_first(tmp_path):
    """The file opens on what most features share rather than on whichever region the set order
    happened to enumerate first -- the same "most notable first" the decoupling table uses."""
    path = tmp_path / 'venn.tsv'
    plotsVenn.write_membership_table(path, SETS, provenance={})
    df = pd.read_csv(path, sep='\t', comment='#')

    assert list(df['n']) == sorted(df['n'], reverse=True)
    assert df.iloc[0]['rank'] == 1


def test_the_table_mirrors_the_figure_path_rather_than_the_build_directory():
    """These tables used to be filed under the directory recorded in the object's build
    provenance, and skipped outright when it was gone -- which it usually is: the demo object
    records a scratch path and the hg38 object a server path, so the tables went missing exactly
    when they were most wanted. They now sit in the results/ twin of the figure's own output
    path, so a reader finds the table by the route they found the plot."""
    from trnagraph.modules import toolsTG

    assert toolsTG.results_mirror('/out/graphs') == '/out/results'


