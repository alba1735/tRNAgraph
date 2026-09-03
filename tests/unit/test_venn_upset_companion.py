"""Every complex Venn gets an UpSet plot reporting the same numbers.

Two reasons it exists. A four- or five-circle Venn is drawable but hard to READ -- finding the
region for a particular combination means tracing overlapping ellipses -- while an UpSet plot
makes each intersection a labelled row. And beyond five sets the Venn cannot be drawn honestly at
all, so UpSet is the only representation left.

Because they are two renderings of one computation, the risk is that they drift. The assertion
that matters is therefore not that the plot is produced, but that its intersection sizes equal
the Venn's regions exactly.

UpSet is cited work (Lex et al. 2014, IEEE TVCG 20(12):1983-1992), which is also why it was
chosen over drawing something bespoke.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from trnagraph.modules import plotsVenn


SETS = {
    'D0 o60': ['a', 'b', 'c', 'd'],
    'D70 o60': ['b', 'c', 'e'],
    'D0 u60': ['c', 'd', 'f'],
    'D70 u60': ['c', 'g'],
}


def test_upset_intersections_equal_the_venn_regions():
    regions = plotsVenn.exclusive_regions(SETS)

    sizes = plotsVenn.upset_intersection_sizes(SETS)

    assert sizes == {region: len(members) for region, members in regions.items()}


def test_every_feature_is_counted_once():
    sizes = plotsVenn.upset_intersection_sizes(SETS)

    assert sum(sizes.values()) == len(set().union(*(set(v) for v in SETS.values())))


def test_the_plot_is_produced_for_a_size_no_venn_can_draw():
    """Six sets: no honest Venn exists, and UpSet is the whole answer."""
    sets = {f'S{i}': [f'f{i}', 'shared'] for i in range(6)}

    fig = plotsVenn.draw_upset(sets, title='six sets')

    try:
        assert fig is not None
        assert plotsVenn.venn_layout(len(sets)) == 'upset_only'
    finally:
        plt.close(fig)
