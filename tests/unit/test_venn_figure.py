"""The drawn Venn reports the same numbers the object stores.

A figure computed independently of the stored membership is a figure that can disagree with the
table beside it, and the disagreement would be invisible -- both look plausible. So the drawing
takes the SAME exclusive-region computation the TSV does, and this pins that they cannot drift.

matplotlib-venn labels its regions by binary id: '10' is the left set only, '01' the right only,
'11' the intersection. Those labels are what a reader actually counts off the figure, so they
are what is asserted, rather than any internal state.

Two sets are drawn area-proportionally, which is the reason for using matplotlib-venn at 2-3
sets rather than the fixed ellipses needed at 4+: the overlap is legible without reading a
number at all.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from trnagraph.modules import plotsVenn


SETS = {'Fragment (u60)': ['t1', 't2', 't3'], 'Full length (o60)': ['t2', 't3', 't4']}


def _labels(sets):
    fig, ax = plt.subplots()
    try:
        diagram = plotsVenn.draw_venn(ax, sets)
        return {region: diagram.get_label_by_id(region).get_text()
                for region in ('10', '01', '11')}
    finally:
        plt.close(fig)


def test_drawn_counts_match_the_exclusive_regions():
    regions = plotsVenn.exclusive_regions(SETS)
    labels = _labels(SETS)

    assert labels['10'] == str(len(regions['Fragment (u60)']))
    assert labels['01'] == str(len(regions['Full length (o60)']))
    assert labels['11'] == str(len(regions['Fragment (u60) & Full length (o60)']))


def test_the_counts_are_the_ones_a_reader_would_expect():
    """Stated independently of the implementation, so the test is not merely self-consistent."""
    labels = _labels(SETS)

    assert (labels['10'], labels['01'], labels['11']) == ('1', '1', '2')


def test_a_set_with_no_unique_members_still_draws():
    """A fully-nested pair is a real result -- every fragment tRNA also seen full length -- and
    must not raise."""
    labels = _labels({'A': ['t1', 't2'], 'B': ['t1', 't2', 't3']})

    assert labels['10'] == '0'
    assert labels['11'] == '2'
