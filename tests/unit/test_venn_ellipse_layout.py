"""Four- and five-set diagrams are drawn in-repo, and must not lose a region.

There is no maintained library for this. `venn`, the package the roadmap named, was last
released in 2018; `matplotlib-venn` stops at three sets; `matplotlib_set_diagrams` draws Euler
diagrams whose own documentation warns that at four or more sets "it can be impossible to
arrange the circles such that all non-empty subsets are shown" -- a silently missing intersection
is a scientific error, not a cosmetic one. So the ellipse arrangement lives here, and this is the
test that makes it trustworthy.

Region labels are placed NUMERICALLY rather than at hardcoded coordinates: the plane is sampled,
each sample point is assigned to the combination of ellipses containing it, and a label goes at
that region's centroid. That means the drawing cannot disagree with the geometry, and -- the
part that matters -- a region with no representable area is DETECTED rather than silently
dropped, so it can be reported instead of vanishing.

Counts come from exclusive_regions(), the same computation feeding the TSV and the 2-set path.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from trnagraph.modules import plotsVenn


def _four_sets():
    """One feature per region, so all 15 regions of a 4-set diagram are non-empty."""
    labels = ['A', 'B', 'C', 'D']
    sets = {label: [] for label in labels}
    for mask in range(1, 16):
        members = [labels[i] for i in range(4) if mask & (1 << i)]
        feature = '+'.join(members)
        for label in members:
            sets[label].append(feature)
    return sets


def _draw(sets):
    fig, ax = plt.subplots()
    try:
        return plotsVenn.draw_ellipse_venn(ax, sets), ax
    finally:
        plt.close(fig)


def test_every_region_of_a_four_set_diagram_is_placed():
    sets = _four_sets()

    placed, unplaceable = _draw(sets)[0]

    assert len(placed) == 15, 'all 2**4 - 1 regions carry a label'
    assert unplaceable == [], 'no region was dropped'


def test_the_drawn_counts_match_the_exclusive_regions():
    sets = _four_sets()
    regions = plotsVenn.exclusive_regions(sets)

    placed, _ = _draw(sets)[0]

    assert {region: count for region, (count, _, _) in placed.items()} == \
           {region: len(members) for region, members in regions.items()}


def test_an_empty_region_gets_no_label():
    """15 regions exist geometrically; only the populated ones are worth drawing."""
    sets = {'A': ['x'], 'B': ['x'], 'C': ['x'], 'D': ['x']}

    placed, _ = _draw(sets)[0]

    assert list(placed) == ['A & B & C & D']


def test_five_sets_place_every_populated_region():
    labels = ['A', 'B', 'C', 'D', 'E']
    sets = {label: [] for label in labels}
    for mask in range(1, 32):
        members = [labels[i] for i in range(5) if mask & (1 << i)]
        for label in members:
            sets[label].append('+'.join(members))

    placed, unplaceable = _draw(sets)[0]

    assert len(placed) + len(unplaceable) == 31
    assert len(unplaceable) == 0, f'unrepresentable regions: {unplaceable}'
