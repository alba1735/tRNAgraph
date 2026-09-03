"""One colour per set, shared by the Venn and its UpSet companion.

matplotlib-venn's stock colours are red/green/blue, which is the worst possible default for a
figure a colourblind reader has to interpret -- and the two automatic diagrams are two-set, where
red/green is exactly the pair that collapses. The palette therefore comes from tRNAgraph's own
Okabe-Ito set, ordered so the first two are the blue/orange pair that stays distinct under every
common form of colour vision deficiency.

The UpSet plot beside a complex diagram shows the same sets as rows. Reading the pair together
means matching a row to a circle, so a row's total bar and its filled matrix dots take that set's
colour. Absent dots stay grey: they mean "not in this set", and colouring them would say the
opposite.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba

from trnagraph.modules import plotsPalette, plotsVenn
from trnagraph.modules.toolsSchemas import VennSet


SETS4 = {'A': ['a', 'x'], 'B': ['b', 'x'], 'C': ['c', 'x'], 'D': ['d', 'x']}


def _sets(*labels, level=None, level_column=None, tags=None):
    """VennSets carrying only what venn_palette reads: label, level and variant tag."""
    tags = tags or ['full'] * len(labels)
    return [VennSet(label=label, readtype='nreads_total_unique_norm', tag=tag,
                    level=level, level_column=level_column)
            for label, tag in zip(labels, tags)]



def test_the_palette_starts_with_a_colourblind_safe_pair():
    colours = plotsPalette.venn_colors(2)

    assert len(colours) == 2
    assert len(set(colours)) == 2
    assert colours == list(plotsPalette.VENN_SET_COLORS[:2])


def test_the_palette_covers_the_largest_drawable_diagram():
    assert len(plotsPalette.VENN_SET_COLORS) >= plotsVenn.MAX_VENN_SETS


def test_a_two_set_venn_uses_the_palette():
    fig, ax = plt.subplots()
    try:
        diagram = plotsVenn.draw_venn(ax, {'A': ['a', 'x'], 'B': ['b', 'x']})
        drawn = to_rgba(diagram.get_patch_by_id('10').get_facecolor())
    finally:
        plt.close(fig)

    assert drawn[:3] == to_rgba(plotsPalette.venn_colors(2)[0])[:3]


def test_the_ellipse_layout_uses_the_same_palette():
    fig, ax = plt.subplots()
    try:
        plotsVenn.draw_ellipse_venn(ax, SETS4)
        faces = [to_rgba(p.get_facecolor())[:3] for p in ax.patches
                 if p.get_facecolor()[3] > 0]
    finally:
        plt.close(fig)

    expected = [to_rgba(c)[:3] for c in plotsPalette.venn_colors(4)]
    assert faces == expected


def test_upset_rows_match_their_venn_circles():
    """The point of the pairing: a row and a circle are the same set, so they are one colour."""
    fig = plotsVenn.draw_upset(SETS4, title='t')
    try:
        axes = fig._trnagraph_upset_axes
        labels = [t.get_text() for t in axes['matrix'].get_yticklabels()]
        by_label = dict(zip(SETS4, plotsPalette.venn_colors(len(SETS4))))

        bars = {label: to_rgba(patch.get_facecolor())[:3]
                for label, patch in zip(labels, axes['totals'].patches)}
        assert bars == {label: to_rgba(by_label[label])[:3] for label in labels}

        dots = axes['matrix'].collections[0]
        offsets = np.asarray(dots.get_offsets())
        faces = np.asarray(dots.get_facecolors())
        for (_, row), colour in zip(offsets, faces):
            if colour[3] < 0.9:
                continue
            assert tuple(colour[:3]) == to_rgba(by_label[labels[int(round(row))]])[:3]
    finally:
        plt.close(fig)


def test_absent_dots_stay_neutral():
    fig = plotsVenn.draw_upset(SETS4, title='t')
    try:
        faces = np.asarray(fig._trnagraph_upset_axes['matrix'].collections[0].get_facecolors())
        faint = [c for c in faces if c[3] < 0.9]
        assert faint, 'this fixture has absent cells'
        for colour in faint:
            assert colour[0] == colour[1] == colour[2], 'absent means grey, not a set colour'
    finally:
        plt.close(fig)


def test_set_labels_do_not_overlap_each_other():
    """Every ellipse label sat at (centre_x, centre_y + 0.42), and two ellipses in the 4-set
    arrangement share a centre_y -- so their labels printed on top of one another
    ("Day 70 Dasym0.n6o0rm:u60"). Placement is now derived from each ellipse's own outward
    boundary."""
    fig, ax = plt.subplots(figsize=(6, 6))
    try:
        plotsVenn.draw_ellipse_venn(ax, {'Day 0 o60': ['a'], 'Day 70 o60': ['b'],
                                         'Day 0 u60': ['c'], 'Day 70 u60': ['d']})
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        boxes = [t.get_window_extent(renderer) for t in ax.texts
                 if any(word in t.get_text() for word in ('Day 0', 'Day 70'))]
        assert len(boxes) == 4, 'one label per circle'
        for i, a in enumerate(boxes):
            for b in boxes[i + 1:]:
                assert not a.overlaps(b), 'set labels must be readable side by side'
    finally:
        plt.close(fig)


def test_a_label_sits_outside_its_own_ellipse_centre():
    """Placed at the outward extremity so it does not land on the region counts."""
    fig, ax = plt.subplots(figsize=(6, 6))
    try:
        plotsVenn.draw_ellipse_venn(ax, {'A': ['a'], 'B': ['b'], 'C': ['c'], 'D': ['d']})
        labels = {t.get_text(): t.get_position() for t in ax.texts if t.get_text() in 'ABCD'}
    finally:
        plt.close(fig)

    centre = (0.5, 0.5)
    for name, (x, y) in labels.items():
        assert (x - centre[0]) ** 2 + (y - centre[1]) ** 2 > 0.04, f'{name} is too central'


def test_the_palette_is_the_default_when_the_style_file_names_nothing():
    assert plotsVenn.venn_palette(_sets('A', 'B', 'C')) == plotsPalette.venn_colors(3)


def test_a_style_file_entry_replaces_a_set_colour():
    """`colors.venn` is keyed by SET LABEL, so a diagram's circles can be coloured to match
    whatever the rest of the figure set uses for those conditions."""
    palette = plotsVenn.venn_palette(_sets('A', 'B'), {'A': '#123456', 'B': '#654321'})

    assert palette == ['#123456', '#654321']


def test_an_unnamed_set_keeps_its_default_colour():
    """Partial overrides are allowed: the file is a substitution into the default palette, not a
    replacement for it, so naming one circle does not force you to name all five."""
    palette = plotsVenn.venn_palette(_sets('A', 'B', 'C'), {'B': '#123456'})

    default = plotsPalette.venn_colors(3)
    assert palette == [default[0], '#123456', default[2]]


def test_set_order_not_file_order_decides_the_palette():
    """The mapping is by label, so reordering the JSON cannot silently recolour the diagram."""
    colormap = {'B': '#654321', 'A': '#123456'}

    assert plotsVenn.venn_palette(_sets('A', 'B'), colormap) == ['#123456', '#654321']


def test_the_upset_companion_honours_the_same_style_entries():
    """The pairing only works while a row and its circle agree, so both read one resolution."""
    colormap = {'A': '#123456'}
    fig = plotsVenn.draw_upset(SETS4, title='t',
                               colors=plotsVenn.venn_palette(_sets(*SETS4), colormap))
    try:
        axes = fig._trnagraph_upset_axes
        labels = [t.get_text() for t in axes['matrix'].get_yticklabels()]
        bars = dict(zip(labels, axes['totals'].patches))
        assert to_rgba(bars['A'].get_facecolor())[:3] == to_rgba('#123456')[:3]
    finally:
        plt.close(fig)


# --- deriving a level's variants from one base colour -------------------------------------

def test_a_level_seen_through_two_variants_derives_a_ramp_from_its_own_colour():
    """The four-circle case crosses two timepoints with two read-length variants. Naming the
    timepoint once in `colors.timepoint` is now enough: the full-length circle keeps that
    colour and the fragment circle takes a darker, hue-rotated tier of it, so a style file no
    longer has to carry hand-lightened hexes for the split half of the diagram."""
    colors = {'timepoint': {'Day 0': '#007FFF'}}
    sets = _sets('Day 0 o60', 'Day 0 u60', level='Day 0', level_column='timepoint',
                 tags=['o60', 'u60'])

    full, fragment = plotsVenn.venn_palette(sets, None, colors=colors)

    assert to_rgba(full)[:3] == to_rgba('#007FFF')[:3], 'full length keeps the named colour'
    assert to_rgba(fragment)[:3] != to_rgba('#007FFF')[:3]


def test_the_derived_fragment_tier_is_darker_not_lighter():
    """The whole point of the change. Lightening was what the style file did by hand, and it
    washes out twice on a Venn -- once in the tint, again in the alpha the circles blend with.
    """
    from matplotlib.colors import rgb_to_hsv

    colors = {'timepoint': {'Day 70': '#FFDC49'}}
    sets = _sets('Day 70 o60', 'Day 70 u60', level='Day 70', level_column='timepoint',
                 tags=['o60', 'u60'])

    full, fragment = plotsVenn.venn_palette(sets, None, colors=colors)
    full_value = rgb_to_hsv(np.array(to_rgba(full)[:3]))[2]
    fragment_value = rgb_to_hsv(np.array(to_rgba(fragment)[:3]))[2]

    assert fragment_value < full_value, 'the derived tier must darken'


def test_the_derivation_matches_the_agreement_volcano():
    """One relationship, one encoding. A reader who has learnt the agreement volcano's tiers
    should not have to learn a second rule for the Venn, which is why both call through to
    plotsPalette.related_ramp rather than each scaling in its own way."""
    from trnagraph.modules import plotsAgreement

    colors = {'timepoint': {'Day 0': '#007FFF'}}
    sets = _sets('Day 0 o60', 'Day 0 u60', level='Day 0', level_column='timepoint',
                 tags=['o60', 'u60'])

    venn = plotsVenn.venn_palette(sets, None, colors=colors)
    agreement = plotsAgreement.direction_ramp('Day 0', 2, colormap={'Day 0': '#007FFF'})

    assert [to_rgba(c)[:3] for c in venn] == [to_rgba(c)[:3] for c in reversed(agreement)]


def test_an_explicit_venn_entry_still_beats_the_derivation():
    """`colors.venn` is the most specific source and stays authoritative, so a diagram that was
    hand-coloured before this change renders exactly as it did."""
    colors = {'timepoint': {'Day 0': '#007FFF'}}
    sets = _sets('Day 0 o60', 'Day 0 u60', level='Day 0', level_column='timepoint',
                 tags=['o60', 'u60'])

    palette = plotsVenn.venn_palette(sets, {'Day 0 u60': '#123456'}, colors=colors)

    assert to_rgba(palette[1])[:3] == to_rgba('#123456')[:3]


def test_a_level_with_no_named_colour_falls_back_to_the_palette():
    """Nothing to derive from is not an error -- the positional palette still colours the
    diagram, which is what an object with no style file gets."""
    sets = _sets('Day 0 o60', 'Day 0 u60', level='Day 0', level_column='timepoint',
                 tags=['o60', 'u60'])

    assert plotsVenn.venn_palette(sets, None, colors={}) == plotsPalette.venn_colors(2)


def test_the_two_automatic_diagrams_are_untouched():
    """fragment_vs_full_length and fiveprime_vs_threeprime carry no level at all, so there is
    no base to derive from and they keep the colourblind-safe positional pair."""
    sets = _sets('Fragment (u60)', 'Full length (o60)', tags=['u60', 'o60'])

    assert plotsVenn.venn_palette(sets, None, colors={'timepoint': {'Day 0': '#007FFF'}}) \
        == plotsPalette.venn_colors(2)


# --- alpha ---------------------------------------------------------------------------------

def test_venn_alpha_is_overridable_from_the_style_file():
    """Exposed rather than retuned: alpha trades the readability of an overlap against the
    readability of a single circle, and which side wins depends on how many circles there are.
    """
    fig, ax = plt.subplots()
    try:
        plotsVenn.draw_ellipse_venn(ax, SETS4, alpha=0.6)
        alphas = {round(p.get_facecolor()[3], 3) for p in ax.patches
                  if p.get_facecolor()[3] > 0}
    finally:
        plt.close(fig)

    assert alphas == {0.6}


def test_venn_alpha_defaults_to_the_tuned_value():
    fig, ax = plt.subplots()
    try:
        plotsVenn.draw_ellipse_venn(ax, SETS4)
        alphas = {round(p.get_facecolor()[3], 3) for p in ax.patches
                  if p.get_facecolor()[3] > 0}
    finally:
        plt.close(fig)

    assert alphas == {round(plotsVenn.VENN_ALPHA, 3)}
