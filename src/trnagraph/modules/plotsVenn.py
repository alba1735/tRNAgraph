#!/usr/bin/env python3
'''
Set-overlap diagrams over feature populations.

Two Venns are drawn automatically whenever `graph -g venn` runs and the data supports them:
fragment vs full-length (the u<N>/o<N> read-length split) and 5' vs 3' (the two end-specific
read types). Both answer the same shape of question -- is a given tRNA present as one species,
the other, or both -- which is why they are worth drawing for any dataset rather than being
declared per experiment.

The plot FAMILY is still gated behind a --config block: these are deliberate analyses, and a
user should not have them produced automatically and read conclusions off them.
'''

import logging
import os

import matplotlib.pyplot as plt
import numpy as np

from . import plotsPalette
from . import toolsTG

from .toolsSchemas import VennPlan, VennSet

plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

logger = logging.getLogger(__name__)


class InvalidVennError(ValueError):
    '''A Venn that cannot be drawn as specified -- wrong number of sets.'''

#: The BARE read type the fragment-vs-full-length Venn is measured on. Total rather than an
#: end-specific count: the split already partitions by length, so the question is whether the
#: tRNA is seen at all within each length class.
#:
#: Bare, and resolved through toolsTG.resolve_readtype() against the run's read basis, so these
#: diagrams follow --allreads like every other graph type. Hardcoding a basis here would let a
#: Venn rest on a different denominator than the volcano beside it with nothing saying so.
SPLIT_READTYPE = 'total'

#: The two bare end-specific read types the second automatic Venn contrasts.
END_READTYPES = ('fiveprime', 'threeprime')


#: The largest diagram that can be drawn HONESTLY. Four and five sets have known-good ellipse
#: arrangements; six does not -- the reference implementations switch to triangles there, and a
#: six-ellipse figure cannot represent all 63 regions, so counts would be missing from the
#: picture with nothing saying so.
MAX_VENN_SETS = 5

#: Opacity of a Venn circle. The circles are drawn semi-transparent so overlaps blend into a
#: third shade -- that blending IS the diagram, which is why this is a tuned default rather
#: than 1.0. Overridable from --style's `venn.alpha`, because it trades the readability of an
#: overlap against the readability of a single circle and the right side of that trade depends
#: on how many circles a diagram carries. matplotlib-venn's own default is 0.4.
VENN_ALPHA = 0.35
#: matplotlib-venn's area-proportional layouts.
PROPORTIONAL_MAX_SETS = 3


def venn_layout(n_sets):
    '''
    Which layout a diagram of `n_sets` circles gets: 'proportional' or 'ellipse'.

    Two and three are drawn area-proportionally, where the overlap is legible without reading a
    number. Four to six cannot be laid out proportionally with circles at all and use fixed
    ellipses, where the areas mean nothing and every region carries its count.

    Beyond five, returns 'upset_only'. Emitting a figure that cannot show every region is worse
    than not emitting one, and the UpSet companion drawn beside each complex diagram is exact at
    any size, so nothing is lost but the shape.
    '''
    if n_sets < 2:
        raise InvalidVennError(
            f'A Venn diagram needs at least 2 sets to show an overlap, got {n_sets}. '
            f'One set is a count, not a comparison.')
    if n_sets > MAX_VENN_SETS:
        # Not an error: the analysis is fine, only the Venn SHAPE cannot carry it. The UpSet
        # companion written beside every complex diagram is exact at any number of sets, so the
        # question still gets answered -- in the representation that can answer it.
        return 'upset_only'
    return 'proportional' if n_sets <= PROPORTIONAL_MAX_SETS else 'ellipse'


#: Ellipse arrangements for the diagram sizes matplotlib-venn cannot draw, as
#: (centre x, centre y, width, height, rotation degrees) in axes coordinates. These are the
#: classical 4- and 5-set Venn constructions -- every one of the 2**n - 1 regions exists as a
#: distinct area, which is the property that makes them usable and that no arrangement of six
#: ellipses has.
ELLIPSE_LAYOUTS = {
    4: [(0.350, 0.400, 0.72, 0.45, 140.0),
        (0.450, 0.500, 0.72, 0.45, 140.0),
        (0.544, 0.500, 0.72, 0.45, 40.0),
        (0.644, 0.400, 0.72, 0.45, 40.0)],
    5: [(0.428, 0.449, 0.87, 0.50, 155.0),
        (0.469, 0.543, 0.87, 0.50, 82.0),
        (0.558, 0.523, 0.87, 0.50, 10.0),
        (0.578, 0.432, 0.87, 0.50, 118.0),
        (0.489, 0.383, 0.87, 0.50, 46.0)],
}

#: Samples per axis when locating region centroids. 600 puts roughly 30 points in the smallest
#: region of the 5-set arrangement, which is enough for a stable centroid without being slow.
_PLACEMENT_GRID = 600


def _membership_masks(layout, grid=_PLACEMENT_GRID):
    '''
    A grid over the unit square, labelled with which ellipses contain each point.

    Returns (xs, ys, masks) where `masks` is an integer bitmask per point: bit i set when the
    point lies inside ellipse i. Regions are then found by grouping points by mask, which is
    what lets label placement follow the ACTUAL geometry rather than a table of coordinates that
    could drift out of step with it.
    '''
    axis = np.linspace(0.0, 1.0, grid)
    xs, ys = np.meshgrid(axis, axis)
    masks = np.zeros(xs.shape, dtype=np.int32)
    for i, (cx, cy, width, height, angle) in enumerate(layout):
        theta = np.radians(angle)
        dx, dy = xs - cx, ys - cy
        # Rotate into the ellipse's own frame, then apply the unit-circle test.
        rx = dx * np.cos(theta) + dy * np.sin(theta)
        ry = -dx * np.sin(theta) + dy * np.cos(theta)
        inside = (rx / (width / 2)) ** 2 + (ry / (height / 2)) ** 2 <= 1.0
        masks |= inside.astype(np.int32) << i
    return xs, ys, masks


#: How far beyond its ellipse's outward extremity a set label sits, in axes coordinates.
_LABEL_MARGIN = 0.04


def _label_anchor(ellipse, centre=(0.5, 0.5), samples=720):
    '''
    Where one ellipse's set label goes: just past the boundary point farthest from the middle.

    A fixed offset from each ellipse's CENTRE collides whenever two ellipses share a centre
    coordinate, which the 4-set arrangement does by construction -- two pairs sit at the same
    height, so their labels printed on top of one another. The outward extremity is distinct for
    every ellipse in both arrangements, because that is what a Venn layout is: one shape rotated
    about a common middle. It also keeps the label clear of the region counts, which live between
    the ellipses rather than at their tips.
    '''
    cx, cy, width, height, angle = ellipse
    theta = np.radians(angle)
    t = np.linspace(0.0, 2 * np.pi, samples, endpoint=False)
    # The ellipse boundary, parameterised and rotated into the diagram's frame.
    xs = cx + (width / 2) * np.cos(t) * np.cos(theta) - (height / 2) * np.sin(t) * np.sin(theta)
    ys = cy + (width / 2) * np.cos(t) * np.sin(theta) + (height / 2) * np.sin(t) * np.cos(theta)
    far = int(np.argmax((xs - centre[0]) ** 2 + (ys - centre[1]) ** 2))
    outward = np.array([xs[far] - centre[0], ys[far] - centre[1]])
    outward /= np.linalg.norm(outward)
    return (float(xs[far] + _LABEL_MARGIN * outward[0]),
            float(ys[far] + _LABEL_MARGIN * outward[1]))


def draw_ellipse_venn(ax, sets, colors=None, alpha=None):
    '''
    Draw a 4- or 5-set Venn with fixed ellipses, returning (placed, unplaceable).

    `placed` maps each populated region to (count, x, y); `unplaceable` lists any populated
    region the geometry cannot show. That second return value is the point of the whole approach:
    a region with no representable area is REPORTED, never silently absent from a figure the
    reader will take as complete.

    Areas carry no meaning at this size -- unlike the 2- and 3-set proportional layouts -- so
    every region is labelled with its count.
    '''
    from matplotlib.patches import Ellipse

    labels = list(sets)
    layout = ELLIPSE_LAYOUTS.get(len(labels))
    if layout is None:
        raise InvalidVennError(
            f'No ellipse arrangement for {len(labels)} sets; {sorted(ELLIPSE_LAYOUTS)} are '
            f'available. venn_layout() decides which sizes reach here.')

    palette = colors or plotsPalette.venn_colors(len(labels))
    for (cx, cy, width, height, angle), colour in zip(layout, palette):
        ax.add_patch(Ellipse((cx, cy), width, height, angle=angle, facecolor=colour,
                             edgecolor='none', alpha=VENN_ALPHA if alpha is None else alpha))
        ax.add_patch(Ellipse((cx, cy), width, height, angle=angle, facecolor='none',
                             edgecolor=plotsPalette.REFERENCE_LINE, linewidth=0.8))

    xs, ys, masks = _membership_masks(layout)
    regions = exclusive_regions(sets)
    placed, unplaceable = {}, []
    for region, members in regions.items():
        wanted = 0
        for i, label in enumerate(labels):
            if label in region.split(' & '):
                wanted |= 1 << i
        selection = masks == wanted
        if not selection.any():
            unplaceable.append(region)
            continue
        x, y = float(xs[selection].mean()), float(ys[selection].mean())
        placed[region] = (len(members), x, y)
        ax.text(x, y, str(len(members)), ha='center', va='center', fontsize=8)

    anchors = [_label_anchor(ellipse) for ellipse in layout]
    for (x, y), label in zip(anchors, labels):
        ax.text(x, y, label, ha='center', va='center', fontsize=9)

    # Square limits, sized to the labels rather than to the ellipses, so the aspect stays equal
    # (the arrangement is only a Venn while it is undistorted) and no label sits on the frame.
    span = max(abs(v - 0.5) for anchor in anchors for v in anchor) + 0.08
    ax.set_xlim(0.5 - span, 0.5 + span)
    ax.set_ylim(0.5 - span, 0.5 + span)
    ax.set_aspect('equal')
    ax.axis('off')
    if unplaceable:
        logger.warning(f'{len(unplaceable)} populated region(s) have no representable area in a '
                       f'{len(labels)}-set diagram and are absent from the figure: {unplaceable}. '
                       f'The UpSet plot beside it carries every region.')
    return placed, unplaceable


def upset_intersection_sizes(sets):
    '''
    {region label: size} as an UpSet plot would report them.

    Derived from exclusive_regions() rather than computed independently, which is the whole
    point: the UpSet plot and the Venn beside it are two renderings of ONE computation, so they
    cannot drift. Independent derivations would be two chances to be wrong and no way to notice.
    '''
    return {region: len(members) for region, members in exclusive_regions(sets).items()}


#: How "main" a read-length variant is, deciding which circle keeps the level's own colour.
#: The unsplit population ranks above full-length, which ranks above fragment -- so a
#: Day 0 x {o60, u60} pair draws full-length in Day 0's own blue and the fragment in a
#: darker, rotated shade of it, rather than two blues a reader has to squint at.
def _variant_rank(tag):
    if tag == 'full':
        return 2
    return 1 if str(tag)[:1] == 'o' else 0


def _derived_level_colors(plan_sets, colors):
    '''
    Colour per set label for circles that are one grouping level seen through several variants.

    Returns {} unless a level appears more than once AND --style names a colour for it, which
    is the only case where there is a base to derive from and something to tell apart. Sets
    with no level -- the two automatic diagrams -- are left to the positional palette.

    The derivation is plotsPalette.related_ramp, the same one the agreement volcano's tiers
    use: darken and rotate hue off the level's own colour rather than lighten it. Lightening
    is what a style file had to do by hand, and it washes out twice on a Venn, once in the
    tint and again in the alpha the circles are blended with.
    '''
    if not colors:
        return {}
    by_level = {}
    for spec in plan_sets:
        if spec.level is None or spec.level_column is None:
            continue
        by_level.setdefault((spec.level_column, spec.level), []).append(spec)

    derived = {}
    for (column, level), specs in by_level.items():
        if len(specs) < 2:
            continue
        base = (colors.get(column) or {}).get(level)
        if base is None:
            continue
        # Strongest tier LAST out of related_ramp, and strongest goes to the most "main"
        # variant, so the colour the style file actually named still appears on the figure.
        ordered = sorted(specs, key=lambda spec: _variant_rank(spec.tag))
        ramp = plotsPalette.related_ramp(base, len(ordered))
        for spec, colour in zip(ordered, ramp):
            derived[spec.label] = colour
    return derived


def venn_palette(plan_sets, colormap=None, colors=None):
    '''
    The colour for each set in `plan_sets`, in order.

    Three sources, most specific first. An explicit `colors.venn` entry keyed by SET LABEL
    always wins -- keyed by label rather than by grouping level because a diagram can carry the
    same level more than once (the four-circle case crosses two timepoints with two read-length
    variants, so "Day 0" names two different circles), and keying on the level would hand both
    the same colour and destroy the diagram.

    Failing that, a level seen through several variants derives a ramp from that level's own
    entry elsewhere in `colors` -- so naming `colors.timepoint` once is enough, and the four
    hand-tinted `colors.venn` hexes a crossed diagram used to need are no longer required.

    Failing both, the positional default palette. A partial mapping is a SUBSTITUTION into it,
    not a replacement: naming one circle leaves the rest where they were, which means a partial
    override can repeat a colour the file also named explicitly -- predictable beats clever
    here, since the alternative is a palette that shifts under you when you add an entry.
    '''
    plan_sets = list(plan_sets)
    labels = [spec.label for spec in plan_sets]
    defaults = plotsPalette.venn_colors(len(labels))
    derived = _derived_level_colors(plan_sets, colors)
    return [(colormap or {}).get(label) or derived.get(label) or default
            for label, default in zip(labels, defaults)]


def draw_upset(sets, title=None, settings=None, colors=None):
    '''
    The UpSet companion for one diagram, returned as a figure.

    Drawn beside every complex Venn. At four and five sets a Venn is drawable but hard to read --
    locating a particular combination means tracing overlapping ellipses -- while UpSet makes
    each intersection a labelled row. Beyond five it is the only honest option, since no ellipse
    arrangement can show every region.

    UpSet is Lex et al. 2014 (IEEE TVCG 20(12):1983-1992); using the cited implementation rather
    than drawing something bespoke is deliberate.
    '''
    from upsetplot import UpSet, from_contents

    contents = {label: set(members) for label, members in sets.items()}
    if not any(contents.values()):
        raise InvalidVennError('Every set is empty; there is no intersection to plot.')
    data = from_contents(contents)
    fig = plt.figure(figsize=toolsTG.figsize_for(settings, (9, 5)))
    upset = UpSet(data, show_counts=True, sort_by='cardinality')

    # One colour per set, shared with the Venn beside it: reading the pair together means
    # matching a row to a circle, and two palettes would make that a lookup rather than a glance.
    by_label = dict(zip(sets, colors or plotsPalette.venn_colors(len(sets))))
    for label, colour in by_label.items():
        upset.style_categories(label, bar_facecolor=colour)

    axes = upset.plot(fig=fig)
    _color_upset_dots(axes, by_label)
    if title:
        fig.suptitle(title, fontsize=10)
    # plot()'s axes dict is the only handle on the drawn matrix, and it is not attached anywhere.
    fig._trnagraph_upset_axes = axes
    return fig


def _color_upset_dots(axes, by_label):
    '''
    Tint each filled matrix dot with its row's set colour.

    upsetplot draws every dot in one PathCollection and separates present from absent by ALPHA
    (1.0 against 0.18) rather than by colour, so the filled ones can be recoloured while absent
    ones stay grey -- an absent dot means "not in this set", and giving it the set's colour would
    say the opposite. Rows are matched by tick label, not by input order, because upsetplot sorts
    categories itself.
    '''
    from matplotlib.colors import to_rgba

    matrix = axes['matrix']
    labels = [text.get_text() for text in matrix.get_yticklabels()]
    if not matrix.collections:
        return
    dots = matrix.collections[0]
    faces = np.asarray(dots.get_facecolors()).copy()
    for i, (_, row) in enumerate(np.asarray(dots.get_offsets())):
        if faces[i][3] < 0.9:
            continue
        label = labels[int(round(row))]
        if label in by_label:
            faces[i] = to_rgba(by_label[label])
    dots.set_facecolors(faces)


def draw_venn(ax, sets, colors=None, alpha=None):
    '''
    Draw a Venn onto `ax`, dispatching on set count (see venn_layout).

    Region counts come from exclusive_regions(), the same computation the TSV is written from,
    so the picture and the table beside it cannot report different numbers. matplotlib-venn is
    used here rather than the fixed-ellipse drawing needed at 4+ sets because at two and three
    it lays the circles out AREA-PROPORTIONALLY: the overlap is legible without reading a count.
    '''
    from matplotlib_venn import venn2, venn3

    labels = list(sets)
    layout = venn_layout(len(labels))
    if layout == 'upset_only':
        raise InvalidVennError(
            f'{len(labels)} sets cannot be drawn as a Venn; use draw_upset(). '
            f'venn_layout() decides this.')
    if layout == 'ellipse':
        return draw_ellipse_venn(ax, sets, colors=colors, alpha=alpha)
    members = [set(sets[label]) for label in labels]
    draw = venn2 if len(labels) == 2 else venn3
    palette = colors or plotsPalette.venn_colors(len(labels))
    return draw(members, set_labels=labels, ax=ax, set_colors=tuple(palette),
                alpha=VENN_ALPHA if alpha is None else alpha)


def exclusive_regions(sets):
    '''
    {region label: [feature, ...]} where every feature appears in exactly ONE region.

    A region is named by the full combination of sets containing the feature, joined with " & ",
    so a 2-set diagram yields "A", "B" and "A & B". This is the same membership an UpSet plot
    encodes -- which is what will let a complex Venn and its UpSet companion report identical
    numbers rather than two computations that merely ought to agree.

    Empty regions are omitted: a row reading `A & B  0` says nothing a reader needs, and with N
    sets there are 2**N - 1 possible regions, most of them empty on real data.
    '''
    labels = list(sets)
    members = {label: set(features) for label, features in sets.items()}
    regions = {}
    for feature in sorted(set().union(*members.values()) if members else set()):
        containing = tuple(label for label in labels if feature in members[label])
        if containing:
            regions.setdefault(' & '.join(containing), []).append(feature)
    return regions


def write_membership_table(path, sets, provenance):
    '''
    Write one Venn's membership to a TSV: which features fall in each exclusive region.

    The figure shows the counts; this is what a reader can act on. A commented provenance header
    names the object and every parameter behind the numbers, so a file left behind in results/
    after a rebuild identifies itself rather than quietly disagreeing with the object. results/
    is a convenience for people and papers -- the AnnData object remains the source of truth and
    nothing reads this back.
    '''
    regions = exclusive_regions(sets)
    lines = ['# tRNAgraph Venn membership']
    lines += [f'# {key}: {value}' for key, value in provenance.items()]
    lines.append('\t'.join(['rank', 'region', 'n', 'features']))
    # Largest region first, so the file opens on what most features share rather than on
    # whichever region the set order happened to enumerate first. Ties keep enumeration order,
    # which is the order the diagram's own labels run in.
    ordered = sorted(regions.items(), key=lambda item: -len(item[1]))
    for rank, (region, features) in enumerate(ordered, start=1):
        lines.append('\t'.join([str(rank), region, str(len(features)), ','.join(features)]))
    path = str(path)
    with open(path, 'w') as handle:
        handle.write('\n'.join(lines) + '\n')
    return path


def require_multivariate_config(config):
    '''
    Return the `multivariate` config block, or refuse the run naming what is missing.

    These plots are opt-in by design (see adataGraph.GRAPH_TYPES_ALL). Refusing here, by name,
    is the difference between a user learning they need to declare the analysis and a user
    getting an empty output directory and no explanation.
    '''
    block = getattr(config, 'multivariate', None) if config is not None else None
    if block is None:
        raise InvalidVennError(
            "`-g venn` needs a `multivariate` block in your --config file, which declares the "
            "grouping column and thresholds the sets are built from. These analyses are "
            "deliberately opt-in rather than part of `-g all`: the sets and cutoffs are choices "
            "about your experiment, and figures produced from unchosen ones invite wrong "
            "conclusions. Minimal example: "
            '{"name": "my_analysis", "multivariate": {"grouping": "group"}}')
    return block


def _split_pairs(adata):
    '''Matched (under, over) split tags, keyed by cutoff -- u60 with o60, not with o50.'''
    tags = list(adata.uns.get('size_splits', {}))
    cutoffs = {}
    for tag in tags:
        if tag[:1] in ('u', 'o') and tag[1:].isdigit():
            cutoffs.setdefault(tag[1:], {})[tag[:1]] = tag
    return {cutoff: (pair['u'], pair['o'])
            for cutoff, pair in sorted(cutoffs.items()) if {'u', 'o'} <= set(pair)}


def declared_venn_plans(adata, block, read_basis=None, variant_tag='full'):
    '''
    Resolve `multivariate.venn` declarations into drawable plans.

    Turns each declared circle into a VennSet: the variant string becomes a tag (validated by
    parse_variant, so an unknown split fails with its message rather than a KeyError later), the
    bare read type becomes an obs column for the run's basis, and the grouping column is attached
    so the set can filter on its own level.

    Labels are generated from whatever actually DISTINGUISHES the circles, so a four-way diagram
    does not end up with four labels reading "tRNAs". An explicit `label` always wins.
    '''
    basis = read_basis or toolsTG.READ_BASIS_UNIQUE
    plans, overridden = [], []
    for declaration in (block.venn or []):
        sets = []
        for spec in declaration.sets:
            # A declared variant wins; --variant only fills the gaps. A diagram contrasting u60
            # with o60 cannot be forced onto one variant, which is exactly what it is for.
            if spec.variant:
                tag = toolsTG.parse_variant(adata, spec.variant).tag
                if tag != variant_tag:
                    overridden.append(spec.variant)
            else:
                tag = variant_tag
            readtype = toolsTG.resolve_readtype(spec.readtype or SPLIT_READTYPE, basis, adata)
            sets.append(VennSet(
                label=spec.label or _distinguishing_label(spec, declaration.sets),
                readtype=readtype, tag=tag, level=spec.level,
                level_column=block.grouping if spec.level is not None else None))
        plans.append(VennPlan(
            name=declaration.name,
            title=declaration.title or declaration.name.replace('_', ' '),
            sets=sets))
    if overridden:
        # Stated rather than silent: a flag that was disregarded is worse unmentioned than
        # refused, since the user has no other way to learn the figure ignored it.
        logger.info(f'--variant {variant_tag!r} is ignored for Venn sets that name their own '
                    f'variant ({sorted(set(overridden))}); each set uses the one it declares.')
    return plans


def _distinguishing_label(spec, siblings):
    '''
    A label naming only the fields that VARY across a diagram's circles.

    Including a field every circle shares adds noise ("D0 norm:u60 total" x4 differing in two
    words); including none leaves them indistinguishable. So the label is built from exactly the
    fields that differ.
    '''
    parts = []
    for field in ('level', 'variant', 'readtype'):
        values = {getattr(sibling, field) for sibling in siblings}
        value = getattr(spec, field)
        if len(values) > 1 and value is not None:
            parts.append(str(value))
    return ' '.join(parts) or 'all features'


def simple_venn_plans(adata, read_basis=None, readtype=SPLIT_READTYPE):
    '''
    The automatic Venns this object can support, and a message for each it cannot.

    Returns (plans, skipped). A missing prerequisite is never an error: a plain build has no
    read-length split and a build without end-specific counts has no 5'/3' contrast, and in both
    cases the other Venn is still worth drawing. But the skip is REPORTED by name, because a
    figure that silently fails to appear is easily mistaken for a biological absence.
    '''
    plans, skipped = [], []
    basis = read_basis or toolsTG.READ_BASIS_UNIQUE
    readtype = toolsTG.resolve_readtype(readtype, basis, adata)
    end_readtypes = [toolsTG.resolve_readtype(rt, basis, adata) for rt in END_READTYPES]

    pairs = _split_pairs(adata)
    if not pairs:
        skipped.append(
            "Skipping the fragment_vs_full_length Venn: this object has no read-length split "
            "variant. Add one with `trnagraph analyze addsplit -c <cutoff>`, or build with "
            "--readlengthsplit.")
    elif readtype not in adata.obs.columns:
        skipped.append(f"Skipping the fragment_vs_full_length Venn: obs has no '{readtype}'.")
    else:
        for cutoff, (under, over) in pairs.items():
            name = ('fragment_vs_full_length' if len(pairs) == 1
                    else f'fragment_vs_full_length_{cutoff}')
            plans.append(VennPlan(
                name=name,
                title=f'Fragment (<{cutoff}nt) vs full-length (>={cutoff}nt) tRNAs',
                sets=[VennSet(label=f'Fragment ({under})', readtype=readtype, tag=under),
                      VennSet(label=f'Full length ({over})', readtype=readtype, tag=over)]))

    missing = [rt for rt in end_readtypes if rt not in adata.obs.columns]
    if missing:
        skipped.append(
            f"Skipping the fiveprime_vs_threeprime Venn: obs has no {missing}. End-specific "
            f"counts are written by `analyze build`; an object built before they existed, or "
            f"one filtered down to other read types, will not carry them.")
    else:
        plans.append(VennPlan(
            name='fiveprime_vs_threeprime',
            title="5' vs 3' tRNA fragments",
            sets=[VennSet(label="5' counts", readtype=end_readtypes[0]),
                  VennSet(label="3' counts", readtype=end_readtypes[1])]))

    return plans, skipped


if __name__ == '__main__':
    pass


def _set_members(adata, spec, cutoff):
    '''
    The features present in one circle.

    Presence is pooled across every sample, not split by group: these two automatic Venns ask
    whether a tRNA is seen at all as a fragment / as a 5' end, which is a property of the
    dataset rather than of one condition. A synthetic single-level grouping reuses
    presence_sets() rather than duplicating its cutoff rule.
    '''
    obs = toolsTG.variant_obs(adata, spec.tag)
    if spec.level is not None:
        if not spec.level_column:
            raise InvalidVennError(
                f"Venn set {spec.label!r} names level {spec.level!r} but no column to look it up "
                f"in. Set `multivariate.grouping` to the obs column the levels belong to.")
        obs = obs[obs[spec.level_column].astype(str) == str(spec.level)]
    obs = obs.assign(_pooled='all')
    return toolsTG.presence_sets(obs, '_pooled', spec.readtype, cutoff)['all']


def visualizer(adata, block, output, config_name='default', settings=None,
               read_basis=None, variant_tag='full', threaded=False, colormap=None,
               colors=None, cutoff=20, results_dir=None):
    '''
    Draw every Venn this object supports, store the membership, and write the tables.

    Order matters: membership is computed once, stored on the object, and BOTH the figure and
    the TSV are rendered from that one computation, so the three can never disagree.
    '''
    individual_output = f'{output}individual/'
    messages = []
    plans, skipped = simple_venn_plans(adata, read_basis=read_basis)
    plans += declared_venn_plans(adata, block, read_basis=read_basis, variant_tag=variant_tag)
    for message in skipped:
        messages.append(message)
        logger.warning(message)
    if not plans:
        return '\n'.join(messages) + '\n' if threaded else None

    if results_dir:
        messages.append(toolsTG.builder(results_dir))

    messages.append(toolsTG.builder(individual_output))
    provenance_base = {'config': config_name, 'cutoff': cutoff,
                       'membership': 'expression-presence'}

    for plan in plans:
        sets = {spec.label: _set_members(adata, spec, cutoff)
                for spec in plan.sets}
        provenance = dict(provenance_base, plan=plan.name,
                          readtypes=sorted({s.readtype for s in plan.sets}),
                          tags=sorted({s.tag for s in plan.sets}))
        # Keyed by the DIAGRAM, not by a variant tag: a Venn spans variants, so every one
        # in a run would otherwise share the same key and overwrite its predecessor.
        toolsTG.write_multivariate(adata, config_name, 'venn', plan.name, sets, provenance)

        subtitle = f'present at mean >= {cutoff:g} normalized reads'
        # Resolved ONCE per diagram and handed to both renderings: the pairing only works while
        # a row and its circle are the same colour.
        palette = venn_palette(plan.sets, colormap, colors=colors)
        layout = venn_layout(len(sets))
        if layout == 'upset_only':
            messages.append(f'{plan.name}: {len(sets)} sets cannot be drawn as a Venn; '
                            f'writing the UpSet plot only.')
        else:
            fig, ax = plt.subplots(figsize=toolsTG.figsize_for(settings, (6, 6)))
            draw_venn(ax, sets, colors=palette, alpha=(settings or {}).get('alpha'))
            ax.set_title(f'{plan.title}\n({subtitle})', fontsize=10)
            path = f'{individual_output}{plan.name}_venn.pdf'
            fig.savefig(path, bbox_inches='tight')
            plt.close(fig)
            messages.append(f'Saving Venn: {path}')

        # Every diagram past the proportional sizes gets an UpSet companion: at four and five
        # sets the Venn is drawable but hard to read, and beyond that it is the only honest
        # representation. Both are rendered from one computation, so they cannot disagree.
        if layout != 'proportional':
            upset_fig = draw_upset(sets, title=f'{plan.title} ({subtitle})', settings=settings,
                                   colors=palette)
            upset_path = f'{individual_output}{plan.name}_upset.pdf'
            upset_fig.savefig(upset_path, bbox_inches='tight')
            plt.close(upset_fig)
            messages.append(f'Saving UpSet: {upset_path}')

        if results_dir:
            table = os.path.join(results_dir, f'{plan.name}_venn.tsv')
            write_membership_table(table, sets, dict(provenance, object=adata.uns.get(
                'trnagraphruninfo', {}).get('expname', 'unknown')))
            messages.append(f'Saving Venn membership table: {table}')

    for message in messages:
        if not threaded:
            logger.info(message)
    return '\n'.join(messages) + '\n' if threaded else None
