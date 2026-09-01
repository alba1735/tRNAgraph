#!/usr/bin/env python3
'''
Every color tRNAgraph draws with, in one place.

Colors used to be literals scattered across the plots*.py modules, in several unrelated
conventions: volcano hardcoded hex, heatmap built a diverging palette inline, seqlogo carried
its own nucleotide dict, and eight modules each called `sns.husl_palette(n)` (or
`sns.color_palette('husl', n)`, the same thing spelled differently) as their categorical
fallback. Making the look uniform meant finding and matching all of them by hand.

This module holds the values only -- no plotting logic and no matplotlib import at module
scope, so it stays cheap to import and safe to read from anywhere.

The defaults have since been retuned for colorblind safety and perceptual uniformity. Two sets
were deliberately NOT changed, because they are inherited conventions rather than choices this
project gets to make: the tRNA arm/loop bands come from gtRNAdb's own palette, and the
nucleotide colors follow the sequence-logo convention that Logomaker calls `classic` (G orange,
T/U red, C blue, A green), which readers decode by memory. Both contain red/green pairs that an
accessibility pass will flag; leave them alone. A --style file can override every ramp and the
categorical fallback per run, so a user who needs different values is not blocked on these.

Naming is by ROLE rather than by appearance -- `REFERENCE_LINE` not `BLACK` -- so that a later
pass can change what a role looks like without every call site reading as a lie.
'''

# --- Categorical fallback ------------------------------------------------------------
# Used whenever a plot needs N visually distinct colors for unordered categories and the user
# supplied no palette for that column.
#
# Okabe-Ito extended to ten, the same values as seaborn's 'colorblind', written out rather than
# named so the default is legible here (this module being the documented control point) and so
# it cannot shift if seaborn retunes that palette. Ten is not arbitrary: the columns tRNAgraph
# actually colors by are `group` (3), `sample` (9-10), `amino` (22-24) and `iso`/anticodon
# (33-54), so ten spans every one of them except the two that no qualitative palette can serve.
CATEGORICAL_COLORS = (
    '#0173b2', '#de8f05', '#029e73', '#d55e00', '#cc78bc',
    '#ca9161', '#fbafe4', '#949494', '#ece133', '#56b4e9',
)

# Above CATEGORICAL_COLORS' length no colorblind-safe qualitative set exists, so the fallback
# switches to husl. This is a trade, not a win: husl is not colorblind-safe, but cycling ten
# colors across 22 amino acids would give four groups the same color and make the plot actively
# misleading, whereas husl at least keeps them nominally distinct. The real answer at that count
# is to stop encoding the dimension in color at all -- a plot-design change, tracked separately.
CATEGORICAL_OVERFLOW_PALETTE = 'husl'


def categorical_palette(n):
    '''N distinct colors for unordered categories -- the shared fallback when --style sets none.'''
    import logging

    import seaborn as sns

    n = int(n)
    if n <= len(CATEGORICAL_COLORS):
        return list(CATEGORICAL_COLORS[:n])
    logging.getLogger(__name__).warning(
        f'{n} categories exceeds the {len(CATEGORICAL_COLORS)}-color colorblind-safe palette; '
        f'falling back to {CATEGORICAL_OVERFLOW_PALETTE}, which is not colorblind-safe at this '
        f'count. Consider grouping the categories or distinguishing them by something other '
        f'than color.'
    )
    return sns.color_palette(CATEGORICAL_OVERFLOW_PALETTE, n)


# --- Continuous / ordered scales -----------------------------------------------------
# Named by what they encode, so a future pass can retune one without disturbing the others.
# Only two pairs of these are ever drawn in one figure, and those are the only pairs that must
# be mutually distinguishable: `lfc` beside `significance` (heatmap's two panels), and `score`
# beside `sequence` (seqlogo). `correlation` and `ordered` each appear alone.
SEQUENTIAL_CORRELATION = 'crest'    # correlation R^2 heatmap
SEQUENTIAL_SIGNIFICANCE = 'rocket_r'  # -log10(p) heatmap panel; warm, against lfc's cool diverging
SEQUENTIAL_SCORE = 'mako_r'         # seqlogo per-position score heatmap and its colorbar
SEQUENTIAL_SEQUENCE = 'Greys'       # seqlogo sequence background: deliberately neutral, not data
# Ordered specificity scale (coverage partition) and cluster numeric coloring: light = least,
# dark = most, which is why it is a sequential ramp rather than categorical hues.
SEQUENTIAL_ORDERED = 'mako_r'

# Diverging scale for log2 fold change, centered on zero. Blue-white-red: the previous
# blue-to-green ramp put the two ends of the scale on a deuteranopia confusion pair, which made
# the sign of a fold change -- the whole point of the plot -- the hardest thing to read on it.
DIVERGING_LFC = 'vlag'


# --- Differential-expression direction -----------------------------------------------
# Volcano point colors when --style names no color for the compared groups.
DE_UP = '#d62728'
DE_DOWN = '#1f77b4'
DE_NONSIGNIFICANT = '#7f7f7f'

# --- Consensus mismatch emphasis -----------------------------------------------------
# The seqlogo's consensus row prints a base in a highlight color where the read consensus
# disagrees with the reference. This is a SEMANTIC emphasis color, unrelated to the nucleotide
# convention below -- it marks disagreement, not identity, so it is free to change.
#
# Two values because the text sits on the score heatmap, which runs light at low scores and
# dark at high ones; a single color cannot stay readable across both, which is why the plain
# consensus text already switches between black and white at the same threshold. Orange rather
# than red: red on a dark ramp is close to unreadable and is the worst hue to pair with the
# green in the nucleotide set for a colorblind reader.
CONSENSUS_MISMATCH_ON_LIGHT = '#b35900'
CONSENSUS_MISMATCH_ON_DARK = '#ffa600'

# --- Nucleotides ---------------------------------------------------------------------
# Seqlogo letter colors. U shares T's color: the same base, differing only by whether the
# plot is rendering RNA or DNA (--logornamode).
NUCLEOTIDE_COLORS = {
    'T': '#ea4335',
    'U': '#ea4335',
    'A': '#34a853',
    'C': '#4285f4',
    'G': '#fbbc05',
    '-': '#ffffff',
}

# --- tRNA structural annotation ------------------------------------------------------
# Arm/loop bands drawn behind the seqlogo, and the acceptor/anticodon stem shading.
ARM_DLOOP = '#ff9999'
ARM_ANTICODON = '#99ff99'
ARM_TLOOP = '#99ffff'
ARM_STEM = '#98c0c0'
# Laid over an arm band to lighten its interior, distinguishing loop from stem.
ARM_INNER_WASH = 'white'
ARM_INNER_WASH_ALPHA = 0.8

# --- Read-end emphasis ---------------------------------------------------------------
# Read starts and ends, drawn as semi-transparent bars over a coverage trace. They must stand
# out against the fill behind them AND against each other, but the previous magenta/cyan pair
# was loud at the cost of being unpleasant and low-contrast once alpha-blended onto white.
# ColorBrewer PRGn's endpoints instead: a colorblind-safe pair that stays legible at alpha, and
# that is deliberately outside CATEGORICAL_COLORS, since the trace behind them is colored from
# that set and the overlay must not be mistaken for another group.
READSTART_COLOR = "#2dc0ff"  # 5' end
READEND_COLOR = "#9bff38"    # 3' end

# --- Neutral furniture ---------------------------------------------------------------
# Non-data ink: reference lines, shaded regions, bar outlines, placeholder swatches.
REFERENCE_LINE = 'black'
BAR_EDGE = 'black'
GRID_LINE = 'lightgray'
AXIS_TEXT = 'black'
AXIS_TEXT_MUTED = 'dimgrey'
COVERAGE_TRACE_MUTED = 'dimgrey'
# Shaded background band marking a region of the tRNA (coverage plots).
REGION_SHADE = '#cacaca'
REGION_SHADE_ALPHA = 0.35
# The collapsed non-tRNA band in stacked count plots: background against the tRNA categories,
# which are the subject of that plot.
NONTRNA_COLLAPSED = '#b0b0b0'
# Marker fill used in legends to show a SHAPE rather than a color (tRNA circle vs non-tRNA
# square), so the swatch must not imply a category color.
LEGEND_SHAPE_SWATCH = 'grey'
BACKGROUND = 'white'


# --- Resolution: turning a --style file's palette blocks into drawable values --------
# Everything above is a DEFAULT. The helpers below let a --style file override the ordered
# scales and the categorical fallback without any plot module learning where the override
# came from: modules call gradient()/categorical() with the `settings` dict they are already
# given, and get back a drawable value whether or not the user configured anything.

# Role -> the colormap name this module defaults it to.
GRADIENT_ROLES = ('correlation', 'significance', 'score', 'sequence', 'ordered', 'lfc')

_GRADIENT_DEFAULTS = {
    'correlation': SEQUENTIAL_CORRELATION,
    'significance': SEQUENTIAL_SIGNIFICANCE,
    'score': SEQUENTIAL_SCORE,
    'sequence': SEQUENTIAL_SEQUENCE,
    'ordered': SEQUENTIAL_ORDERED,
    'lfc': DIVERGING_LFC,
}


def _default_gradient(role):
    '''
    The built-in colormap for one role, as a Colormap object.

    Goes through build_colormap() rather than plt.get_cmap() directly because SEQUENTIAL_ORDERED
    is `mako_r`, which seaborn registers on import -- resolving it without that import raises
    'not a valid value for name' on any path that has not already pulled seaborn in.
    '''
    return build_colormap(_GRADIENT_DEFAULTS[role])


def build_colormap(value):
    '''
    Turn one style-file gradient value into a Colormap.

    A string is a registered colormap name (seaborn is imported first, since `mako`/`vlag`
    and friends are registered by seaborn rather than shipped with matplotlib). A list is
    two or more color tokens interpolated into an evenly-spaced ramp -- which is what lets a
    user supply lab or journal colors without knowing matplotlib's colormap API.
    '''
    import seaborn as sns  # noqa: F401  -- registers mako/rocket/crest/flare/vlag/icefire
    import matplotlib.pyplot as plt
    from matplotlib.colors import Colormap, LinearSegmentedColormap

    if isinstance(value, Colormap):
        return value
    if isinstance(value, str):
        return plt.get_cmap(value)
    return LinearSegmentedColormap.from_list('trnagraph_custom', list(value))


def resolve_gradients(spec):
    '''
    All six roles as Colormap objects: the built-in default wherever `spec` sets nothing.

    Every graph type receives every role rather than only the ones it draws with. Mapping
    roles to graph types would have to be maintained by hand as plots are added, and a plot
    reading a role the map forgot to grant it would fail inside a worker process -- the
    failure class tests/unit/test_plot_settings_scope.py exists to catch.
    '''
    spec = spec or {}
    if hasattr(spec, 'model_dump'):
        spec = spec.model_dump(exclude_none=True)
    return {role: (build_colormap(spec[role]) if spec.get(role) is not None
                   else _default_gradient(role))
            for role in GRADIENT_ROLES}


def gradient(settings, role):
    '''
    The Colormap for one role, honouring --style when the caller was given resolved settings.

    Plot modules go through here rather than reading the SEQUENTIAL_*/DIVERGING_* constants
    directly, so a user's override reaches them without any module knowing about style files.
    Falls back to the built-in default when called with no settings, which keeps the modules
    usable from tests and from any path that has not resolved a style.
    '''
    if role not in GRADIENT_ROLES:
        raise KeyError(f'Unknown gradient role {role!r}. Known roles: {list(GRADIENT_ROLES)}')
    resolved = (settings or {}).get('gradients')
    if resolved and resolved.get(role) is not None:
        return resolved[role]
    return _default_gradient(role)


def discrete_colors(cmap, n):
    '''
    N colors sampled from a continuous ramp, matching `sns.color_palette(name, n)` exactly.

    seaborn samples a continuous colormap at `linspace(0, 1, n + 2)[1:-1]`, trimming both
    extremes; sampling at `linspace(0, 1, n)` instead shifts colors by as much as 0.4 per
    channel. Call sites that take discrete colors from a ramp use this so that switching them
    from a colormap NAME to a resolved Colormap object leaves their output unchanged.
    '''
    import numpy as np

    return [tuple(c) for c in cmap(np.linspace(0, 1, int(n) + 2)[1:-1])[:, :3]]


def categorical(settings, n):
    '''
    N colors for unordered categories, honouring a --style `categorical` palette.

    An unset palette keeps the built-in fallback. A named palette is sampled for n. An
    explicit list is CYCLED when it is shorter than n rather than being swapped for a
    generated palette: a user who wrote the list has said which colors they want, and
    substituting hues they never chose would defeat supplying it.
    '''
    import logging

    spec = (settings or {}).get('categorical')
    if spec is None:
        return categorical_palette(n)
    if isinstance(spec, str):
        import seaborn as sns

        return sns.color_palette(spec, n)
    colors = list(spec)
    if n > len(colors):
        logging.getLogger(__name__).warning(
            f'Style file supplies {len(colors)} categorical colors but {n} categories are '
            f'being drawn; the list is being cycled, so some categories share a color.'
        )
    return [colors[i % len(colors)] for i in range(int(n))]
