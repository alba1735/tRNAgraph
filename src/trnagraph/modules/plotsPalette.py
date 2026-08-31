#!/usr/bin/env python3
'''
Every color tRNAgraph draws with, in one place.

Colors used to be literals scattered across the plots*.py modules, in several unrelated
conventions: volcano hardcoded hex, heatmap built a diverging palette inline, seqlogo carried
its own nucleotide dict, and eight modules each called `sns.husl_palette(n)` (or
`sns.color_palette('husl', n)`, the same thing spelled differently) as their categorical
fallback. Making the look uniform meant finding and matching all of them by hand.

This module holds the values only -- no plotting logic and no matplotlib import at module
scope, so it stays cheap to import and safe to read from anywhere. The values here are exactly
what the modules used before centralization: this is a refactor, not a restyle. Changing the
defaults is a separate, deliberate piece of work (see docs/roadmap.md).

Naming is by ROLE rather than by appearance -- `REFERENCE_LINE` not `BLACK` -- so that a later
pass can change what a role looks like without every call site reading as a lie.
'''

# --- Categorical fallback ------------------------------------------------------------
# Used whenever a plot needs N visually distinct colors for unordered categories and the user
# supplied no palette for that column. One name for what were two spellings of the same call.
CATEGORICAL_PALETTE = 'husl'


def categorical_palette(n):
    '''N distinct colors for unordered categories -- the shared fallback when --style sets none.'''
    import seaborn as sns

    return sns.husl_palette(n)


# --- Continuous / ordered scales -----------------------------------------------------
# Named by what they encode, so a future pass can retune one without disturbing the others.
SEQUENTIAL_CORRELATION = 'Blues'    # correlation R^2 heatmap
SEQUENTIAL_SIGNIFICANCE = 'Greens'  # -log10(p) heatmap panel
SEQUENTIAL_SCORE = 'Blues'          # seqlogo per-position score heatmap and its colorbar
SEQUENTIAL_SEQUENCE = 'Greys'       # seqlogo sequence heatmap background
# Ordered specificity scale (coverage partition) and cluster numeric coloring: light = least,
# dark = most, which is why it is a sequential ramp rather than categorical hues.
SEQUENTIAL_ORDERED = 'mako_r'

# Diverging scale for log2 fold change, centered on zero.
DIVERGING_LFC_KWARGS = dict(h_neg=255, h_pos=85, s=255, l=70, sep=20)


def diverging_lfc_cmap():
    '''Diverging colormap for log2 fold change, centered on zero.'''
    import seaborn as sns

    return sns.diverging_palette(as_cmap=True, **DIVERGING_LFC_KWARGS)


# --- Differential-expression direction -----------------------------------------------
# Volcano point colors when --style names no color for the compared groups.
DE_UP = '#d62728'
DE_DOWN = '#1f77b4'
DE_NONSIGNIFICANT = '#7f7f7f'

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
# Read starts and ends overlaid on a coverage trace; deliberately loud, since they mark
# fragment boundaries against the coverage fill behind them.
READSTART_COLOR = 'magenta'
READEND_COLOR = 'cyan'

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

# Role -> the default this module defines for it. `lfc` is built rather than named, so the
# mapping holds a callable there; _default_gradient() collapses the difference.
GRADIENT_ROLES = ('correlation', 'significance', 'score', 'sequence', 'ordered', 'lfc')

_GRADIENT_DEFAULTS = {
    'correlation': SEQUENTIAL_CORRELATION,
    'significance': SEQUENTIAL_SIGNIFICANCE,
    'score': SEQUENTIAL_SCORE,
    'sequence': SEQUENTIAL_SEQUENCE,
    'ordered': SEQUENTIAL_ORDERED,
    'lfc': diverging_lfc_cmap,
}


def _default_gradient(role):
    '''
    The built-in colormap for one role, as a Colormap object.

    Goes through build_colormap() rather than plt.get_cmap() directly because SEQUENTIAL_ORDERED
    is `mako_r`, which seaborn registers on import -- resolving it without that import raises
    'not a valid value for name' on any path that has not already pulled seaborn in.
    '''
    default = _GRADIENT_DEFAULTS[role]
    if callable(default):
        return default()
    return build_colormap(default)


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
