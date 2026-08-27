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
