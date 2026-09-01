"""Colors live in `plotsPalette`, not scattered through the plot modules.

They used to be literals in several unrelated conventions -- volcano hardcoded hex, heatmap
built its diverging palette inline, seqlogo carried its own nucleotide dict, and eight modules
each reached for `sns.husl_palette(n)` (or `sns.color_palette('husl', n)`, the same call
spelled differently). Making the look uniform meant finding and matching every one by hand.

Centralizing was a pure refactor: 305 of 321 rendered figures came back byte-identical, and
the 16 that did not were the mismatch position plots, which differ between two *identical*
runs because seaborn's stripplot jitter draws from an unseeded global RNG. That is fixed
separately and pinned below.

These tests keep new colors from drifting back into the plot modules.
"""
import ast
import pathlib

import matplotlib.colors as mplcolors
import pytest

from trnagraph.modules import plotsPalette

PLOT_MODULES = sorted(pathlib.Path('src/trnagraph/modules').glob('plots*.py'))
# The palette module is where literals are supposed to live.
CONVERTED = [p for p in PLOT_MODULES if p.name not in ('plotsPalette.py',)
             and not p.name.startswith('plotsLegacy')]


def _string_literals(path):
    tree = ast.parse(path.read_text())
    return [n.value for n in ast.walk(tree)
            if isinstance(n, ast.Constant) and isinstance(n.value, str)]


@pytest.mark.parametrize('path', CONVERTED, ids=lambda p: p.name)
def test_no_plot_module_carries_a_raw_hex_color(path):
    """A new hex literal here means the next uniformity pass has to hunt for it again."""
    offenders = [v for v in _string_literals(path)
                 if v.startswith('#') and len(v) in (4, 7)
                 and all(c in '0123456789abcdefABCDEF' for c in v[1:])]

    assert not offenders, (
        f'{path.name} hardcodes {offenders}. Add a named role to plotsPalette and reference '
        f'it, so a later restyle has one place to change.'
    )


# Saturated CSS color names are always a category or semantic color and belong in the palette.
# Neutrals (black/white/greys) are furniture: several are still literals, which is untidy but
# not the failure this guards against.
SATURATED_NAMES = frozenset({
    'magenta', 'cyan', 'red', 'blue', 'green', 'yellow', 'orange', 'purple', 'pink',
    'lime', 'navy', 'teal', 'olive', 'maroon', 'fuchsia', 'aqua',
})


@pytest.mark.parametrize('path', CONVERTED, ids=lambda p: p.name)
def test_no_plot_module_hardcodes_a_saturated_named_color(path):
    """
    The hex check above misses colors spelled as NAMES, and that gap shipped a real bug:
    plotsCoverage drew read starts/ends from plotsPalette but built its legend swatches from a
    separate hardcoded ['magenta', 'cyan'], so retuning the palette recolored the bars and left
    the legend asserting the old colors. A duplicated color is worse than a scattered one --
    it can disagree with itself.
    """
    offenders = sorted({v for v in _string_literals(path) if v.lower() in SATURATED_NAMES})

    assert not offenders, (
        f'{path.name} hardcodes {offenders}. Add a named role to plotsPalette and reference it, '
        f'so the value cannot drift out of step with the same color used elsewhere.'
    )


@pytest.mark.parametrize('path', CONVERTED, ids=lambda p: p.name)
def test_no_plot_module_calls_the_categorical_palette_directly(path):
    """Two spellings of one concept is what made the fallback hard to change."""
    source = path.read_text()

    assert 'sns.husl_palette(' not in source, f'{path.name}: use plotsPalette.categorical_palette()'
    assert "color_palette('husl'" not in source, f'{path.name}: use plotsPalette.categorical_palette()'


# --- the palette's own contract ------------------------------------------------------

def test_categorical_palette_returns_the_requested_number_of_colors():
    assert len(plotsPalette.categorical_palette(7)) == 7


def test_the_categorical_fallback_is_the_colorblind_safe_set():
    """
    Okabe-Ito extended to ten -- the same values as seaborn's 'colorblind', written out in
    plotsPalette so the default is legible there and cannot shift if seaborn retunes it. This
    pin is what makes the next restyle a deliberate edit rather than a drift.
    """
    import seaborn as sns

    assert plotsPalette.categorical_palette(5) == [
        '#0173b2', '#de8f05', '#029e73', '#d55e00', '#cc78bc']
    assert list(plotsPalette.CATEGORICAL_COLORS) == [
        mplcolors.to_hex(c) for c in sns.color_palette('colorblind', 10)]


def test_the_fallback_switches_to_husl_above_the_safe_palette(caplog):
    """
    Ten is where the colorblind-safe set runs out. Cycling it across 22 amino acids would give
    four groups the same color and make the plot actively misleading, so husl takes over --
    not colorblind-safe, but at least nominally distinct. The warning names the count.
    """
    import seaborn as sns

    n = len(plotsPalette.CATEGORICAL_COLORS) + 1
    assert plotsPalette.categorical_palette(n) == sns.color_palette('husl', n)
    assert 'not colorblind-safe at this count' in caplog.text


def test_no_color_repeats_at_the_boundary():
    """Off-by-one here would silently duplicate a color at exactly ten categories."""
    exact = plotsPalette.categorical_palette(len(plotsPalette.CATEGORICAL_COLORS))

    assert len(set(exact)) == len(plotsPalette.CATEGORICAL_COLORS)


def test_the_diverging_lfc_ramp_is_blue_white_red():
    """
    The previous blue-to-green ramp put the two ends on a deuteranopia confusion pair, making
    the sign of a fold change the hardest thing to read on the plot. Asserted by CHANNEL rather
    than by name: the property that matters is that the ends separate on red vs blue and the
    centre stays neutral, which is what keeps zero readable.
    """
    cmap = plotsPalette.resolve_gradients(None)['lfc']
    low, mid, high = cmap(0.0), cmap(0.5), cmap(1.0)

    assert low[2] > low[0], 'the negative end should be blue-dominant'
    assert high[0] > high[2], 'the positive end should be red-dominant'
    assert min(mid[:3]) > 0.8, 'the centre should be near-neutral so zero reads as zero'


@pytest.mark.parametrize('role,expected', [
    ('correlation', 'crest'), ('significance', 'rocket_r'), ('score', 'mako_r'),
    ('sequence', 'Greys'), ('ordered', 'mako_r'), ('lfc', 'vlag'),
])
def test_each_ramp_is_pinned_to_its_chosen_colormap(role, expected):
    """Changing one of these should be a deliberate edit to this list, not a silent drift."""
    assert plotsPalette.resolve_gradients(None)[role].name == expected


def test_the_two_co_drawn_ramp_pairs_are_distinguishable():
    """
    Only two pairs are ever drawn in one figure, and only those must not be confusable: lfc
    beside significance (the heatmap's two panels), and score beside sequence (the seqlogo).
    """
    import numpy as np

    resolved = plotsPalette.resolve_gradients(None)
    for a, b in [('lfc', 'significance'), ('score', 'sequence')]:
        samples = np.linspace(0.15, 0.85, 6)
        left = np.array([resolved[a](x)[:3] for x in samples])
        right = np.array([resolved[b](x)[:3] for x in samples])
        assert np.abs(left - right).max() > 0.25, f'{a} and {b} are too close to tell apart'


def test_nucleotide_colors_cover_every_base_plus_gap():
    assert set(plotsPalette.NUCLEOTIDE_COLORS) == {'A', 'C', 'G', 'T', 'U', '-'}


def test_u_and_t_share_a_color():
    """Same base; they differ only by whether the plot renders RNA or DNA."""
    assert plotsPalette.NUCLEOTIDE_COLORS['U'] == plotsPalette.NUCLEOTIDE_COLORS['T']


def test_de_direction_colors_are_distinguishable_from_the_nonsignificant_grey():
    assert len({plotsPalette.DE_UP, plotsPalette.DE_DOWN, plotsPalette.DE_NONSIGNIFICANT}) == 3


def _resolve_as_ramp_or_palette(value, name):
    """
    A non-color string constant is either a colormap name or a seaborn palette name. Those are
    two different namespaces -- 'mako_r' is a colormap and not a palette, 'husl' is the reverse
    -- and the style schema already treats them separately, so the check does too. Either
    resolving is a pass; neither means the constant names nothing real.
    """
    import seaborn as sns

    try:
        plotsPalette.build_colormap(value)
        return
    except (ValueError, KeyError):
        pass
    try:
        sns.color_palette(value, 2)
    except Exception:
        raise AssertionError(
            f'{name} = {value!r} is neither a color, a registered colormap, nor a seaborn '
            f'palette name.') from None


def test_every_exported_color_is_parseable_by_matplotlib():
    import matplotlib.colors as mplcolors

    # Private module tables (e.g. _GRADIENT_DEFAULTS, which maps a role to a colormap NAME
    # rather than to a color) are not exported values and are not what this asserts.
    # Note '_NAME'.isupper() is True, so the underscore check has to be explicit.
    # GRADIENT_ROLES lists the role KEYS a style file may set, not colors -- the same category
    # as the private _GRADIENT_DEFAULTS table, which the underscore check already drops.
    # READTYPE_MARKERS/_FALLBACK map a read type to a SHAPE (unicode glyph, matplotlib
    # marker code, legend label) -- shapes are used precisely because colour is already
    # spoken for wherever they are drawn, so there is no colour in them to parse.
    non_color_tables = {'GRADIENT_ROLES', 'READTYPE_MARKERS', 'READTYPE_MARKER_FALLBACK'}
    names = [n for n in dir(plotsPalette)
             if n.isupper() and not n.startswith('_') and n not in non_color_tables
             and not n.endswith(('_ALPHA', '_KWARGS'))]
    for name in names:
        value = getattr(plotsPalette, name)
        if isinstance(value, str):
            # A value here is either a color or the name of a ramp/palette. Resolving whichever
            # it is beats a hardcoded skip list, which went stale the moment a ramp was retuned
            # and would have hidden a genuinely broken name.
            if not mplcolors.is_color_like(value):
                _resolve_as_ramp_or_palette(value, name)
        elif isinstance(value, (dict, tuple, list)):
            for color in (value.values() if isinstance(value, dict) else value):
                mplcolors.to_rgb(color)


def test_palette_module_imports_without_matplotlib_at_module_scope():
    """It is read from anywhere, so it must stay cheap to import."""
    source = pathlib.Path('src/trnagraph/modules/plotsPalette.py').read_text()
    tree = ast.parse(source)

    top_level_imports = [n for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom))]
    assert top_level_imports == [], 'plotsPalette should import lazily inside its helpers'


# --- reproducibility ------------------------------------------------------------------

def test_mismatch_jitter_is_seeded():
    """Unseeded jitter produced a different figure on every regeneration of the same data."""
    from trnagraph.modules import plotsMismatch
    import numpy as np

    before = np.random.get_state()
    with plotsMismatch._seeded_jitter():
        first = np.random.uniform(size=5).tolist()
    with plotsMismatch._seeded_jitter():
        second = np.random.uniform(size=5).tolist()
    after = np.random.get_state()

    assert first == second, 'the same seed must produce the same jitter'
    assert np.array_equal(before[1], after[1]), 'the caller\'s RNG state must be restored'
