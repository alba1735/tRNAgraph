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


@pytest.mark.parametrize('path', CONVERTED, ids=lambda p: p.name)
def test_no_plot_module_calls_the_categorical_palette_directly(path):
    """Two spellings of one concept is what made the fallback hard to change."""
    source = path.read_text()

    assert 'sns.husl_palette(' not in source, f'{path.name}: use plotsPalette.categorical_palette()'
    assert "color_palette('husl'" not in source, f'{path.name}: use plotsPalette.categorical_palette()'


# --- the palette's own contract ------------------------------------------------------

def test_categorical_palette_returns_the_requested_number_of_colors():
    assert len(plotsPalette.categorical_palette(7)) == 7


def test_categorical_palette_matches_what_the_modules_used_before():
    """Refactor, not restyle: the fallback must be the same colors as before."""
    import seaborn as sns

    assert plotsPalette.categorical_palette(5) == sns.husl_palette(5)


def test_diverging_lfc_cmap_matches_the_previous_inline_definition():
    import seaborn as sns

    expected = sns.diverging_palette(255, 85, s=255, l=70, sep=20, as_cmap=True)
    actual = plotsPalette.diverging_lfc_cmap()

    assert actual(0.0) == expected(0.0)
    assert actual(0.5) == expected(0.5)
    assert actual(1.0) == expected(1.0)


def test_nucleotide_colors_cover_every_base_plus_gap():
    assert set(plotsPalette.NUCLEOTIDE_COLORS) == {'A', 'C', 'G', 'T', 'U', '-'}


def test_u_and_t_share_a_color():
    """Same base; they differ only by whether the plot renders RNA or DNA."""
    assert plotsPalette.NUCLEOTIDE_COLORS['U'] == plotsPalette.NUCLEOTIDE_COLORS['T']


def test_de_direction_colors_are_distinguishable_from_the_nonsignificant_grey():
    assert len({plotsPalette.DE_UP, plotsPalette.DE_DOWN, plotsPalette.DE_NONSIGNIFICANT}) == 3


def test_every_exported_color_is_parseable_by_matplotlib():
    import matplotlib.colors as mplcolors

    names = [n for n in dir(plotsPalette) if n.isupper() and not n.endswith(('_ALPHA', '_KWARGS'))]
    for name in names:
        value = getattr(plotsPalette, name)
        if isinstance(value, str) and value not in ('husl', 'Blues', 'Greens', 'Greys', 'mako_r'):
            mplcolors.to_rgb(value)
        elif isinstance(value, dict):
            for color in value.values():
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
