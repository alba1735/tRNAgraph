"""Layout properties that only break on real data shapes.

Both defects fixed here were invisible on the demo dataset and only appeared on hg38, because
each is a function of how much is being plotted rather than of the plotting code alone: the
heatmap's title gap grew as the heatmap got SHORTER, and the volcano's labels piled up as the
plot got DENSER. Neither is expressible as "does this run without raising", so both are pinned
as measured geometry.
"""
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import pathlib

import numpy as np
import pandas as pd
import pytest

from trnagraph.modules import plotsHeatmap, plotsPalette, plotsVolcano, toolsTG


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close('all')


# --- heatmap title placement ----------------------------------------------------------

def _heatmap_gap(n_rows, n_cols=3):
    """Vertical space between the panels' top and the figure title, in figure coordinates."""
    cols = [f'g{i}_log2' for i in range(n_cols)] + [f'g{i}_pval' for i in range(n_cols)]
    df = pd.DataFrame(np.random.RandomState(0).rand(n_rows, 2 * n_cols), columns=cols,
                      index=[f'feat{i}' for i in range(n_rows)])
    gradients = plotsPalette.resolve_gradients(None)
    plotsHeatmap.heatmap_plot(df, cols[0], gradients['lfc'], gradients['significance'], 25)

    fig = plt.gcf()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    panels_top = max(a.get_tightbbox(renderer).transformed(inv).y1 for a in fig.axes[:2])
    title = next(t for t in fig.texts if 'Heatmap of log2FC' in t.get_text())
    return title.get_position()[1] - panels_top


@pytest.mark.parametrize('n_rows', [4, 10, 20, 50])
def test_the_heatmap_title_sits_just_above_the_panels(n_rows):
    """
    `square=True` sizes the panels from the row count, so in a fixed-height figure a short
    heatmap left the title stranded far above them -- and bbox_inches='tight' cannot trim
    interior white, because the title pins the top edge. Measured before the fix, that band ran
    to 44% of the saved image at 4 rows against 5% at 20.
    """
    assert 0 <= _heatmap_gap(n_rows) < 0.05


def test_the_title_gap_does_not_depend_on_how_much_is_plotted():
    """The actual complaint: the spacing changed with the number of features in the heatmap."""
    gaps = [_heatmap_gap(n) for n in (4, 10, 20, 50)]

    assert max(gaps) - min(gaps) < 0.02, f'gap varies with row count: {gaps}'


# --- volcano label culling ------------------------------------------------------------

def _labelled_boxes(ax, names):
    fig = ax.get_figure()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    wanted = {str(n) for n in names}
    return [t.get_window_extent(renderer) for t in ax.texts if t.get_text() in wanted]


def test_labels_that_cannot_be_placed_are_dropped_rather_than_stacked():
    """
    adjust_text repels; it cannot create space. Asked to fit 100 labels into a crowded volcano
    it converged on a stack of overlapping text that was worse than no labels at all, so any
    label still colliding after the solve is now removed instead of drawn.
    """
    fig, ax = plt.subplots(figsize=(4, 4))
    # 40 features crammed into a tiny region: far more labels than can possibly fit.
    names = [f'tRNA-Feature-Name-{i}' for i in range(40)]
    x = pd.Series(np.linspace(0.0, 0.05, 40), index=names)
    y = pd.Series(np.linspace(0.0, 0.05, 40), index=names)
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)

    kept = plotsVolcano._place_labels(ax, list(names), x, y)

    assert 0 < len(kept) < len(names), f'expected culling, kept {len(kept)} of {len(names)}'


def test_surviving_labels_do_not_overlap_each_other():
    fig, ax = plt.subplots(figsize=(4, 4))
    names = [f'tRNA-Feature-Name-{i}' for i in range(40)]
    x = pd.Series(np.linspace(0.0, 0.05, 40), index=names)
    y = pd.Series(np.linspace(0.0, 0.05, 40), index=names)
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)

    kept = plotsVolcano._place_labels(ax, list(names), x, y)
    boxes = _labelled_boxes(ax, kept)

    overlaps = [(i, j) for i in range(len(boxes)) for j in range(i + 1, len(boxes))
                if boxes[i].overlaps(boxes[j])]
    assert not overlaps, f'{len(overlaps)} pairs of labels still overlap'


def test_culling_keeps_the_highest_priority_labels():
    """`names` arrives ordered most-interesting-first, so the survivors must come off the top."""
    fig, ax = plt.subplots(figsize=(4, 4))
    names = [f'tRNA-Feature-Name-{i}' for i in range(40)]
    x = pd.Series(np.linspace(0.0, 0.05, 40), index=names)
    y = pd.Series(np.linspace(0.0, 0.05, 40), index=names)
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)

    kept = plotsVolcano._place_labels(ax, list(names), x, y)

    assert kept[0] == names[0]
    assert kept == [n for n in names if n in set(kept)], 'priority order was not preserved'


def test_a_sparse_plot_keeps_every_label():
    """Culling must not fire when there is room; only density should cost a label."""
    fig, ax = plt.subplots(figsize=(6, 6))
    names = ['tRNA-Alpha-1', 'tRNA-Beta-2', 'tRNA-Gamma-3']
    x = pd.Series([-3.0, 0.0, 3.0], index=names)
    y = pd.Series([8.0, 4.0, 1.0], index=names)
    ax.set_xlim(-5, 5)
    ax.set_ylim(0, 10)

    assert plotsVolcano._place_labels(ax, list(names), x, y) == names


def test_a_label_that_did_not_move_gets_no_connector():
    """A line drawn under a label sitting on its own point is noise, not guidance."""
    fig, ax = plt.subplots(figsize=(6, 6))
    names = ['tRNA-Alpha-1']
    x = pd.Series([0.0], index=names)
    y = pd.Series([5.0], index=names)
    ax.set_xlim(-5, 5)
    ax.set_ylim(0, 10)

    plotsVolcano._place_labels(ax, list(names), x, y)
    moved = np.hypot(*(ax.transData.transform(ax.texts[0].get_position())
                       - ax.transData.transform((0.0, 5.0))))

    if moved <= plotsVolcano.MIN_CONNECTOR_PX:
        assert not [a for a in ax.texts if a.get_text() == ''], 'connector drawn for a still label'


# --- figsize actually reaching individual plots ---------------------------------------

def test_figsize_for_prefers_the_configured_value():
    assert toolsTG.figsize_for({'figsize': (3.25, 3.25)}, (6, 6)) == (3.25, 3.25)


@pytest.mark.parametrize('settings', [None, {}, {'figsize': None}, {'dpi': 300}])
def test_figsize_for_falls_back_to_the_modules_default(settings):
    assert toolsTG.figsize_for(settings, (6, 6)) == (6, 6)


def test_a_configured_figsize_changes_the_rendered_figure():
    """
    style_context() puts figsize into rcParams, but every plot module passes an explicit
    figsize= to plt.figure/plt.subplots and an explicit argument beats the rcParam -- so the
    setting silently did nothing, for every graph type, not just the one it was noticed on.
    """
    import seaborn as sns

    plt.figure(figsize=toolsTG.figsize_for({'figsize': (3.25, 3.25)}, (6, 6)))
    sns.heatmap(pd.DataFrame(np.random.RandomState(0).rand(4, 4)))

    assert tuple(plt.gcf().get_size_inches()) == (3.25, 3.25)


INDIVIDUAL_PLOT_MODULES = ['plotsVolcano.py', 'plotsCorrelation.py', 'plotsPca.py',
                           'plotsCount.py', 'plotsCompare.py', 'plotsMismatch.py',
                           'plotsCluster.py', 'plotsRadar.py', 'plotsHeatmap.py',
                           'plotsCoverage.py', 'plotsSeqlogo.py']


@pytest.mark.parametrize('name', INDIVIDUAL_PLOT_MODULES)
def test_every_module_with_an_individual_plot_consults_figsize_for(name):
    """
    A module that sizes a single-plot figure must go through the helper. Combined and
    multi-page pages deliberately keep their computed geometry, so this only asserts that the
    helper is referenced at all -- enough to catch a module that was never wired up.
    """
    source = pathlib.Path('src/trnagraph/modules').joinpath(name).read_text()

    assert 'figsize_for(' in source, (
        f'{name} hardcodes every figsize, so a --style figsize cannot reach it. Individual '
        f'plots should size themselves with toolsTG.figsize_for(settings, <default>).'
    )


# --- where square is preserved and where an explicit size wins ------------------------

def test_a_combined_page_keeps_square_panels_even_with_a_figsize():
    """A grid of stretched panels reads badly, and the page computes its own geometry anyway."""
    assert plotsVolcano._keep_square(standalone=False, settings={'figsize': (9.75, 6.5)})


def test_a_standalone_volcano_with_an_explicit_size_is_not_forced_square():
    """Otherwise the extra width of a wide figure is discarded, then cropped away on save."""
    assert not plotsVolcano._keep_square(standalone=True, settings={'figsize': (9.75, 6.5)})


@pytest.mark.parametrize('settings', [None, {}, {'figsize': None}, {'dpi': 300}])
def test_a_standalone_volcano_stays_square_by_default(settings):
    assert plotsVolcano._keep_square(standalone=True, settings=settings)


@pytest.mark.parametrize('name', ['plotsPca.py', 'plotsCorrelation.py'])
def test_pca_and_correlation_stay_square_unconditionally(name):
    """
    figsize sets how big the square is for these two, it does not unsquare them: a correlation
    matrix is square because it is symmetric, not as a preference.
    """
    source = pathlib.Path('src/trnagraph/modules').joinpath(name).read_text()

    assert 'set_box_aspect(1)' in source
    assert '_keep_square' not in source, f'{name} must not make its square conditional'
