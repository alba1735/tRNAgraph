"""The specificity overview's legend was anchored to the figure corner, not the plots.

`fig.legend(..., bbox_to_anchor=(0.995, 0.995))` puts the legend at the top of the canvas,
but the 4x4 subplot grid ends at `subplotpars.top` (0.88 by default). The legend was
therefore stranded in a tall empty band, which `bbox_inches='tight'` preserved because the
legend itself was the topmost artist. Anchoring to the grid instead keeps it beside the
plots it labels.
"""
import matplotlib

matplotlib.use("Agg")

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

from trnagraph.modules.plotsCoverage import LEGEND_TITLE_CLEARANCE


def _legend_anchor(fig):
    """Reproduce the anchor the overview uses, so the rule is asserted in one place."""
    return (fig.subplotpars.right, fig.subplotpars.top + LEGEND_TITLE_CLEARANCE)


def test_clearance_keeps_the_legend_above_the_axes_but_inside_the_canvas():
    assert LEGEND_TITLE_CLEARANCE > 0, "must clear the top row's axis titles"
    assert LEGEND_TITLE_CLEARANCE < 0.05, "a larger gap reopens the dead band it replaced"


def test_anchor_tracks_the_subplot_grid_not_the_figure_corner():
    fig = plt.figure(figsize=(24, 22))
    fig.add_subplot(4, 4, 1)
    try:
        x, y = _legend_anchor(fig)

        assert y < 0.95, "0.995 was the stranded-in-a-dead-band placement"
        assert y > fig.subplotpars.top, "must sit above the axes, not on top of them"
        assert x == fig.subplotpars.right, "right edge follows the grid, not the canvas"
    finally:
        plt.close(fig)


def test_legend_sits_close_above_the_topmost_axes():
    """The regression: a large vertical gap between the legend and the first plot row."""
    fig = plt.figure(figsize=(24, 22))
    axes = [fig.add_subplot(4, 4, i + 1) for i in range(16)]
    handles = [mpatches.Patch(color="C0", label=f"cat{i}") for i in range(4)]
    try:
        legend = fig.legend(handles=handles, loc="lower right", frameon=False,
                            ncol=len(handles), bbox_to_anchor=_legend_anchor(fig))
        fig.canvas.draw()

        transform = fig.transFigure.inverted()
        legend_bottom = transform.transform(legend.get_window_extent())[0][1]
        top_axes = max(ax.get_position().y1 for ax in axes)

        gap = legend_bottom - top_axes
        assert 0 <= gap < 0.05, f"legend floats {gap:.3f} of the figure above the plots"
    finally:
        plt.close(fig)


def test_legend_is_one_horizontal_row():
    """Four categories stacked vertically at the corner is what made it read as detached."""
    fig = plt.figure(figsize=(24, 22))
    fig.add_subplot(4, 4, 1)
    handles = [mpatches.Patch(color="C0", label=f"cat{i}") for i in range(4)]
    try:
        legend = fig.legend(handles=handles, loc="lower right", frameon=False,
                            ncol=len(handles), bbox_to_anchor=_legend_anchor(fig))
        fig.canvas.draw()
        extent = legend.get_window_extent()

        assert extent.width > extent.height, "expected a single wide row, not a tall stack"
    finally:
        plt.close(fig)
