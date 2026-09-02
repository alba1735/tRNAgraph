"""The heatmap's significance colorbar carried its own copy of the threshold pair.

`cbar.set_ticks([0, 1.3, 3])` with labels `['1', '0.05', '<=0.001']` restated, as pre-computed
-log10 literals, the same two thresholds plotsVolcano held separately -- so the two modules
could drift apart, and both inherited the rounded 1.3 that is not -log10(0.05).

The panel's colour ceiling is the same number again: `sns.heatmap(pval_tdf, ..., vmax=3)` is
-log10(0.001), so the scale saturates exactly at the strong-evidence threshold. Deriving all
three from one source is what stops a future change to the pair moving the ticks without
moving the scale under them.
"""
from trnagraph.modules import plotsHeatmap, plotsThresholds


def test_ticks_are_derived_positions_not_rounded_literals():
    ticks, _ = plotsHeatmap.significance_colorbar_ticks()

    assert ticks == [0.0, 1.3010299956639813, 3.0]


def test_labels_name_the_p_values_those_positions_stand_for():
    _, labels = plotsHeatmap.significance_colorbar_ticks()

    assert labels == ['1', '0.05', '<=0.001']


def test_the_panel_ceiling_is_the_strong_threshold():
    """vmax saturates the scale at p = 0.001; it must move with the threshold, not be a literal."""
    assert plotsHeatmap.SIGNIFICANCE_VMAX == plotsThresholds.neglog10(plotsThresholds.PVAL_STRONG)
