"""Regression tests for roadmap.md's volcano-stage crash: `_draw_volcano()` passed the full
log2FC / -log10(padj) columns straight into `adjustText.adjust_text()`, which builds a
`scipy.spatial.KDTree` over every point (not just labeled ones) to repel labels -- any NaN
(from a gene DESeq2's independent filtering/Cook's-outlier step excluded, i.e. never tested,
not "insignificant") or Inf (from a padj that underflowed to exactly 0.0 on an extremely
significant hit) raised `ValueError: data must be finite`. Reproduced against the real hg38
server run's cluster.h5ad. Fix: drop NaN rows (safe -- they were never given a real p-value)
and clip a zero padj away from zero instead of dropping it (an underflowed p-value is a real,
extremely significant result, not something to silently lose from the plot)."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from trnagraph.modules.plotsVolcano import _draw_volcano

PAIR = ('A', 'B')


def _make_df_pairs():
    # adjustText's internal repel step (the one that crashes on non-finite input) only kicks in
    # once there are enough labeled points to potentially overlap -- a handful of significant
    # rows (not just one) is needed to reliably exercise it, matching the real large-dataset repro.
    rng = np.random.default_rng(0)
    n_sig = 10
    index = ['trnaNaN', 'trnaZeroP'] + [f'trnaSig{i}' for i in range(n_sig)] + ['trnaNormal1', 'trnaNormal2']
    log2fc = [np.nan, 3.0] + list(rng.choice([-4.0, 4.0], size=n_sig)) + [0.2, -2.0]
    pval = [np.nan, 0.0] + list(rng.uniform(1e-6, 1e-4, size=n_sig)) + [0.5, 0.001]
    return pd.DataFrame({'log2_A-B': log2fc, 'pval_A-B': pval}, index=index)


def test_draw_volcano_does_not_crash_on_nan_and_zero_pvalues():
    df_pairs = _make_df_pairs()
    fig, ax = plt.subplots()
    try:
        _draw_volcano(ax, df_pairs, PAIR, colormap=None, toplabels=None, feature_types=None, title='test')
    finally:
        plt.close(fig)


def test_draw_volcano_drops_nan_rows_from_the_plotted_points():
    df_pairs = _make_df_pairs()
    fig, ax = plt.subplots()
    try:
        _draw_volcano(ax, df_pairs, PAIR, colormap=None, toplabels=None, feature_types=None, title='test')
        offsets = ax.collections[-1].get_offsets()
    finally:
        plt.close(fig)

    assert len(offsets) == len(df_pairs) - 1, "the NaN row (untested gene) must not appear among the plotted points"
    assert np.isfinite(offsets).all()


def test_draw_volcano_preserves_zero_pvalue_point_instead_of_dropping_it():
    """A padj of exactly 0.0 is an extremely significant real result (float64 underflow), not
    noise -- it must still be plotted (with a large but finite y), not silently discarded."""
    df_pairs = _make_df_pairs()
    fig, ax = plt.subplots()
    try:
        _draw_volcano(ax, df_pairs, PAIR, colormap=None, toplabels=None, feature_types=None, title='test')
        offsets = ax.collections[-1].get_offsets()
    finally:
        plt.close(fig)

    ys = np.asarray(offsets)[:, 1]
    assert np.isfinite(ys).all()
    assert ys.max() > 100, "the zero-pvalue point should plot with a large finite -log10(p), not be clipped away to near zero"
