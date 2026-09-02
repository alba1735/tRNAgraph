"""Significance thresholds are held once, as p-values, with -log10 derived from them.

p = 0.05 and p = 0.001 are the gold-standard pair this project draws and classifies on. They
used to live as PRE-TRANSFORMED constants in two modules that did not know about each other:
plotsVolcano carried `PVAL_SIG_LINE = 1.3` / `PVAL_STRONG_LINE = 3`, and plotsHeatmap
independently hardcoded matching colorbar ticks `[0, 1.3, 3]` labelled '1'/'0.05'/'<=0.001'.

Two defects followed from that. The pair could drift apart, since nothing tied the two
modules' copies together. And 1.3 is a ROUNDED -log10(0.05) -- the true value is
1.30103 -- so the line labelled p = 0.05 was not drawn at p = 0.05, and points in the sliver
between the two were classified on the wrong side of a threshold the figure claimed.

The expected -log10 values here are written as literals rather than recomputed with math.log10
on purpose: a test that derives the expectation the same way the code does passes by
construction and can never disagree with it.
"""
from trnagraph.modules import plotsThresholds


def test_the_gold_standard_pair_is_held_as_p_values():
    assert plotsThresholds.PVAL_SIG == 0.05
    assert plotsThresholds.PVAL_STRONG == 0.001


def test_log2fc_threshold_is_held_alongside_them():
    assert plotsThresholds.LOG2FC_THRESHOLD == 1.5


def test_neglog10_is_exact_not_rounded():
    assert plotsThresholds.neglog10(0.05) == 1.3010299956639813
    assert plotsThresholds.neglog10(0.001) == 3.0


def test_the_significance_line_is_no_longer_the_rounded_1_point_3():
    """The specific defect this module exists to fix, pinned so it cannot come back."""
    assert plotsThresholds.neglog10(plotsThresholds.PVAL_SIG) != 1.3
