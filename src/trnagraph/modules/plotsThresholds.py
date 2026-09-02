#!/usr/bin/env python3
'''
Every significance threshold tRNAgraph draws or classifies on, in one place.

p = 0.05 and p = 0.001 are the gold-standard pair for this field, and both are drawn on every
figure that shows a p-value. Which of the two CLASSIFIES, as opposed to merely being drawn, is
a per-plot choice and is made at the call site: the volcano classifies at PVAL_SIG and draws
PVAL_STRONG as a reference line, while the set-membership plots classify at PVAL_STRONG,
because they make a stronger claim than a single comparison does.

Thresholds are held as P-VALUES and the axis positions derived, rather than stored
pre-transformed. Before this module they lived as -log10 literals in two modules that did not
know about each other -- plotsVolcano's `PVAL_SIG_LINE = 1.3` / `PVAL_STRONG_LINE = 3`, and
plotsHeatmap's colorbar ticks `[0, 1.3, 3]` -- which let the pair drift apart and, worse,
rounded -log10(0.05) to 1.3 when it is 1.30103. The line labelled p = 0.05 was therefore not
drawn at p = 0.05, and a point in the sliver between the two was classified on the wrong side
of the threshold the figure claimed.

A sibling of plotsPalette rather than part of it: that module is documented as holding colors
and nothing else, and a significance threshold is not a color. Like plotsPalette, this holds
values only -- no plotting logic and no matplotlib import at module scope -- so it stays cheap
to import and safe to read from anywhere.

Note the p-values these are compared against are BH-ADJUSTED (adataLog2FC stores `padj`, not
the raw per-comparison p-value), so p <= 0.001 here is a stricter claim than the bare number
suggests.
'''

import math

#: The conventional significance threshold, and what the volcano classifies on.
PVAL_SIG = 0.05
#: The strong-evidence threshold. Drawn on the volcano; classified on by the set-membership plots.
PVAL_STRONG = 0.001
#: Minimum |log2 fold change| for a feature to count as differentially expressed.
LOG2FC_THRESHOLD = 1.5


def neglog10(pvalue):
    '''
    A p-value's position on a -log10 axis.

    Derived rather than tabulated so a threshold and the line drawn for it cannot disagree.
    '''
    return -math.log10(pvalue)
