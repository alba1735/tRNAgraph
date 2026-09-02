"""The volcano classified points against a rounded threshold, not the one it labelled.

`PVAL_SIG_LINE = 1.3` is a rounded -log10(0.05); the true value is 1.3010299956639813. The
classification test is `-log10(padj) > PVAL_SIG_LINE`, so every feature whose adjusted p-value
fell between 0.05 and 10**-1.3 = 0.050119 was coloured and labelled as significant while
actually failing the p <= 0.05 the figure's own reference line claimed to draw.

The window is narrow but it is on the wrong side of the line a reader is trusting, and on a
human dataset with thousands of features the tail is dense exactly there.

padj = 0.0501 is inside that window and is the discriminating case: significant under the old
rounded line, correctly not significant under the derived one.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from trnagraph.modules import plotsPalette, plotsVolcano


PAIR = ('A', 'B')


def _frame(padj_by_feature):
    """One feature per padj, every one well past the |log2FC| threshold."""
    return pd.DataFrame(
        {f'log2_{PAIR[0]}-{PAIR[1]}': [2.0] * len(padj_by_feature),
         f'pval_{PAIR[0]}-{PAIR[1]}': list(padj_by_feature.values())},
        index=list(padj_by_feature),
    )


def _facecolor_by_feature(df):
    fig, ax = plt.subplots()
    try:
        plotsVolcano._draw_volcano(ax, df, PAIR, colormap=None, toplabels=0,
                                   feature_types=None, title='t', show_legend=False)
        colors = ax.collections[0].get_facecolors()
        return dict(zip(df.index, [tuple(np.round(c[:3], 6)) for c in colors]))
    finally:
        plt.close(fig)


def _nonsignificant_rgb():
    from matplotlib.colors import to_rgb
    return tuple(np.round(to_rgb(plotsPalette.DE_NONSIGNIFICANT), 6))


def test_a_padj_just_over_0_05_is_not_significant():
    """0.0501 fails p <= 0.05, so it must not be coloured as a hit."""
    colors = _facecolor_by_feature(_frame({'in_sliver': 0.0501}))

    assert colors['in_sliver'] == _nonsignificant_rgb()


def test_a_padj_just_under_0_05_is_still_significant():
    """The control: the fix must not push genuine hits out of the call."""
    colors = _facecolor_by_feature(_frame({'genuine_hit': 0.0499}))

    assert colors['genuine_hit'] != _nonsignificant_rgb()
