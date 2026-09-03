"""The agreement volcano itself.

One contrast on the x-axis, read type in the marker shape, agreement count in the colour. The
encodings have to stay separable: a reader recovers "which read type" from the shape and "how
consistent" from the colour, so neither may leak into the other.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from matplotlib.colors import to_hex

from trnagraph.modules import plotsAgreement, plotsPalette


REF, OTHER = 'Day 0', 'Day 70'
DRAWN = (REF, OTHER)
CONTRASTS = [(REF, 'Day 35'), (REF, OTHER)]
COLORS = {REF: '#007FFF', OTHER: '#FFDC49'}


def _frame(features, per_contrast):
    data = {}
    for (a, b), (log2fc, pval) in per_contrast.items():
        data[f'log2_{a}-{b}'] = pd.Series(log2fc, index=features)
        data[f'pval_{a}-{b}'] = pd.Series(pval, index=features)
    return pd.DataFrame(data, index=features)


def _table(readtypes=('fiveprime',)):
    """Three features: consistently up, consistently down, and quiet."""
    features = ['up-both', 'down-both', 'quiet']
    frame = _frame(features, {
        (REF, 'Day 35'): ([4.0, -4.0, 0.1], [1e-8, 1e-8, 0.9]),
        (REF, OTHER): ([5.0, -5.0, 0.2], [1e-9, 1e-9, 0.8]),
    })
    return plotsAgreement.agreement_table({rt: frame for rt in readtypes}, CONTRASTS, DRAWN,
                                          log2fc=1.5, padj=0.001)


def _draw(table, **kw):
    fig, ax = plt.subplots()
    kw.setdefault('colormap', COLORS)
    plotsAgreement.draw_agreement(ax, table, DRAWN, log2fc=1.5, padj=0.001, **kw)
    return fig, ax


def test_each_read_type_gets_its_own_marker_shape():
    """Shape is the read-type channel, taken from the palette the heatmap and legend share."""
    table = _table(('fiveprime', 'threeprime', 'wholecounts'))
    fig, ax = _draw(table)
    try:
        drawn = {tuple(c.get_paths()[0].vertices.round(4).flatten()) for c in ax.collections}
    finally:
        plt.close(fig)

    assert len(drawn) == 3, 'three read types must be three distinguishable shapes'


def test_a_feature_with_no_significance_call_is_drawn_in_the_neutral_colour():
    fig, ax = _draw(_table())
    try:
        colors = [to_hex(c) for coll in ax.collections
                  for c in coll.get_facecolors()]
    finally:
        plt.close(fig)

    assert to_hex(plotsAgreement.UNCALLED_COLOR) in colors


def test_the_two_directions_take_different_ramps():
    """Up and down must never share a colour, or the figure says nothing about direction."""
    up = plotsAgreement.direction_ramp(OTHER, 2, colormap=COLORS)
    down = plotsAgreement.direction_ramp(REF, 2, colormap=COLORS)

    assert not set(map(to_hex, up)) & set(map(to_hex, down))


def test_the_strongest_tier_wears_the_level_colour_on_the_figure():
    fig, ax = _draw(_table())
    try:
        colors = {to_hex(c) for coll in ax.collections for c in coll.get_facecolors()}
    finally:
        plt.close(fig)

    assert to_hex(COLORS[OTHER]) in colors, 'up-both is 2 of 2 and should wear Day 70'
    assert to_hex(COLORS[REF]) in colors, 'down-both is 2 of 2 and should wear Day 0'


def test_both_p_value_lines_and_both_fold_change_lines_are_drawn():
    fig, ax = _draw(_table())
    try:
        horizontals = sorted(round(line.get_ydata()[0], 4) for line in ax.lines
                             if len(set(line.get_ydata())) == 1)
        verticals = sorted(round(line.get_xdata()[0], 4) for line in ax.lines
                           if len(set(line.get_xdata())) == 1)
    finally:
        plt.close(fig)

    from trnagraph.modules import plotsThresholds
    assert round(plotsThresholds.neglog10(0.05), 4) in horizontals
    assert round(plotsThresholds.neglog10(0.001), 4) in horizontals
    assert [-1.5, 1.5] == [v for v in verticals if abs(v) == 1.5]


def test_both_thresholds_are_printed_and_distinguished():
    """D26, sharpened by the two-threshold split: printing only 0.001 while colouring points
    down to 0.05 is a caption that contradicts the figure. Both appear, each saying its job."""
    fig, ax = _draw(_table())
    try:
        text = ' '.join(t.get_text() for t in ax.texts) + ' '.join(
            t.get_text() for t in fig.texts)
        legend = ax.get_legend()
        if legend:
            text += ' '.join(t.get_text() for t in legend.get_texts())
    finally:
        plt.close(fig)

    assert '0.001' in text, 'the agreement threshold'
    assert '0.05' in text, 'the threshold that decides which points are coloured at all'
    assert '1.5' in text


def test_a_feature_is_labelled_once_even_when_it_appears_in_three_read_types():
    """D15: one tRNA must not consume the label budget three times over."""
    table = _table(('fiveprime', 'threeprime', 'wholecounts'))
    fig, ax = _draw(table, toplabels=100)
    try:
        labels = [t.get_text() for t in ax.texts if t.get_text() in
                  {'up-both', 'down-both', 'quiet'}]
    finally:
        plt.close(fig)

    assert sorted(labels) == ['down-both', 'up-both'], 'one label per feature, significant only'


def test_the_legend_names_tiers_as_literal_counts():
    """D14: "2 of 2" rather than an invented tier name that stops being true at N=5."""
    fig, ax = _draw(_table())
    try:
        legend = ax.get_legend()
        text = ' '.join(t.get_text() for t in legend.get_texts())
    finally:
        plt.close(fig)

    assert '2 of 2' in text


def test_an_empty_table_is_refused_rather_than_drawn_blank():
    with pytest.raises(plotsAgreement.InvalidAgreementError):
        _draw(_table().iloc[0:0])


def test_the_read_type_legend_uses_the_shared_marker_vocabulary():
    table = _table(('fiveprime', 'threeprime'))
    fig, ax = _draw(table)
    try:
        text = ' '.join(t.get_text() for t in ax.get_legend().get_texts())
    finally:
        plt.close(fig)

    assert plotsPalette.READTYPE_MARKERS['fiveprime']['label'] in text
    assert plotsPalette.READTYPE_MARKERS['threeprime']['label'] in text
