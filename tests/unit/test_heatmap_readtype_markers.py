"""The combined heatmap marks each row with its read type.

`visualizer` stacks every --diffrts read type except total into one frame, so a tRNA appears
once per read type under the SAME index label. Sorting is by effect size, which is the plot's
whole purpose, so the copies scatter and nothing says which row came from which read type.

Two properties matter beyond "a glyph appears". The lookup must be POSITIONAL: an earlier
attempt used `tdf.loc[row, 'readtype']`, which returns a Series rather than a scalar exactly
when the index is duplicated -- that is, precisely when the feature is needed. And the
per-read-type heatmaps must be untouched, since a column of identical glyphs teaches nothing.
"""
import numpy as np
import pandas as pd
import pytest

from trnagraph.modules import plotsHeatmap


FEATURES = ['tRNA-Ala-AGC-1', 'tRNA-Gly-GCC-2', 'tRNA-Val-CAC-5']
READTYPES = ['nreads_fiveprime_unique_norm', 'nreads_threeprime_unique_norm',
             'nreads_wholecounts_unique_norm']


def _frame(readtypes=READTYPES):
    """The shape df_combine has: one block of features per read type, stacked."""
    blocks = []
    for offset, readtype in enumerate(readtypes):
        block = pd.DataFrame(
            {'log2_A-B': np.linspace(-2, 2, len(FEATURES)) + offset,
             'pval_A-B': np.linspace(0.001, 0.4, len(FEATURES))},
            index=FEATURES,
        )
        block['readtype'] = readtype
        blocks.append(block)
    return pd.concat(blocks, axis=0)


def _labels(monkeypatch, frame, orientation='vertical'):
    """The tick labels each panel was handed, in draw order."""
    import matplotlib
    matplotlib.use('Agg')

    captured = []
    real_heatmap = plotsHeatmap.sns.heatmap

    def spy(data, *args, **kwargs):
        captured.append({'y': kwargs.get('yticklabels'), 'x': kwargs.get('xticklabels')})
        return real_heatmap(data, *args, **kwargs)

    monkeypatch.setattr(plotsHeatmap.sns, 'heatmap', spy)
    plotsHeatmap.heatmap_plot(frame, 'log2_A-B', 'vlag', 'rocket_r', 25,
                              orientation=orientation)
    return captured


def test_each_row_is_prefixed_with_its_read_types_glyph(monkeypatch):
    labels = _labels(monkeypatch, _frame())[0]['y']

    assert labels is not None, 'the labelled panel must be given explicit labels'
    glyphs = {str(label)[0] for label in labels}
    assert glyphs == {'●', '■', '◆'}, f'expected one glyph per read type, got {glyphs}'
    assert all('tRNA-' in str(label) for label in labels), 'the feature name is still there'


def test_a_single_read_type_heatmap_is_not_marked(monkeypatch):
    """Every row would carry the same glyph, which is noise -- and the per-read-type files
    already name the read type in their filename and title."""
    captured = _labels(monkeypatch, _frame(readtypes=READTYPES[:1]))

    assert captured[0]['y'] is None, 'no explicit labels means seaborn draws the plain index'


def test_duplicate_index_labels_are_resolved_positionally(monkeypatch):
    """The regression that sank the earlier attempt. Every feature appears under all three
    read types, so a label-based lookup cannot tell the copies apart; a correct implementation
    produces each glyph exactly as many times as its read type has rows."""
    labels = _labels(monkeypatch, _frame())[0]['y']

    counts = {glyph: sum(1 for label in labels if str(label).startswith(glyph))
              for glyph in ('●', '■', '◆')}
    assert counts == {'●': len(FEATURES), '■': len(FEATURES), '◆': len(FEATURES)}


def _legend_entries(monkeypatch, frame, orientation='vertical'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plotsHeatmap.heatmap_plot(frame, 'log2_A-B', 'vlag', 'rocket_r', 25,
                              orientation=orientation)
    figure = plt.gcf()
    return [text.get_text() for legend in figure.legends for text in legend.get_texts()]


def test_a_legend_names_the_read_types_that_are_present(monkeypatch):
    """A glyph a reader cannot decode is decoration. The legend uses the published wording."""
    entries = _legend_entries(monkeypatch, _frame())

    assert set(entries) == {"5' counts", "3' counts", 'Whole counts'}


def test_the_legend_lists_only_what_is_drawn(monkeypatch):
    """Naming the whole vocabulary would advertise read types this plot does not contain."""
    entries = _legend_entries(monkeypatch, _frame(readtypes=READTYPES[:2]))

    assert set(entries) == {"5' counts", "3' counts"}
    assert 'Whole counts' not in entries


def test_an_unmarked_heatmap_has_no_legend(monkeypatch):
    assert _legend_entries(monkeypatch, _frame(readtypes=READTYPES[:1])) == []


def test_the_glyph_follows_the_labels_when_the_layout_is_transposed(monkeypatch):
    """Riding in the label rather than being drawn beside it is what makes this work under
    both orientations without any coordinate maths: transposing moves the labels from the
    left panel's y axis to the bottom panel's x axis, and the glyph goes with them."""
    captured = _labels(monkeypatch, _frame(), orientation='horizontal')

    assert captured[0]['x'] is False, 'the top panel shares the axis and shows no labels'
    labels = captured[1]['x']
    assert labels is not None
    # SUFFIXED here, not prefixed. Rotating a label 90 degrees puts its first character at the
    # bottom -- furthest from the cells it describes -- so the glyph goes on the end, which is
    # the end nearest the axis. The eye then travels cell -> glyph -> name in both layouts.
    assert {str(label)[-1] for label in labels} == {'●', '■', '◆'}
    assert all(str(label).startswith('tRNA-') for label in labels)


def test_the_legend_order_is_canonical_not_whatever_sorted_first(monkeypatch):
    """Entries were coming out in order of first appearance in the sorted frame, so the same
    three read types would appear in a different order on the next plot. The order comes from
    the vocabulary instead, which keeps it stable across figures."""
    entries = _legend_entries(monkeypatch, _frame())

    assert entries == ["5' counts", "3' counts", 'Whole counts']
