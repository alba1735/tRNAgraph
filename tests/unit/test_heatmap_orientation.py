"""`--heatorient` lays the heatmap out for a landscape page.

The default puts features on rows and comparisons on columns, in two `square=True` panels side
by side. That is the wrong aspect once the feature count grows: a 50-row heatmap beside a
second 50-row heatmap is tall and narrow. `horizontal` transposes the data -- features to
columns, comparisons to rows -- and stacks the two panels vertically.

`vertical` stays the default and must keep today's output exactly, filenames included.
"""
import numpy as np
import pandas as pd
import pytest

from trnagraph.modules import plotsHeatmap, toolsTG


FEATURES = [f'tRNA-Ala-AGC-{i}' for i in range(6)]
COMPARISONS = ['A-B', 'A-C', 'B-C']


def _log2fc_frame():
    """The shape toolsTG.adataLog2FC returns: features on the index, log2_/pval_ column pairs."""
    data = {}
    for i, pair in enumerate(COMPARISONS):
        data[f'log2_{pair}'] = np.linspace(-3, 3, len(FEATURES)) + i
        data[f'pval_{pair}'] = np.linspace(0.001, 0.5, len(FEATURES))
    frame = pd.DataFrame(data, index=FEATURES)
    frame['readtype'] = 'nreads_total_unique_norm'
    return frame


def _drawn(monkeypatch, orientation):
    """Render one heatmap, capturing the frames handed to seaborn and the panel geometry."""
    import matplotlib
    matplotlib.use('Agg')

    captured = {'frames': []}
    real_heatmap = plotsHeatmap.sns.heatmap

    def spy(data, *args, **kwargs):
        captured['frames'].append(data)
        captured.setdefault('cbar_kws', []).append(kwargs.get('cbar_kws'))
        captured.setdefault('axes', []).append(kwargs.get('ax'))
        return real_heatmap(data, *args, **kwargs)

    monkeypatch.setattr(plotsHeatmap.sns, 'heatmap', spy)
    plotsHeatmap.heatmap_plot(_log2fc_frame(), f'log2_{COMPARISONS[0]}', 'vlag', 'rocket_r',
                              25, orientation=orientation)
    # The panels seaborn was actually handed. Reading the gridspec instead would report the
    # colorbar's own geometry: sns.heatmap's colorbar rewrites the parent gridspec, so a
    # two-panel figure reports (3, 2).
    first, second = (ax.get_position() for ax in captured['axes'])
    captured['side_by_side'] = second.x0 > first.x1 - 1e-6
    captured['stacked'] = second.y1 < first.y0 + 1e-6
    return captured


def test_the_default_puts_features_on_rows_in_two_panels_side_by_side():
    """Pinned so `vertical` cannot drift while `horizontal` is being built."""
    import matplotlib
    matplotlib.use('Agg')

    captured = _drawn(pytest.MonkeyPatch(), 'vertical')

    assert captured['side_by_side'], 'the pval panel sits to the right of the log2FC panel'
    assert not captured['stacked']
    log_frame = captured['frames'][0]
    # A set, not a list: the rows are sorted by the selected comparison, which is what puts
    # the strongest fold changes at the top and is not what this test is about.
    assert set(log_frame.index) == set(FEATURES)
    assert all('log2_' in str(c) for c in log_frame.columns)


def test_horizontal_transposes_the_data_and_stacks_the_panels():
    captured = _drawn(pytest.MonkeyPatch(), 'horizontal')

    assert captured['stacked'], 'the pval panel sits below the log2FC panel'
    assert not captured['side_by_side']
    log_frame = captured['frames'][0]
    assert set(log_frame.columns) == set(FEATURES), 'features become the x axis'
    assert all('log2_' in str(i) for i in log_frame.index), 'comparisons become the rows'


def test_an_unrecognised_orientation_is_a_usage_error_naming_the_choices():
    with pytest.raises(toolsTG.InvalidParameterError) as raised:
        plotsHeatmap.heatmap_plot(_log2fc_frame(), f'log2_{COMPARISONS[0]}', 'vlag',
                                  'rocket_r', 25, orientation='sideways')
    message = str(raised.value)
    assert '--heatorient' in message
    assert 'vertical' in message and 'horizontal' in message


def _adata():
    """A small real object, so the visualizer runs its actual DE path rather than a stand-in."""
    import anndata as ad

    rng = np.random.default_rng(0)
    samples = [f'{g}_{i}' for g in ('A', 'B') for i in range(3)]
    rows = []
    for trna in FEATURES:
        for sample in samples:
            group = sample.split('_')[0]
            mean = 400 if group == 'A' else 150
            raw = int(rng.negative_binomial(n=10, p=10 / (10 + mean)))
            rows.append({'trna': trna, 'sample': sample, 'group': group,
                         'nreads_total_unique_raw': raw, 'nreads_total_unique_norm': float(raw)})
    obs = pd.DataFrame(rows)
    obs.index = [f'{r.trna}_{r["sample"]}' for _, r in obs.iterrows()]
    return ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs)


def _run_visualizer(tmp_path, orientation):
    import matplotlib
    matplotlib.use('Agg')

    output = f'{tmp_path}/'
    plotsHeatmap.visualizer(_adata(), 'group', ['nreads_total_unique_norm'], 0, 25, True,
                            output, orientation=orientation)
    return sorted(p.name for p in tmp_path.rglob('*') if p.is_file())


def test_the_default_orientation_keeps_todays_filenames(tmp_path):
    written = _run_visualizer(tmp_path, 'vertical')

    assert any(name.endswith('_heatmap.pdf') for name in written)
    assert not any('horizontal' in name for name in written), (
        'adding the option must not relocate the output anyone already has')


def test_a_horizontal_run_records_the_orientation_in_its_pdf_names_but_not_the_csv(tmp_path):
    """The PDFs are figures and the two orientations are different pictures of the same data,
    so both can exist. The CSV is the DATA -- transposing or renaming it would make two runs
    of one analysis produce exports that no longer diff against each other."""
    written = _run_visualizer(tmp_path, 'horizontal')

    pdfs = [name for name in written if name.endswith('.pdf')]
    csvs = [name for name in written if name.endswith('.csv')]
    assert pdfs and csvs
    assert all('_horizontal' in name for name in pdfs)
    assert not any('horizontal' in name for name in csvs)


def test_the_csv_holds_features_on_rows_whatever_the_orientation(tmp_path):
    _run_visualizer(tmp_path, 'horizontal')
    csv = next(p for p in tmp_path.rglob('*.csv'))

    frame = pd.read_csv(csv, index_col=0)
    assert set(frame.index) <= set(FEATURES), 'the export keeps features on the index'
    assert any('log2_' in str(c) for c in frame.columns)


def test_the_cli_exposes_heatorient_and_defaults_to_vertical():
    import inspect

    from trnagraph import cli

    parameter = inspect.signature(cli.graph).parameters['heatorient']
    assert parameter.default.default == 'vertical'


def test_heatorient_is_settable_from_a_config_file():
    from trnagraph.modules.toolsSchemas import GraphFlags

    assert 'heatorient' in GraphFlags.model_fields
    assert GraphFlags(heatorient='horizontal').heatorient == 'horizontal'


def test_the_stacked_panels_sit_close_together(monkeypatch):
    """The first stacked render left a band of white between the panels taller than the panels
    themselves: `square=True` ties panel height to width * rows/columns, so with many features
    and few comparisons a fixed figure height is mostly slack. The height is now derived from
    what is being drawn, so the gap has to stay small relative to the panels."""
    captured = _drawn(monkeypatch, 'horizontal')
    top, bottom = (ax.get_position() for ax in captured['axes'])

    gap = top.y0 - bottom.y1
    panel_height = top.y1 - top.y0
    assert gap < panel_height * 2, f'gap {gap:.3f} against panel height {panel_height:.3f}'


def test_the_stacked_colorbars_are_not_drawn_under_the_panels(monkeypatch):
    """A horizontal colorbar below each panel was the first attempt; the bottom panel's landed
    on top of the rotated feature labels."""
    captured = _drawn(monkeypatch, 'horizontal')

    assert all('orientation' not in (kws or {}) for kws in captured['cbar_kws'])
