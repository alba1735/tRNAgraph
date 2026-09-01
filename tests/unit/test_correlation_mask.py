"""`--corrmask` hides one half of the correlation matrix.

A correlation matrix is symmetric, so half of every one restates the other half. Masking a
triangle is the conventional way to say so. The diagonal is kept: it is R^2 = 1 by
construction and carries no information, but it is the visual anchor tying a row to its axis
label, and dropping it turns the triangle into a stair-step.
"""
import numpy as np
import pandas as pd
import pytest

from trnagraph.modules import plotsCorrelation, toolsTG


def _frame(n=4):
    values = np.eye(n) + 0.1
    return pd.DataFrame(values, index=[f'r{i}' for i in range(n)],
                        columns=[f'r{i}' for i in range(n)])


@pytest.mark.parametrize('requested', [None, 'none'])
def test_no_mask_is_the_default(requested):
    assert plotsCorrelation._triangle_mask(_frame(), requested) is None


def test_upper_hides_above_the_diagonal_and_keeps_the_diagonal():
    mask = plotsCorrelation._triangle_mask(_frame(), 'upper')

    assert mask[0, 1], 'a cell above the diagonal is hidden'
    assert not mask[1, 0], 'the mirrored cell below it is kept'
    assert not np.diagonal(mask).any(), 'the diagonal anchors the axis labels and stays'


def test_lower_hides_below_the_diagonal_and_keeps_the_diagonal():
    mask = plotsCorrelation._triangle_mask(_frame(), 'lower')

    assert mask[1, 0]
    assert not mask[0, 1]
    assert not np.diagonal(mask).any()


def test_an_unrecognised_mask_is_a_usage_error_naming_the_choices():
    """Silently drawing the full matrix for a mistyped value would be the same
    silently-wrong-output failure the label validation removed."""
    with pytest.raises(toolsTG.InvalidParameterError) as raised:
        plotsCorrelation._triangle_mask(_frame(), 'top')
    message = str(raised.value)
    assert '--corrmask' in message
    for choice in ('none', 'upper', 'lower'):
        assert choice in message


def _wide(n=4):
    """A wide (feature x sample) frame whose values clear _plot_corr_matrix's >= 20 threshold."""
    import numpy as np
    rows = np.arange(1, n * 6 + 1, dtype=float).reshape(6, n) * 5
    return pd.DataFrame(rows, columns=[f's{i}' for i in range(n)])


def _draw(monkeypatch, corr_mask):
    import matplotlib
    matplotlib.use('Agg')

    captured = {}
    real_heatmap = plotsCorrelation.sns.heatmap

    def spy(*args, **kwargs):
        captured['mask'] = kwargs.get('mask')
        return real_heatmap(*args, **kwargs)

    def fake_save(path, settings=None):
        captured['path'] = path
        return path

    monkeypatch.setattr(plotsCorrelation.sns, 'heatmap', spy)
    monkeypatch.setattr(plotsCorrelation.toolsTG, 'save_current', fake_save)
    plotsCorrelation._plot_corr_matrix(_wide(), 'pearson', 'sample', 'out/',
                                       'pearson_sample_total_correlation_matrix', 'total',
                                       False, corr_mask=corr_mask)
    return captured


def test_the_default_draws_the_full_matrix_under_its_existing_filename(monkeypatch):
    """`none` must be byte-for-byte what earlier versions produced, filename included --
    otherwise adding the option silently relocates everyone's existing output."""
    captured = _draw(monkeypatch, 'none')

    assert captured['mask'] is None
    assert captured['path'] == 'out/pearson_sample_total_correlation_matrix.pdf'


@pytest.mark.parametrize('corr_mask', ['upper', 'lower'])
def test_a_masked_matrix_is_drawn_masked_and_says_so_in_its_filename(monkeypatch, corr_mask):
    """The token exists so a masked and an unmasked run can sit side by side: the two render
    the same data, so without it the second would silently overwrite the first."""
    captured = _draw(monkeypatch, corr_mask)

    assert captured['mask'] is not None
    assert captured['path'] == f'out/pearson_sample_total_correlation_matrix_{corr_mask}.pdf'


def test_the_cli_exposes_corrmask_and_defaults_to_the_full_matrix():
    import inspect

    from trnagraph import cli

    parameter = inspect.signature(cli.graph).parameters['corrmask']
    assert parameter.default.default == 'none', 'the full matrix stays the default'


def test_corrmask_is_settable_from_a_config_file():
    """Every graph option that is a selection choice belongs in --config's flags block, so one
    file can carry a whole saved analysis. A test in test_graph_flags.py pins the model against
    the command; this one states the intent for this flag specifically."""
    from trnagraph.modules.toolsSchemas import GraphFlags

    assert 'corrmask' in GraphFlags.model_fields
    assert GraphFlags(corrmask='upper').corrmask == 'upper'
