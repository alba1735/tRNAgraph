"""Precedence and output-format handling for `--style`.

The rule the user set: built-in defaults, then the style file's `defaults`, then the graph
type's own block, then any CLI flag. A flag left unset must never mask the file, which is
why `None` overrides are dropped rather than applied.

`save_figure` centralises what 35 call sites used to hardcode as '.pdf', so `--format svg`
reaches every plot without each module knowing the flag exists.
"""
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from trnagraph.modules import toolsTG
from trnagraph.modules.toolsSchemas import OUTPUT_FORMATS, StyleFile


def _style(payload):
    return StyleFile.model_validate(payload)


# --- precedence ---------------------------------------------------------------------

def test_builtin_defaults_apply_with_no_style_file():
    resolved = toolsTG.resolve_plot_style(None, 'volcano')

    assert resolved['format'] == 'pdf'
    assert resolved['dpi'] == 300


def test_file_defaults_beat_builtins():
    resolved = toolsTG.resolve_plot_style(_style({'defaults': {'dpi': 150}}), 'volcano')

    assert resolved['dpi'] == 150


def test_graph_type_block_beats_file_defaults():
    style = _style({'defaults': {'marker_size': 10}, 'volcano': {'marker_size': 40}})

    assert toolsTG.resolve_plot_style(style, 'volcano')['marker_size'] == 40
    assert toolsTG.resolve_plot_style(style, 'pca')['marker_size'] == 10


def test_cli_flag_beats_everything():
    style = _style({'defaults': {'format': 'png'}, 'volcano': {'format': 'svg'}})

    resolved = toolsTG.resolve_plot_style(style, 'volcano', format='pdf')

    assert resolved['format'] == 'pdf'


def test_an_unset_flag_does_not_mask_the_file():
    """The bug this guards: a flag defaulting to None overwriting a configured value."""
    style = _style({'defaults': {'format': 'svg', 'dpi': 150}})

    resolved = toolsTG.resolve_plot_style(style, 'volcano', format=None, dpi=None)

    assert resolved['format'] == 'svg'
    assert resolved['dpi'] == 150


def test_resolution_does_not_mutate_the_builtin_defaults():
    toolsTG.resolve_plot_style(_style({'defaults': {'dpi': 1}}), 'pca')

    assert toolsTG.PLOT_STYLE_DEFAULTS['dpi'] == 300


# --- save_figure --------------------------------------------------------------------

@pytest.mark.parametrize('fmt', OUTPUT_FORMATS)
def test_save_figure_writes_the_configured_format(tmp_path, fmt):
    fig = plt.figure()
    try:
        path = toolsTG.save_figure(fig, str(tmp_path / 'plot'), {'format': fmt, 'dpi': 100})
    finally:
        plt.close(fig)

    assert path.endswith(f'.{fmt}')
    assert (tmp_path / f'plot.{fmt}').exists()


def test_save_figure_defaults_to_pdf(tmp_path):
    fig = plt.figure()
    try:
        path = toolsTG.save_figure(fig, str(tmp_path / 'plot'))
    finally:
        plt.close(fig)

    assert path.endswith('.pdf')


def test_save_figure_takes_the_extension_from_settings_not_the_caller(tmp_path):
    """Call sites pass a stem; appending '.pdf' themselves is what this replaces."""
    fig = plt.figure()
    try:
        path = toolsTG.save_figure(fig, str(tmp_path / 'plot'), {'format': 'png'})
    finally:
        plt.close(fig)

    assert path.endswith('plot.png')
    assert not (tmp_path / 'plot.png.pdf').exists()


def test_save_figure_honours_dpi(tmp_path):
    """A PNG's pixel dimensions are observable, so dpi is checkable rather than assumed."""
    small = plt.figure(figsize=(2, 2))
    large = plt.figure(figsize=(2, 2))
    try:
        toolsTG.save_figure(small, str(tmp_path / 'small'), {'format': 'png', 'dpi': 50},
                            bbox_inches=None)
        toolsTG.save_figure(large, str(tmp_path / 'large'), {'format': 'png', 'dpi': 200},
                            bbox_inches=None)
    finally:
        plt.close(small)
        plt.close(large)

    assert (tmp_path / 'large.png').stat().st_size > (tmp_path / 'small.png').stat().st_size
