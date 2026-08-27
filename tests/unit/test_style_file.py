"""`--style` replaces `--colormap`, carrying palette and presentation in one file.

Preparing a figure set means reusing one styling file across many filter configs, so styling
is a different concern with a different lifetime from `--config` (which selects observations
and whose `name` becomes an output subdirectory). It is one file rather than two because the
palette and the presentation settings are always chosen together.

`defaults` applies everywhere; a graph-type block overrides it for that type; a CLI flag
overrides both.
"""
import pytest
from pydantic import ValidationError

from trnagraph.modules.toolsSchemas import (
    OUTPUT_FORMATS,
    STYLE_GRAPH_TYPES,
    StyleBlock,
    StyleFile,
)


def test_defaults_apply_to_every_graph_type():
    style = StyleFile.model_validate({'defaults': {'dpi': 200, 'font_size': 8}})

    for graph_type in STYLE_GRAPH_TYPES:
        assert style.resolve(graph_type) == {'dpi': 200, 'font_size': 8.0}


def test_a_graph_type_block_overrides_defaults():
    style = StyleFile.model_validate({
        'defaults': {'marker_size': 10, 'dpi': 300},
        'volcano': {'marker_size': 40},
    })

    assert style.resolve('volcano') == {'marker_size': 40.0, 'dpi': 300}
    assert style.resolve('pca') == {'marker_size': 10.0, 'dpi': 300}


def test_unset_keys_are_absent_rather_than_none():
    """A caller must be able to tell 'not configured' from 'configured to the built-in value'."""
    style = StyleFile.model_validate({'defaults': {'dpi': 300}})

    assert 'marker_size' not in style.resolve('pca')


def test_an_empty_file_resolves_to_nothing():
    style = StyleFile.model_validate({})

    assert style.resolve('volcano') == {}
    assert style.colors_for('group') is None


# --- legacy colormap compatibility --------------------------------------------------

def test_a_legacy_colormap_file_is_read_as_colors():
    """The shipped assets/colormap.json shape must keep working under --style."""
    legacy = {'group': {'VC_24h': '#FFAE49', 'VC_log': '#44B7C2'},
              'trimtype': {'Merged': '#44B7C2'}}

    style = StyleFile.model_validate(legacy)

    assert style.colors == legacy
    assert style.colors_for('group')['VC_24h'] == '#FFAE49'


def test_an_explicit_colors_block_is_preferred_over_the_legacy_shape():
    style = StyleFile.model_validate({'colors': {'group': {'A': '#000'}},
                                      'defaults': {'dpi': 100}})

    assert style.colors_for('group') == {'A': '#000'}
    assert style.resolve('pca') == {'dpi': 100}


def test_colors_for_an_unconfigured_column_is_none():
    style = StyleFile.model_validate({'colors': {'group': {'A': '#000'}}})

    assert style.colors_for('amino') is None
    assert style.colors_for(None) is None


# --- validation ---------------------------------------------------------------------

def test_a_misspelled_setting_fails_loudly():
    """Silently ignoring it would leave the user wondering why the figure did not change."""
    with pytest.raises(ValidationError, match='marker_sizes'):
        StyleFile.model_validate({'defaults': {'marker_sizes': 10}})


def test_an_unknown_graph_type_block_fails():
    with pytest.raises(ValidationError):
        StyleFile.model_validate({'defaults': {'dpi': 100}, 'volcanoes': {'marker_size': 5}})


@pytest.mark.parametrize('bad', [
    {'dpi': 0}, {'dpi': -100}, {'marker_size': 0}, {'font_size': -1},
    {'alpha': 1.5}, {'alpha': -0.1}, {'figsize': [0, 5]}, {'figsize': [8, -2]},
    {'format': 'tiff'},
])
def test_out_of_range_values_are_rejected(bad):
    with pytest.raises(ValidationError):
        StyleBlock.model_validate(bad)


@pytest.mark.parametrize('fmt', OUTPUT_FORMATS)
def test_every_supported_format_validates(fmt):
    assert StyleBlock.model_validate({'format': fmt}).format == fmt


def test_style_graph_types_match_the_grapher_dispatch():
    """The style file must not accept a block for a graph type that cannot consume it."""
    from trnagraph.modules.adataGraph import anndataGrapher
    import inspect

    source = inspect.getsource(anndataGrapher.dispatch_plot)
    dispatched = {gt for gt in STYLE_GRAPH_TYPES if f"gt == '{gt}'" in source}
    # 'trimming' is dispatched by preprocess trim, not by the grapher.
    assert dispatched == set(STYLE_GRAPH_TYPES) - {'trimming'}


def test_figsize_accepts_a_two_element_sequence():
    block = StyleBlock.model_validate({'figsize': [8, 6]})

    assert block.figsize == (8.0, 6.0)


# --- per-graph-type key support -----------------------------------------------------

def test_a_key_a_graph_type_cannot_use_is_rejected():
    """Silently ignoring it is the failure mode this schema exists to prevent."""
    with pytest.raises(ValidationError, match='marker_size'):
        StyleFile.model_validate({'count': {'marker_size': 10}})


def test_defaults_may_carry_keys_only_some_types_use():
    """defaults broadcasts, so a type that cannot use a key simply skips it."""
    style = StyleFile.model_validate({'defaults': {'marker_size': 10}})

    assert style.resolve('count') == {'marker_size': 10.0}
    assert style.resolve('volcano') == {'marker_size': 10.0}


def test_universal_keys_are_accepted_by_every_type():
    from trnagraph.modules.toolsSchemas import UNIVERSAL_STYLE_KEYS

    payload = {'dpi': 200, 'font_size': 9, 'figsize': [4, 4], 'format': 'png'}
    assert set(payload) == set(UNIVERSAL_STYLE_KEYS)
    for graph_type in STYLE_GRAPH_TYPES:
        StyleFile.model_validate({graph_type: payload})


def test_structural_alpha_plots_reject_alpha():
    """coverage/radar/logo alpha is shading geometry, not a point-opacity knob."""
    for graph_type in ('coverage', 'radar', 'logo'):
        with pytest.raises(ValidationError, match='alpha'):
            StyleFile.model_validate({graph_type: {'alpha': 0.5}})


def test_scatter_plots_accept_marker_size_and_rasterization():
    for graph_type in ('volcano', 'pca', 'cluster', 'mismatch'):
        style = StyleFile.model_validate(
            {graph_type: {'marker_size': 30, 'alpha': 0.7, 'rasterize_over': 1000}})
        assert style.resolve(graph_type)['marker_size'] == 30
