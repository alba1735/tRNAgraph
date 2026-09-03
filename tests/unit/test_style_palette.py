"""The `--style` palette surface: gradients, the categorical fallback, and the template.

Three distinct things are pinned here.

*Resolution* -- a style file's gradient/categorical values become drawable colormaps and color
lists, and an unset file falls through to plotsPalette's own constants. These tests assert the
resolver reads those constants; the concrete values the constants hold are pinned separately in
test_plots_palette.py, so retuning a default is a one-line edit in one place rather than a
change that has to be mirrored across two test files.

*The template* -- `assets/style.template.json` is generated from the schema, so a test
re-derives it and fails if the shipped file drifts. It also asserts the property that makes
the template safe to ship: passing it as a real `--style` file changes nothing at all.

*Leaks* -- a gradient role the schema exposes but no plot module reads is a knob that goes
nowhere. That is the inverse of test_plots_palette.py's rule against hex literals inside plot
modules; together the two close the loop in both directions.
"""
import ast
import json
import pathlib
import types

import numpy as np
import pytest
import seaborn as sns
from pydantic import ValidationError

from trnagraph.modules import plotsPalette, toolsTG, toolsTemplate
from trnagraph.modules.toolsSchemas import (
    GRADIENT_ROLES, GRAPH_STYLE_SUPPORT, STYLE_GRAPH_TYPES, StyleBlock, StyleFile,
    UNIVERSAL_STYLE_KEYS)

TEMPLATE_PATH = pathlib.Path('src/trnagraph/assets/style.template.json')
PLOT_MODULES = sorted(pathlib.Path('src/trnagraph/modules').glob('plots*.py'))


# --- resolution -----------------------------------------------------------------------

def test_unset_gradients_resolve_to_the_modules_declared_defaults():
    """
    The resolver must read plotsPalette's own constants, so that retuning a default is a
    one-line edit there rather than something that also has to be mirrored here. The concrete
    VALUES those constants hold are pinned separately, in test_plots_palette.py.
    """
    resolved = plotsPalette.resolve_gradients(None)
    declared = {
        'correlation': plotsPalette.SEQUENTIAL_CORRELATION,
        'significance': plotsPalette.SEQUENTIAL_SIGNIFICANCE,
        'score': plotsPalette.SEQUENTIAL_SCORE,
        'sequence': plotsPalette.SEQUENTIAL_SEQUENCE,
        'ordered': plotsPalette.SEQUENTIAL_ORDERED,
        'lfc': plotsPalette.DIVERGING_LFC,
    }
    for role, name in declared.items():
        assert resolved[role].name == name


def test_a_named_gradient_is_resolved():
    style = StyleFile.model_validate({'gradients': {'correlation': 'crest'}})

    assert plotsPalette.resolve_gradients(style.gradients)['correlation'].name == 'crest'


def test_a_color_list_becomes_an_interpolated_ramp():
    style = StyleFile.model_validate({'gradients': {'score': ['#000000', '#ffffff']}})
    cmap = plotsPalette.resolve_gradients(style.gradients)['score']

    assert cmap(0.0)[:3] == (0.0, 0.0, 0.0)
    assert cmap(1.0)[:3] == (1.0, 1.0, 1.0)
    # Interpolated, not banded.
    assert 0.4 < cmap(0.5)[0] < 0.6


def test_overriding_one_role_leaves_the_others_at_their_defaults():
    style = StyleFile.model_validate({'gradients': {'lfc': 'icefire'}})
    resolved = plotsPalette.resolve_gradients(style.gradients)

    assert resolved['lfc'].name == 'icefire'
    assert resolved['significance'].name == plotsPalette.SEQUENTIAL_SIGNIFICANCE


def test_every_graph_type_receives_every_role():
    """A per-type subset would need a role->graph-type map, maintained by hand as plots grow."""
    for graph_type in STYLE_GRAPH_TYPES:
        settings = toolsTG.resolve_plot_style(None, graph_type)
        assert set(settings['gradients']) == set(GRADIENT_ROLES)


def test_gradient_falls_back_when_called_without_settings():
    """Plot modules must stay usable from tests and from paths that resolved no style."""
    assert plotsPalette.gradient(None, 'ordered').name == plotsPalette.SEQUENTIAL_ORDERED
    assert plotsPalette.gradient({}, 'correlation').name == plotsPalette.SEQUENTIAL_CORRELATION


def test_an_unknown_role_is_a_loud_error():
    with pytest.raises(KeyError, match='Unknown gradient role'):
        plotsPalette.gradient({}, 'significants')


def test_resolved_settings_survive_a_pool_boundary():
    """Coverage and correlation plots are reached through multiprocessing pools."""
    import pickle

    settings = toolsTG.resolve_plot_style(
        StyleFile.model_validate({'gradients': {'ordered': ['#000000', '#ffffff']}}), 'coverage')
    restored = pickle.loads(pickle.dumps(settings))

    assert restored['gradients']['ordered'](0.5) == settings['gradients']['ordered'](0.5)


# --- discrete sampling ----------------------------------------------------------------

@pytest.mark.parametrize('name', ['mako_r', 'crest', 'Blues'])
@pytest.mark.parametrize('n', [2, 3, 4, 8])
def test_discrete_colors_matches_seaborn_exactly(name, n):
    """
    plotsCoverage takes N discrete colors from the `ordered` ramp. seaborn samples a continuous
    colormap at linspace(0, 1, n + 2)[1:-1], trimming both ends; sampling at linspace(0, 1, n)
    instead shifts colors by up to 0.4 per channel. Switching that call site from a colormap
    NAME to a resolved Colormap object therefore has to reproduce seaborn's spacing, or the
    coverage partition silently changes appearance for users who set nothing.
    """
    sampled = plotsPalette.discrete_colors(plotsPalette.build_colormap(name), n)

    assert np.allclose(sampled, sns.color_palette(name, n))


def test_naive_sampling_would_have_drifted():
    """Guard the guard: prove the spacing correction is load-bearing, not cargo-culted."""
    cmap = plotsPalette.build_colormap('mako_r')
    naive = cmap(np.linspace(0, 1, 5))[:, :3]

    assert not np.allclose(naive, sns.color_palette('mako_r', 5))


# --- categorical fallback -------------------------------------------------------------

def test_unset_categorical_uses_the_built_in_fallback():
    """Values are pinned in test_plots_palette.py; this pins that the plumbing reaches them."""
    assert plotsPalette.categorical({}, 6) == plotsPalette.categorical_palette(6)
    assert plotsPalette.categorical(None, 3) == plotsPalette.categorical_palette(3)


def test_a_named_categorical_palette_is_used():
    settings = {'categorical': 'colorblind'}

    assert np.allclose(plotsPalette.categorical(settings, 4), sns.color_palette('colorblind', 4))


def test_an_explicit_list_is_used_in_order():
    settings = {'categorical': ['#ff0000', '#00ff00', '#0000ff']}

    assert plotsPalette.categorical(settings, 3) == ['#ff0000', '#00ff00', '#0000ff']


def test_a_short_list_cycles_rather_than_being_replaced(caplog):
    """
    A user who wrote a list has said which colors they want; swapping in generated hues they
    never chose would defeat supplying it. The built-in default falls back to husl instead --
    that asymmetry is deliberate, and test_the_builtin_fallback_never_repeats_a_color pins it.
    """
    settings = {'categorical': ['#ff0000', '#00ff00']}

    assert plotsPalette.categorical(settings, 5) == [
        '#ff0000', '#00ff00', '#ff0000', '#00ff00', '#ff0000']
    assert 'being cycled' in caplog.text


def test_the_builtin_fallback_never_repeats_a_color():
    """
    Below the cap the colorblind-safe set is sliced; above it husl generates as many distinct
    hues as asked for. Neither path may repeat, which is what separates the built-in fallback
    from an explicit user list (that one is deliberately cycled).
    """
    for n in (5, len(plotsPalette.CATEGORICAL_COLORS), 22):
        assert len({tuple(c) if not isinstance(c, str) else c
                    for c in plotsPalette.categorical({}, n)}) == n


# --- schema ---------------------------------------------------------------------------

@pytest.mark.parametrize('payload,message', [
    ({'gradients': {'correlation': 'nosuchmap'}}, 'not registered'),
    ({'gradients': {'significants': 'Blues'}}, 'Extra inputs'),
    ({'gradients': {'lfc': ['#ffffff']}}, 'at least two colors'),
    ({'gradients': {'score': ['#zzzzzz', '#ffffff']}}, 'does not recognize as colors'),
    ({'categorical': 'nosuchpalette'}, 'does not recognize'),
    ({'categorical': []}, 'at least one color'),
])
def test_bad_palette_values_fail_at_load_not_at_draw(payload, message):
    """A style file is validated once; a bad value must not survive to a pooled worker."""
    with pytest.raises(ValidationError, match=message):
        StyleFile.model_validate(payload)


def test_seaborn_registered_colormaps_are_accepted():
    """mako/vlag/crest are registered by seaborn, not shipped with matplotlib."""
    style = StyleFile.model_validate(
        {'gradients': {'ordered': 'mako_r', 'lfc': 'vlag', 'correlation': 'crest'}})

    assert style.gradients.ordered == 'mako_r'


def test_a_legacy_colormap_file_still_loads_alongside_the_new_keys():
    """The bare-colormap shape is detected by absence of known keys, which now include these."""
    legacy = StyleFile.model_validate({'group': {'VC_log': '#44B7C2'}})

    assert legacy.colors == {'group': {'VC_log': '#44B7C2'}}
    assert legacy.gradients is None


def test_matplotlib_color_names_are_accepted_as_stops():
    style = StyleFile.model_validate({'gradients': {'sequence': ['tab:blue', 'white', 'red']}})

    assert plotsPalette.resolve_gradients(style.gradients)['sequence'] is not None


# --- the shipped template -------------------------------------------------------------

def _template():
    return json.loads(TEMPLATE_PATH.read_text())


def test_the_template_validates():
    StyleFile.model_validate(_template())


def test_the_template_enumerates_exactly_the_schema_top_level():
    expected = {'colors', 'gradients', 'categorical', 'defaults', *STYLE_GRAPH_TYPES}

    assert set(_template()) == expected


def test_the_template_enumerates_exactly_the_gradient_roles():
    assert set(_template()['gradients']) == set(GRADIENT_ROLES)


def test_the_template_defaults_block_enumerates_every_setting():
    assert set(_template()['defaults']) == set(StyleBlock.model_fields)


@pytest.mark.parametrize('graph_type', STYLE_GRAPH_TYPES)
def test_each_template_graph_block_lists_exactly_what_that_type_supports(graph_type):
    """
    A template that offered `alpha` under `coverage` would advertise a key the schema rejects,
    and one that omitted a supported key would hide it. Both directions are checked.
    """
    expected = UNIVERSAL_STYLE_KEYS | GRAPH_STYLE_SUPPORT.get(graph_type, frozenset())

    assert set(_template()[graph_type]) == expected


def test_every_template_value_is_null():
    template = _template()
    for key, value in template.items():
        if isinstance(value, dict):
            assert all(v is None for v in value.values()), f'{key} carries a non-null value'
        else:
            assert value is None, f'{key} carries a non-null value'


@pytest.mark.parametrize('graph_type', STYLE_GRAPH_TYPES)
def test_the_template_is_a_no_op(graph_type):
    """
    The property that makes an all-null template safe to ship: passing it as a real --style
    file resolves to exactly what passing no style file at all resolves to, so nothing in it
    is accidentally load-bearing.
    """
    blank = toolsTG.resolve_plot_style(StyleFile.model_validate(_template()), graph_type)
    unset = toolsTG.resolve_plot_style(None, graph_type)

    assert blank['gradients'].keys() == unset['gradients'].keys()
    for role in unset['gradients']:
        assert blank['gradients'][role](0.5) == unset['gradients'][role](0.5)
    def strip(s):
        return {k: v for k, v in s.items() if k != 'gradients'}
    assert strip(blank) == strip(unset)


# --- leaks ----------------------------------------------------------------------------

def _roles_read_by_plot_modules():
    """Every string constant passed to a plotsPalette.gradient(...) call, across the modules."""
    found = set()
    for path in PLOT_MODULES:
        for node in ast.walk(ast.parse(path.read_text())):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            if isinstance(func, ast.Attribute) and func.attr == 'gradient' and len(node.args) == 2:
                role = node.args[1]
                if isinstance(role, ast.Constant):
                    found.add(role.value)
    return found


def test_every_exposed_gradient_role_is_actually_drawn_with():
    """
    A role the style file accepts but no plot reads is a knob that silently does nothing --
    the failure mode the whole schema exists to prevent, pointed at itself.
    """
    unread = set(GRADIENT_ROLES) - _roles_read_by_plot_modules()

    assert not unread, (
        f'{sorted(unread)} can be set in a style file but no plot module reads them. Either '
        f'wire the role into the plot that should honour it, or remove it from GRADIENT_ROLES '
        f'and the template.'
    )


def test_no_plot_module_reads_a_role_the_schema_does_not_expose():
    undeclared = _roles_read_by_plot_modules() - set(GRADIENT_ROLES)

    assert not undeclared, (
        f'{sorted(undeclared)} are read by plot modules but are not declared in '
        f'GRADIENT_ROLES, so a style file cannot set them and gradient() will raise.'
    )


def _color_palette_calls_with_non_literal_first_arg():
    """Every `sns.color_palette(X, ...)` in a plot module where X is not a string literal."""
    offenders = []
    # plotsPalette is where a name legitimately arrives in a variable: it resolves the
    # built-in overflow palette and whatever name a --style file supplied. Everywhere else a
    # variable in that position is a resolved Colormap and therefore a bug.
    for path in (p for p in PLOT_MODULES if p.name != 'plotsPalette.py'):
        for node in ast.walk(ast.parse(path.read_text())):
            if not isinstance(node, ast.Call) or not node.args:
                continue
            func = node.func
            name = func.attr if isinstance(func, ast.Attribute) else getattr(func, 'id', None)
            if name != 'color_palette':
                continue
            first = node.args[0]
            if not (isinstance(first, ast.Constant) and isinstance(first.value, str)):
                offenders.append(f'{path.name}:{node.lineno}')
    return offenders


def test_no_plot_module_passes_a_resolved_colormap_to_color_palette():
    """
    `sns.color_palette` takes a NAME; a resolved gradient is a Colormap OBJECT, and passing one
    raises `TypeError: 'ListedColormap' object is not iterable`.

    This is not hypothetical. plotsCluster kept `numericcolormap` for two consumers at once --
    `ScalarMappable(cmap=...)`, which needs the object, and `color_palette`, which needs the
    name -- so resolving the role to an object fixed one call site and broke the other. It
    escaped every render check because the demo object carries no `cluster_runinfo`, so
    `-g cluster` is skipped on it entirely and only a real clustered object reaches the code.
    A static check costs nothing and does not depend on having the right data to hand.

    Take discrete colors from a resolved gradient with plotsPalette.discrete_colors() instead.
    """
    offenders = _color_palette_calls_with_non_literal_first_arg()

    assert not offenders, (
        f'{offenders} pass a non-literal to sns.color_palette(). If that value is a resolved '
        f'gradient it will raise at draw time; use plotsPalette.discrete_colors(cmap, n).'
    )


def test_the_check_would_have_caught_the_original_bug(tmp_path):
    """Guard the guard: the exact shape that shipped must be rejected by this detector."""
    import ast as _ast

    bad = _ast.parse('pal = sns.color_palette(self.numericcolormap, n)')
    found = [n for n in _ast.walk(bad)
             if isinstance(n, _ast.Call) and getattr(n.func, 'attr', None) == 'color_palette'
             and not isinstance(n.args[0], _ast.Constant)]

    assert found, 'the detector no longer recognises the shape it was written for'


@pytest.mark.parametrize('role', GRADIENT_ROLES)
def test_a_resolved_gradient_works_with_both_of_its_consumers(role):
    """
    A resolved role is handed to two different APIs across the plot modules: matplotlib's
    cmap= parameter (which wants the object) and discrete sampling (which does not). Both must
    accept what resolve_gradients returns, for every role.
    """
    import matplotlib.pyplot as plt

    cmap = plotsPalette.resolve_gradients(None)[role]

    assert plt.cm.ScalarMappable(cmap=cmap) is not None
    assert len(plotsPalette.discrete_colors(cmap, 4)) == 4


# --- tools template -------------------------------------------------------------------

def test_the_writer_emits_the_template(tmp_path):
    args = types.SimpleNamespace(style=True, output=str(tmp_path), overwrite=False)
    written = toolsTemplate.TemplateWriter(args).run()

    assert [pathlib.Path(p).name for p in written] == ['style.template.json']
    assert json.loads(pathlib.Path(written[0]).read_text()) == _template()


def test_no_selector_writes_everything(tmp_path):
    """A user who does not know the file names cannot be expected to name one."""
    args = types.SimpleNamespace(style=False, output=str(tmp_path), overwrite=False)

    assert len(toolsTemplate.TemplateWriter(args).run()) == len(toolsTemplate.TEMPLATES)


def test_an_existing_file_is_not_clobbered(tmp_path):
    """Emit, edit, re-emit after an upgrade is the workflow; silent replacement loses work."""
    (tmp_path / 'style.template.json').write_text('{"edited": true}')
    args = types.SimpleNamespace(style=True, output=str(tmp_path), overwrite=False)

    with pytest.raises(FileExistsError, match='--overwrite'):
        toolsTemplate.TemplateWriter(args).run()
    assert (tmp_path / 'style.template.json').read_text() == '{"edited": true}'


def test_overwrite_replaces_it(tmp_path):
    (tmp_path / 'style.template.json').write_text('{"edited": true}')
    args = types.SimpleNamespace(style=True, output=str(tmp_path), overwrite=True)
    toolsTemplate.TemplateWriter(args).run()

    assert json.loads((tmp_path / 'style.template.json').read_text()) == _template()


def test_the_template_resolves_inside_the_installed_package():
    """
    Walking directories up from __file__ landed outside the package under a non-editable
    install; assets_dir() goes through importlib.resources instead.
    """
    assert (pathlib.Path(toolsTG.assets_dir()) / 'style.template.json').is_file()


# --- line width ------------------------------------------------------------------------

def test_line_width_is_accepted_by_the_graph_types_that_draw_lines_or_bar_edges():
    """coverage/radar/mismatch/count are the modules with a data trace or a bar edge. The
    scatter types are served by marker_size, and heatmap/correlation/logo have no stroke a
    user would want to set. `venn` takes alpha instead: its circles have no stroke either, and
    what a reader needs to control there is how far an overlap darkens."""
    for graph_type in ('coverage', 'radar', 'mismatch', 'count'):
        assert 'line_width' in GRAPH_STYLE_SUPPORT[graph_type]


def test_line_width_is_rejected_where_it_would_do_nothing():
    """extra='forbid' plus this table is what stops a style file silently doing nothing --
    a key accepted under `heatmap` would advertise a setting that never reaches a figure."""
    for graph_type in ('heatmap', 'correlation', 'logo'):
        assert 'line_width' not in GRAPH_STYLE_SUPPORT[graph_type]
        with pytest.raises(ValidationError):
            StyleFile.model_validate({graph_type: {'line_width': 1.0}})


def test_line_width_must_be_positive():
    with pytest.raises(ValidationError):
        StyleBlock(line_width=0)
    with pytest.raises(ValidationError):
        StyleBlock(line_width=-1)
    assert StyleBlock(line_width=0.5).line_width == 0.5
