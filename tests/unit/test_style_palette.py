"""The `--style` palette surface: gradients, the categorical fallback, and the template.

Three distinct things are pinned here.

*Resolution* -- a style file's gradient/categorical values become drawable colormaps and color
lists, and an UNSET file resolves to exactly what the built-in defaults were. Commit-scope
matters: adding user overrides must not restyle anything for a user who configured nothing,
so the no-override paths are asserted against the previous behavior rather than against
themselves.

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

def test_unset_gradients_resolve_to_the_built_in_defaults():
    """Adding the override machinery must not restyle a user who configured nothing."""
    resolved = plotsPalette.resolve_gradients(None)
    previous = {
        'correlation': 'Blues', 'significance': 'Greens', 'score': 'Blues',
        'sequence': 'Greys', 'ordered': 'mako_r',
    }
    for role, name in previous.items():
        assert resolved[role].name == name


def test_unset_lfc_matches_the_previous_inline_diverging_palette():
    resolved = plotsPalette.resolve_gradients(None)['lfc']
    expected = sns.diverging_palette(as_cmap=True, **plotsPalette.DIVERGING_LFC_KWARGS)

    for point in (0.0, 0.5, 1.0):
        assert resolved(point) == expected(point)


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
    style = StyleFile.model_validate({'gradients': {'lfc': 'vlag'}})
    resolved = plotsPalette.resolve_gradients(style.gradients)

    assert resolved['lfc'].name == 'vlag'
    assert resolved['significance'].name == 'Greens'


def test_every_graph_type_receives_every_role():
    """A per-type subset would need a role->graph-type map, maintained by hand as plots grow."""
    for graph_type in STYLE_GRAPH_TYPES:
        settings = toolsTG.resolve_plot_style(None, graph_type)
        assert set(settings['gradients']) == set(GRADIENT_ROLES)


def test_gradient_falls_back_when_called_without_settings():
    """Plot modules must stay usable from tests and from paths that resolved no style."""
    assert plotsPalette.gradient(None, 'ordered').name == 'mako_r'
    assert plotsPalette.gradient({}, 'correlation').name == 'Blues'


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

def test_unset_categorical_is_unchanged_from_the_built_in_fallback():
    assert np.allclose(plotsPalette.categorical({}, 6), sns.husl_palette(6))
    assert np.allclose(plotsPalette.categorical(None, 3), sns.husl_palette(3))


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
    that asymmetry is deliberate, and test_the_builtin_fallback_is_not_cycled pins the pair.
    """
    settings = {'categorical': ['#ff0000', '#00ff00']}

    assert plotsPalette.categorical(settings, 5) == [
        '#ff0000', '#00ff00', '#ff0000', '#00ff00', '#ff0000']
    assert 'being cycled' in caplog.text


def test_the_builtin_fallback_is_not_cycled():
    """husl generates as many distinct hues as asked for, so nothing repeats."""
    assert len({tuple(c) for c in plotsPalette.categorical({}, 22)}) == 22


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
    strip = lambda s: {k: v for k, v in s.items() if k != 'gradients'}
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
