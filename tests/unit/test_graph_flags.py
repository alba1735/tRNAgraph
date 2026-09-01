"""`graph` options pinned from a --config file.

A config file already names a run and scopes which data it sees, so "which grouping column,
which readtypes, which cutoff" is the same kind of statement -- one file can carry a whole
saved analysis instead of a shell line that has to be retyped correctly each time. These
tests pin the three things that can quietly go wrong.

*Drift* -- GraphFlags is written out by hand, because importing cli.py into toolsSchemas
would invert the package's dependency direction and defeat lazy_imports. The cost of that
choice is that the model can fall out of step with the command, so it is checked against the
real click signature here, and against the shipped template.

*Precedence* -- a flag typed on the command line must beat the file. Nearly every graph option
has a real default rather than a None sentinel, so this cannot be done by comparing values;
cli.py captures click's ParameterSource instead and the merge consults it.

*Ordering* -- almost everything in anndataGrapher.__init__ is derived from options this block
can set (the variant view, the read basis, covtype, the resolved grouping columns, the
graphtypes expansion). The config used to be read after all of them, so the merge had to move
ahead of that work; a test pins the order rather than trusting it to survive edits.
"""
import ast
import inspect
import json
import pathlib
import types

import pytest
import typer
from pydantic import ValidationError

from trnagraph import cli
from trnagraph.modules import adataGraph, toolsTemplate
from trnagraph.modules.toolsSchemas import (
    GRAPH_FLAG_EXCLUSIONS, GRAPH_LIST_FLAGS, GraphFilterConfig, GraphFlags)

TEMPLATE_PATH = pathlib.Path('src/trnagraph/assets/config.template.json')


def _graph_options():
    """Real CLI options on `graph` -- excludes the injected click Context, which is not one."""
    return {name: p for name, p in inspect.signature(cli.graph).parameters.items()
            if p.annotation is not typer.Context}


def _merge(flags, typed=None, **args):
    """Run the merge in isolation, without constructing a real grapher (which needs an h5ad)."""
    grapher = object.__new__(adataGraph.anndataGrapher)
    grapher.logger = __import__('logging').getLogger('test')
    grapher.config = GraphFilterConfig(name='t', flags=GraphFlags(**flags))
    grapher.args = types.SimpleNamespace(cli_specified=typed, **args)
    grapher._apply_config_flags()
    return grapher.args


# --- drift ----------------------------------------------------------------------------

def test_the_model_covers_exactly_the_eligible_graph_options():
    """The hand-written model's only real risk is falling out of step with the command."""
    eligible = set(_graph_options()) - GRAPH_FLAG_EXCLUSIONS

    assert set(GraphFlags.model_fields) == eligible


def test_every_exclusion_is_a_real_graph_option():
    """A stale exclusion would silently widen the model without anyone noticing."""
    assert GRAPH_FLAG_EXCLUSIONS <= set(_graph_options())


def test_paths_and_process_controls_are_not_settable():
    """A config that redirects its own output or reopens itself is a footgun, not a feature."""
    for name in ('anndata', 'output', 'config', 'style', 'threads', 'quiet', 'verbose'):
        assert name not in GraphFlags.model_fields


def test_format_belongs_to_style_not_config():
    """The one key both files could have claimed; excluding it leaves no precedence rule."""
    assert 'format' in GRAPH_FLAG_EXCLUSIONS
    with pytest.raises(ValidationError):
        GraphFlags(format='png')


def test_an_unknown_flag_is_rejected():
    with pytest.raises(ValidationError, match='Extra inputs'):
        GraphFilterConfig.model_validate({'name': 'x', 'flags': {'heatcuttoff': 80}})


def test_every_declared_list_flag_really_is_a_list_field():
    for name in GRAPH_LIST_FLAGS:
        annotation = str(GraphFlags.model_fields[name].annotation)
        assert 'List[str]' in annotation, f'{name} is declared a list flag but is {annotation}'


# --- values ---------------------------------------------------------------------------

@pytest.mark.parametrize('name', GRAPH_LIST_FLAGS)
def test_an_empty_list_is_rejected(name):
    """A list replaces the default, so an empty one selects nothing and reports success."""
    with pytest.raises(ValidationError, match='replaces the default'):
        GraphFlags(**{name: []})


def test_a_list_replaces_rather_than_extends():
    """Replacing is the only semantics that lets a config NARROW a list, which is the point."""
    args = _merge({'diffrts': ['fiveprime']},
                  diffrts=['wholecounts', 'fiveprime', 'threeprime', 'other', 'total'])

    assert args.diffrts == ['fiveprime']


def test_unset_flags_are_left_alone():
    """model_dump(exclude_none=True) is what keeps an unmentioned key from being nulled out."""
    args = _merge({'heatcutoff': 200}, heatcutoff=80, covgrp='group')

    assert (args.heatcutoff, args.covgrp) == (200, 'group')


def test_a_false_boolean_is_applied_not_treated_as_unset():
    """False is a real value; only None means 'not configured'."""
    args = _merge({'lfcshrink': False}, lfcshrink=True)

    assert args.lfcshrink is False


# --- precedence -----------------------------------------------------------------------

def test_a_typed_flag_beats_the_file():
    args = _merge({'heatcutoff': 200}, typed=frozenset({'heatcutoff'}), heatcutoff=999)

    assert args.heatcutoff == 999


def test_an_untyped_flag_takes_the_file_value():
    args = _merge({'heatcutoff': 200}, typed=frozenset({'covgrp'}), heatcutoff=80)

    assert args.heatcutoff == 200


def test_the_file_applies_in_full_when_nothing_was_typed():
    """Namespaces built directly -- the Python API and tests -- carry no cli_specified."""
    grapher = object.__new__(adataGraph.anndataGrapher)
    grapher.logger = __import__('logging').getLogger('test')
    grapher.config = GraphFilterConfig(name='t', flags=GraphFlags(volcutoff=5))
    grapher.args = types.SimpleNamespace(volcutoff=80)
    grapher._apply_config_flags()

    assert grapher.args.volcutoff == 5


def test_no_flags_block_changes_nothing():
    grapher = object.__new__(adataGraph.anndataGrapher)
    grapher.logger = __import__('logging').getLogger('test')
    grapher.config = GraphFilterConfig(name='t')
    grapher.args = types.SimpleNamespace(heatcutoff=80)
    grapher._apply_config_flags()

    assert grapher.args.heatcutoff == 80


def test_cli_specified_params_reports_only_typed_options():
    """Comparing against defaults cannot work: heatcutoff=80 is both a default and typeable."""
    seen = {}
    app = typer.Typer()

    @app.command()
    def probe(ctx: typer.Context,
              alpha: int = typer.Option(80, '--alpha'),
              beta: str = typer.Option('group', '--beta')):
        seen.update(typed=cli.cli_specified_params(ctx))

    from typer.testing import CliRunner
    # --alpha is passed its own default value, which must still count as explicitly typed.
    CliRunner().invoke(app, ['--alpha', '80'])

    assert 'alpha' in seen['typed'] and 'beta' not in seen['typed']


# --- ordering -------------------------------------------------------------------------

def test_flags_are_merged_before_anything_derived_from_them():
    """
    The variant view, read basis, covtype, resolved grouping columns and the graphtypes
    expansion are all computed from options this block can set. The config used to be read
    after all of them, so a `flags` block set there would have been silently too late.
    """
    import textwrap

    # getsource on a method returns it still indented under its class.
    tree = ast.parse(textwrap.dedent(inspect.getsource(adataGraph.anndataGrapher.__init__)))

    order = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            name = getattr(node.func, 'attr', None)
            if name in ('_apply_config_flags', 'parse_variant', 'read_basis',
                        'resolve_covtype', '_validate_label_args'):
                order.append((node.lineno, name))
    order.sort()
    names = [n for _, n in order]

    assert names[0] == '_apply_config_flags', (
        f'_apply_config_flags must run before {names[:3]}; got order {names}')


# --- the shipped template -------------------------------------------------------------

def _template():
    return json.loads(TEMPLATE_PATH.read_text())


def test_the_config_template_validates():
    GraphFilterConfig.model_validate(_template())


def test_the_config_template_enumerates_every_flag():
    assert set(_template()['flags']) == set(GraphFlags.model_fields)


def test_the_config_template_enumerates_every_filter_key():
    assert set(_template()) == {'name', 'obs', 'obs_r', 'var', 'var_r', 'flags'}


def test_every_config_template_flag_is_null():
    assert all(v is None for v in _template()['flags'].values())


def test_the_config_template_sets_nothing():
    """
    `name` is required and becomes an output subdirectory, so unlike the style template this
    one cannot be entirely null. Everything that is optional is, and the flags block resolves
    to nothing at all.
    """
    parsed = GraphFilterConfig.model_validate(_template())

    assert parsed.flags.model_dump(exclude_none=True) == {}
    assert (parsed.obs, parsed.obs_r, parsed.var, parsed.var_r) == (None, None, None, None)


def test_both_templates_are_emitted_by_default(tmp_path):
    args = types.SimpleNamespace(style=False, config=False, output=str(tmp_path), overwrite=False)
    written = {pathlib.Path(p).name for p in toolsTemplate.TemplateWriter(args).run()}

    assert written == {'style.template.json', 'config.template.json'}


def test_the_config_template_can_be_written_alone(tmp_path):
    args = types.SimpleNamespace(style=False, config=True, output=str(tmp_path), overwrite=False)
    written = toolsTemplate.TemplateWriter(args).run()

    assert [pathlib.Path(p).name for p in written] == ['config.template.json']
