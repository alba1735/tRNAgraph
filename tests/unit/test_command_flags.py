"""The `--config` file's per-command `flags` blocks.

`--config` used to carry a single bare `flags` block, which was `graph`'s options and only
`graph`'s. It now takes one block per command -- `flags.graph`, `flags.build`, `flags.map`,
and so on -- so a single file can drive a whole run from trimming through graphing instead of
every stage's settings living in a shell script that has to be retyped correctly.

The models are written out by hand rather than derived from the typer commands, because
importing cli.py into toolsSchemas would invert the package's dependency direction
(cli -> modules) and defeat lazy_imports by dragging the whole CLI in whenever a schema is
touched. These drift tests are what makes that safe: they assert each model's fields are
exactly its command's eligible options, so a flag added to the CLI and forgotten here fails.
"""
import inspect

import pytest
import typer
from pydantic import ValidationError

from trnagraph import cli
from trnagraph.modules.toolsSchemas import (
    COMMAND_FLAG_EXCLUSIONS, COMMAND_FLAG_MODELS, CommandFlags, RunConfig,
)


def _command_callbacks():
    """{command name: its typer callback}, for every command a config block exists for."""
    found = {'graph': cli.graph}
    for app in (cli.preprocess_app, cli.analyze_app, cli.tools_app):
        for info in app.registered_commands:
            found[info.name] = info.callback
    return found


def _eligible(command):
    """The options of one command a config file may set."""
    params = inspect.signature(_command_callbacks()[command]).parameters
    return {name for name, p in params.items()
            if p.annotation is not typer.Context} - COMMAND_FLAG_EXCLUSIONS


@pytest.mark.parametrize('command', sorted(COMMAND_FLAG_MODELS))
def test_each_model_covers_exactly_its_commands_eligible_options(command):
    assert set(COMMAND_FLAG_MODELS[command].model_fields) == _eligible(command)


@pytest.mark.parametrize('command', sorted(COMMAND_FLAG_MODELS))
def test_no_block_can_set_an_output_path_or_a_process_control(command):
    """Input paths ARE settable -- for build, map and makedb they are the analysis, and
    withholding them would leave those blocks nearly empty. What stays excluded is where the
    run WRITES and how it is driven: an output directory, thread count, verbosity, and the
    config and style files themselves, which a file that could name them could redirect."""
    fields = set(COMMAND_FLAG_MODELS[command].model_fields)
    for name in ('output', 'threads', 'quiet', 'config', 'style'):
        assert name not in fields


def test_the_blocks_are_exactly_the_commands_worth_saving():
    """csv, merge, info, test and template are deliberately absent: their options are paths,
    output destinations and one-shot selectors, with nothing worth carrying between runs."""
    assert set(CommandFlags.model_fields) == {
        'graph', 'build', 'map', 'trim', 'makedb', 'cluster', 'addsplit', 'log2fc',
    }
    assert set(COMMAND_FLAG_MODELS) == set(CommandFlags.model_fields)


def test_every_block_name_is_a_real_command():
    assert set(CommandFlags.model_fields) <= set(_command_callbacks())


def test_an_unknown_command_block_is_rejected():
    with pytest.raises(ValidationError, match='Extra inputs'):
        RunConfig.model_validate({'name': 'x', 'flags': {'grpah': {}}})


def test_an_unknown_key_inside_a_block_is_rejected():
    with pytest.raises(ValidationError, match='Extra inputs'):
        RunConfig.model_validate({'name': 'x', 'flags': {'build': {'maxmismatchs': '1'}}})


def test_the_filters_stay_at_the_top_level():
    """obs/var filters are not CLI flags and are not scoped to one command: they say which
    subset of the object every AnnData-consuming command sees."""
    config = RunConfig.model_validate({
        'name': 'x', 'obs_r': {'amino': ['Und']}, 'flags': {'graph': {'heatcutoff': 20}},
    })
    assert config.obs_r == {'amino': ['Und']}
    assert config.flags.graph.heatcutoff == 20


def _rejection_lines(payload):
    from trnagraph.modules.toolsSchemas import explain_rejected_keys
    with pytest.raises(ValidationError) as caught:
        RunConfig.model_validate(payload)
    return '\n'.join(explain_rejected_keys(caught.value, 'config'))


def test_an_old_style_bare_flags_block_says_where_its_keys_now_go():
    """`flags` used to BE graph's options. It is now a mapping of command to options, so every
    existing config file has keys one level too shallow -- and pydantic reports only 'Extra
    inputs are not permitted', which is true and useless. Named rather than aliased: the tool
    is pre-publication with effectively one user, and a permanent alias for a key that has
    changed meaning is worse than a clear error."""
    lines = _rejection_lines({'name': 'x', 'flags': {'heatcutoff': 20, 'covgrp': 'timepoint'}})

    assert 'flags.graph' in lines, lines
    assert 'heatcutoff' in lines and 'covgrp' in lines, lines


def test_a_key_in_the_wrong_command_block_is_pointed_at_the_right_one():
    """`gtf` is a build option; putting it under `flags.map` is an easy slip once there are
    eight blocks, and 'did you mean' over map's own keys would answer with nonsense."""
    lines = _rejection_lines({'name': 'x', 'flags': {'map': {'gtf': 'x.gtf'}}})

    assert 'flags.build' in lines, lines


def test_a_style_key_in_a_command_block_still_points_at_the_style_file():
    lines = _rejection_lines({'name': 'x', 'flags': {'graph': {'line_width': 0.5}}})

    assert '--style' in lines, lines


# --- what a config file has to be able to supply -------------------------------------------

def test_options_a_config_can_supply_are_not_required_by_the_parser():
    """`-i`, `-d` and friends used to be typer-required, which made them impossible to set
    from a file: click rejected the invocation before the config was ever read, so a config
    that named a manifest still had to have it retyped on the command line. They are optional
    at the parser now and checked after the merge instead."""
    import typer as _typer

    for command, names in [('trim', ['input']), ('map', ['database', 'input']),
                           ('build', ['input', 'database']),
                           ('makedb', ['genome', 'trnaout', 'trnafa', 'namemap']),
                           ('addsplit', ['readlengthsplit'])]:
        params = inspect.signature(_command_callbacks()[command]).parameters
        for name in names:
            default = params[name].default
            assert default.default is None, (
                f'{command} --{name} is still parser-required, so a config cannot supply it')
            assert default is not _typer.models.ArgumentInfo


def test_a_missing_required_option_names_both_ways_to_supply_it():
    from types import SimpleNamespace

    from trnagraph.modules.toolsTG import InvalidParameterError, require_options

    with pytest.raises(InvalidParameterError) as caught:
        require_options(SimpleNamespace(input=None, database='db'), 'build',
                        ['input', 'database'])

    message = str(caught.value)
    assert '--input' in message and 'flags.build' in message
    assert '--database' not in message, 'only the missing option should be named'


def test_the_typed_merge_bookkeeping_never_reaches_the_h5ad():
    """`cli_specified` rides on the args namespace so the merge knows what the user typed, and
    adataBuild records that namespace into uns['trnagraphruninfo']['flags']. It is a frozenset,
    which h5ad cannot write -- recording it failed the entire build at the moment of saving."""
    import inspect as _inspect

    from trnagraph.modules import adataBuild

    source = _inspect.getsource(adataBuild.AnnDataBuilder.__init__)
    assert "if k != 'cli_specified'" in source
