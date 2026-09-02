"""The CLI help / documentation split.

The convention is that `--help` says what an option does in about a line, and the full
explanation -- the why, the interactions, the defaults' reasoning -- lives in exactly one
place, `docs/cli_reference.md`. Help strings had drifted into carrying whole paragraphs: one
reached 603 characters, and 26 were over 180, restating text the reference already had in more
depth.

These two tests are what makes the convention hold rather than decay. The cap stops a help
string growing back into an essay, and the coverage test is what makes the cap safe: trimming
is only lossless while every option is genuinely documented somewhere else.
"""
import inspect
import pathlib

import pytest

from trnagraph import cli

#: Roughly two terminal lines at the width typer wraps to. The longest surviving string is 174.
MAX_HELP = 180

REFERENCE = pathlib.Path('docs/cli_reference.md')


def _commands():
    found = {'graph': cli.graph}
    for app in (cli.preprocess_app, cli.analyze_app, cli.tools_app):
        for info in app.registered_commands:
            found[info.name] = info.callback
    return found


def _options():
    """(command, parameter name, long flags, help text) for every documented CLI option."""
    for command, fn in sorted(_commands().items()):
        for name, param in inspect.signature(fn).parameters.items():
            decls = getattr(param.default, 'param_decls', None) or []
            flags = [flag for decl in decls for flag in str(decl).split('/')
                     if flag.startswith('--')]
            help_text = getattr(param.default, 'help', None)
            if flags:
                yield command, name, flags, help_text


@pytest.mark.parametrize('option', list(_options()), ids=lambda o: f'{o[0]}.{o[1]}')
def test_help_strings_stay_about_a_line(option):
    command, _, flags, help_text = option
    if not help_text:
        return
    assert len(help_text) <= MAX_HELP, (
        f'`{command} {flags[0]}` help is {len(help_text)} characters. Keep it to roughly a '
        f'line and put the full explanation in docs/cli_reference.md, which is the one place '
        f'it belongs.'
    )


@pytest.mark.parametrize('option', list(_options()), ids=lambda o: f'{o[0]}.{o[1]}')
def test_every_option_appears_in_the_reference(option):
    """What makes the cap safe: a short help string is only acceptable because the reference
    carries the detail. An option documented nowhere would fail here rather than quietly
    becoming the sole description of itself."""
    command, _, flags, _help = option
    body = REFERENCE.read_text()

    assert any(f'`{flag}`' in body for flag in flags), (
        f'`{command} {flags[0]}` is not documented in {REFERENCE}. Help text is deliberately '
        f'brief, so the reference has to carry the explanation.'
    )
