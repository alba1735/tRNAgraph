#!/usr/bin/env python3
"""
Regenerate the `--config` flags reference in docs/advanced_usage.md.

The table is written between marker comments in the doc and committed, so the documentation
reads correctly on GitHub with no build step. A unit test regenerates it and fails when the
committed text no longer matches the models, which is what stops it going stale the way a
hand-written table would.

Two sources, each already authoritative for its half: the pydantic models in toolsSchemas say
which keys exist and what TYPE each takes, and the typer commands in cli.py say what each one
DEFAULTS to and what it means. Importing cli here is fine -- this is a dev script, not the
package's own import graph, which is why toolsSchemas itself still declares its fields by hand.

Usage:  python scripts/generate_config_reference.py [--check]
"""
from __future__ import annotations

import argparse
import inspect
import pathlib
import re
import sys
import typing

BEGIN = '<!-- BEGIN GENERATED: config-flags (scripts/generate_config_reference.py) -->'
END = '<!-- END GENERATED: config-flags -->'

DOC = pathlib.Path(__file__).resolve().parents[1] / 'docs' / 'advanced_usage.md'

#: The order blocks appear in, which is pipeline order rather than alphabetical: someone
#: reading this is following a run from one end to the other.
BLOCK_ORDER = ['makedb', 'trim', 'map', 'build', 'addsplit', 'cluster', 'graph', 'log2fc']

BLOCK_TITLES = {
    'makedb': '`preprocess makedb`',
    'trim': '`preprocess trim`',
    'map': '`preprocess map`',
    'build': '`analyze build`',
    'addsplit': '`analyze addsplit`',
    'cluster': '`analyze cluster`',
    'graph': '`graph`',
    'log2fc': '`tools log2fc`',
}


def _unwrap_optional(annotation):
    """`Optional[X]` -> `X`; every field on these models is Optional by construction."""
    if typing.get_origin(annotation) is typing.Union:
        args = [a for a in typing.get_args(annotation) if a is not type(None)]
        if len(args) == 1:
            return args[0]
    return annotation


def _plural(rendered: str) -> str:
    return {'string': 'strings', 'integer': 'integers', 'number': 'numbers'}.get(
        rendered, rendered)


def render_type(annotation) -> str:
    """A JSON-shaped description of what a key accepts.

    This is the whole reason the table exists: a template showing `"lfcshrink": null` says
    nothing about whether it wants true/false, a number, or one of a set of words.
    """
    inner = _unwrap_optional(annotation)
    origin = typing.get_origin(inner)

    if origin is typing.Literal:
        return ' \\| '.join(f'`"{value}"`' for value in typing.get_args(inner))
    if origin is list:
        (item,) = typing.get_args(inner) or (str,)
        return f'list of {_plural(render_type(item))}'
    return {
        bool: '`true` \\| `false`',
        int: 'integer',
        float: 'number',
        str: 'string',
    }.get(inner, f'`{getattr(inner, "__name__", inner)}`')


def render_default(param) -> str:
    """What the command does when the key is absent, taken from the typer signature."""
    if param is None:
        return '—'
    default = getattr(param.default, 'default', param.default)
    if default is inspect.Parameter.empty or default is Ellipsis or default is None:
        return '—'
    if isinstance(default, bool):
        return f'`{str(default).lower()}`'
    if isinstance(default, (list, tuple)):
        return '`' + ', '.join(str(v) for v in default) + '`' if default else '—'
    return f'`{default}`'


def render_help(param) -> str:
    """The command's own one-line description, trimmed to its first sentence."""
    if param is None:
        return ''
    text = (getattr(param.default, 'help', None) or '').strip()
    if not text:
        return ''
    # First sentence only: several help strings carry a second clause aimed at the terminal.
    text = re.split(r'(?<=[a-z0-9)\'"])\. ', text)[0].rstrip('.')
    return text.replace('|', '\\|').replace('\n', ' ')


def _command_params():
    """{command: {option name: inspect.Parameter}} for every command with a block."""
    from trnagraph import cli

    commands = {'graph': cli.graph}
    for app in (cli.preprocess_app, cli.analyze_app, cli.tools_app):
        for info in app.registered_commands:
            commands[info.name] = info.callback
    return {name: dict(inspect.signature(fn).parameters) for name, fn in commands.items()}


def build_table() -> str:
    from trnagraph.modules.toolsSchemas import COMMAND_FLAG_MODELS

    params = _command_params()
    out = []
    for command in BLOCK_ORDER:
        model = COMMAND_FLAG_MODELS[command]
        out.append(f'#### `flags.{command}` — {BLOCK_TITLES[command]}\n')
        out.append('| Key | Accepts | Default | Meaning |')
        out.append('| --- | --- | --- | --- |')
        for key in sorted(model.model_fields):
            field = model.model_fields[key]
            param = params[command].get(key)
            out.append(f'| `{key}` | {render_type(field.annotation)} | '
                       f'{render_default(param)} | {render_help(param)} |')
        out.append('')
    return '\n'.join(out).rstrip() + '\n'


def render_document(existing: str, table: str) -> str:
    if BEGIN not in existing or END not in existing:
        raise SystemExit(f'{DOC} is missing the {BEGIN!r} / {END!r} markers.')
    head, rest = existing.split(BEGIN, 1)
    _, tail = rest.split(END, 1)
    return f'{head}{BEGIN}\n\n{table}\n{END}{tail}'


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true',
                        help='exit non-zero if the committed table is out of date')
    options = parser.parse_args(argv)

    existing = DOC.read_text()
    updated = render_document(existing, build_table())
    if options.check:
        if updated != existing:
            print(f'{DOC} is out of date; run: python scripts/generate_config_reference.py')
            return 1
        return 0
    DOC.write_text(updated)
    print(f'Wrote the config flags reference into {DOC}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
