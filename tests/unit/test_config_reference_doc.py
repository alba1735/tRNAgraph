"""The generated `--config` flags reference in docs/advanced_usage.md.

The template can only ever show `"lfcshrink": null`, which says nothing about whether the key
wants true/false, a number, or one of a set of words. The reference tables answer that, and
are generated from the pydantic models (which own the types) and the typer commands (which own
the defaults and the one-line meanings) so neither becomes a second source of truth.

The table is committed rather than built on demand, so the doc reads correctly on GitHub with
no build step. This test is what stops the committed copy going stale: it regenerates and
compares, exactly as `--check` does, and points at the script when they differ. It is the same
guarantee the schema already has against the CLI via the drift tests.
"""
import importlib.util
import pathlib

import pytest

SCRIPT = pathlib.Path('scripts/generate_config_reference.py')


@pytest.fixture(scope='module')
def generator():
    spec = importlib.util.spec_from_file_location('generate_config_reference', SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_committed_reference_is_up_to_date(generator):
    existing = generator.DOC.read_text()
    regenerated = generator.render_document(existing, generator.build_table())

    assert regenerated == existing, (
        'docs/advanced_usage.md no longer matches the schema. Regenerate it with:\n'
        '    python scripts/generate_config_reference.py'
    )


def test_check_mode_agrees_with_the_committed_file(generator):
    assert generator.main(['--check']) == 0


def test_every_block_has_a_table(generator):
    body = generator.DOC.read_text()
    for command in generator.BLOCK_ORDER:
        assert f'#### `flags.{command}`' in body


def test_a_boolean_key_states_both_values(generator):
    """The reported motivation for the table: `lfcshrink` gave no clue what it accepted."""
    body = generator.DOC.read_text()
    row = next(line for line in body.splitlines() if line.startswith('| `lfcshrink` |'))

    assert '`true`' in row and '`false`' in row, row


def test_a_choice_key_lists_its_choices(generator):
    body = generator.DOC.read_text()
    row = next(line for line in body.splitlines() if line.startswith('| `heatorient` |'))

    assert '"vertical"' in row and '"horizontal"' in row, row


def test_the_pipe_in_a_choice_is_escaped_so_the_table_still_renders(generator):
    """An unescaped `|` inside a cell ends the column, which silently mangles the row."""
    assert generator.render_type(__import__('typing').Optional[bool]) == '`true` \\| `false`'
