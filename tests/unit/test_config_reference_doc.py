"""The generated `--config` flags reference in docs/advanced_usage.md.

The template can only ever show `"shrink": null`, which says nothing about whether the key
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

import re

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


def _row_for(body, key):
    """The table row whose Key cell is `key`.

    Matched on the parsed cell rather than a `line.startswith('| `key` |')` prefix: the tables
    are column-aligned, so a row's raw prefix carries however much padding that column needs and
    a prefix match silently finds nothing.
    """
    for line in body.splitlines():
        cells = [c.strip() for c in re.split(r'(?<!\\)\|', line)[1:-1]]
        if cells and cells[0] == f'`{key}`':
            return line
    raise AssertionError(f'no table row for {key!r}')


def test_a_boolean_key_states_both_values(generator):
    """The reported motivation for the table: a bare null gave no clue what it accepted."""
    body = generator.DOC.read_text()
    row = _row_for(body, 'savesplitbams')

    assert '`true`' in row and '`false`' in row, row


def test_the_shrink_key_lists_its_methods(generator):
    """`shrink` is the key that prompted the table: a bare null said nothing about whether it
    wanted true/false, a number, or the name of an estimator."""
    body = generator.DOC.read_text()
    row = _row_for(body, 'shrink')

    assert '"apeGLM"' in row and '"none"' in row, row


def test_a_choice_key_lists_its_choices(generator):
    body = generator.DOC.read_text()
    row = _row_for(body, 'heatorient')

    assert '"vertical"' in row and '"horizontal"' in row, row


def test_the_pipe_in_a_choice_is_escaped_so_the_table_still_renders(generator):
    """An unescaped `|` inside a cell ends the column, which silently mangles the row."""
    assert generator.render_type(__import__('typing').Optional[bool]) == '`true` \\| `false`'


def _cells(row):
    """Split a markdown row on UNESCAPED pipes -- `\\|` is literal text inside a cell, and
    splitting on it would mis-count the columns of every row containing a choice type."""
    import re
    return re.split(r'(?<!\\)\|', row)[1:-1]


def test_generated_tables_are_column_aligned(generator):
    """The generator emits the same padding a markdown formatter would.

    The rows used to be written unpadded (`| a | b |`), which reads fine on GitHub but is not
    what a formatter leaves behind. An editor formatting docs/advanced_usage.md on save
    re-aligned every table in the file, and the committed result then no longer matched what
    this generator produced -- so `test_the_committed_reference_is_up_to_date` failed for a
    reason that had nothing to do with the schema it exists to guard.

    Padding here rather than relaxing the comparison keeps regeneration idempotent under a
    formatter: the generator's output IS the formatted output, so neither can dirty the file
    the other wrote.
    """
    for block in generator.build_table().split('#### ')[1:]:
        rows = [line for line in block.splitlines() if line.startswith('|')]
        widths = {tuple(len(cell) for cell in _cells(row)) for row in rows}

        assert len(widths) == 1, f'ragged table: {sorted(widths)}'


def test_an_underscore_that_could_open_emphasis_is_escaped(generator):
    """Matching the markdown formatter's own escaping, for the same idempotence reason as the
    padding above.

    CommonMark lets `_` open emphasis only when it is NOT intraword, so a formatter escapes
    `<output>_bam` (the underscore follows `>`) and leaves `dedup_method`, `regen_uns` and
    `R1_Path` alone. Escaping every underscore would corrupt those; escaping none leaves the
    committed file permanently one save away from failing the up-to-date test.
    """
    assert generator.escape_markdown('processed/<output>_bam') == r'processed/<output>\_bam'
    assert generator.escape_markdown('dedup_method') == 'dedup_method'
    assert generator.escape_markdown('R1_Path') == 'R1_Path'
    assert generator.escape_markdown('_leading') == r'\_leading'
