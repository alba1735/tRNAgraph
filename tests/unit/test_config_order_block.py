"""`--config` gains a top-level `order` block declaring explicit category order.

Every plotting module falls back to `sorted()`, so a timecourse's groups render alphabetically
-- Day 0, Day 35, Day 70 happens to survive that, but D0/D14/D35/D70 does not, and neither
does any label whose chronology is not its alphabet. Order is also what picks the DESeq2
reference level, so it cannot live in a per-command flags block: it describes what the
experiment IS, not how one command ran, which is the same line the schema already draws for
the obs/var filters.

Top-level, therefore, as a sibling of `obs_r` rather than a key under `flags.build`.
"""

from trnagraph.modules.toolsSchemas import RunConfig


def test_order_is_accepted_at_the_top_level():
    config = RunConfig.model_validate({
        'name': 'organoid_timecourse',
        'order': {'timepoint': ['Day 0', 'Day 35', 'Day 70']},
    })

    assert config.order == {'timepoint': ['Day 0', 'Day 35', 'Day 70']}


def test_order_is_optional():
    """Every existing config file predates this key and must keep validating."""
    assert RunConfig.model_validate({'name': 'plain'}).order is None
