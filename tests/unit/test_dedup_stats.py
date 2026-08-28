"""Per-sample deduplication statistics, parsed from umi_tools' own logs.

Whether deduplication behaved sensibly is not visible from the deduplicated bams alone. The
numbers that answer it -- how many reads went in, how many molecules came out, how concentrated
those reads were, and how close the UMI space came to its ceiling -- are all reported by
umi_tools and then thrown away. Parsing them costs nothing (the log is already written next to
each bam) whereas recomputing them would mean a second full pass over every file.

Read together they say whether a sample should be trusted: a low reads-per-molecule with few
distinct positions is a low-complexity library, while a high maximum-UMIs-per-position against
the 4^n ceiling for an n-base UMI means the tag space is being exhausted and counts at the
deepest features are compressed.
"""
import pytest

from trnagraph.modules.toolsDedup import parse_dedup_log, dedup_stats_row

LOG = """\
2026-08-27 14:09:24,765 INFO Parsed 1000000 input reads
2026-08-27 14:26:26,363 INFO Reads: Input Reads: 69351623
2026-08-27 14:26:26,363 INFO Number of reads out: 16078732
2026-08-27 14:26:26,363 INFO Total number of positions deduplicated: 11200496
2026-08-27 14:26:26,363 INFO Mean number of unique UMIs per position: 1.77
2026-08-27 14:26:26,363 INFO Max. number of unique UMIs per position: 4267
"""


def test_parses_every_reported_field(tmp_path):
    p = tmp_path / "D0_1_dedup.log"
    p.write_text(LOG)

    stats = parse_dedup_log(str(p))

    assert stats["input_reads"] == 69351623
    assert stats["output_reads"] == 16078732
    assert stats["positions"] == 11200496
    assert stats["mean_umis_per_position"] == pytest.approx(1.77)
    assert stats["max_umis_per_position"] == 4267


def test_missing_log_yields_no_stats(tmp_path):
    """An older umi_tools, or a log that never got written, must not break the run."""
    assert parse_dedup_log(str(tmp_path / "absent.log")) == {}


def test_derived_rates_are_computed_from_the_parsed_counts(tmp_path):
    p = tmp_path / "D0_1_dedup.log"
    p.write_text(LOG)

    row = dedup_stats_row("D0_1", str(p))

    assert row["sample"] == "D0_1"
    # 16078732 / 69351623
    assert row["retained_pct"] == pytest.approx(23.18, abs=0.01)
    # 69351623 / 16078732
    assert row["reads_per_molecule"] == pytest.approx(4.31, abs=0.01)
    # 69351623 / 11200496
    assert row["reads_per_position"] == pytest.approx(6.19, abs=0.01)


def test_row_survives_an_unparseable_log(tmp_path):
    p = tmp_path / "broken_dedup.log"
    p.write_text("2026-08-27 INFO something else entirely\n")

    row = dedup_stats_row("broken", str(p))

    assert row["sample"] == "broken"
    assert row["input_reads"] is None
    assert row["reads_per_molecule"] is None
