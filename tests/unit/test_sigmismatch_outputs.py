"""Regression tests for the significant-mismatch outputs written by `toolsGetCoverage`.

Two distinct tRAX artifacts had collided on one filename here. tRNAgraph's
`mismatch/<exp>-sigmismatch.txt` was a BED lifted from `tRAX/getgenomicmismatches.py`
(a standalone script `processsamples.py` never invokes), while tRAX's file of that name
is written by `newcoverageplots.R:581` and is a row filter over the size-factor-normalized
coverage table. The BED additionally indexed *alignment columns* where its tRAX source
indexes ungapped sequence positions, so `getbase()` was handed a column number and
produced coordinates past the end of the feature -- measured on vibrChol1, BED ends ran
to 124 on tRNAs whose padded contig is only 116 long.

These tests pin both artifacts: the BED emits in-range sequence coordinates, and the
table reproduces tRAX's shape.
"""
import io

import pytest

from trnagraph.modules import toolsGetCoverage as tgc
from trnagraph.modules import toolsTG


def _coverage(region, values):
    """Build a ReadCoverage over `region` pre-loaded with per-sequence-position counts."""
    cov = tgc.ReadCoverage(region)
    assert len(values) == region.length(), "test fixture must cover every sequence position"
    cov.coverage = list(values)
    cov.totalreads = sum(values)
    return cov


def _coverage_info(region, *, coverage, mismatches, deletions=None, bases=None):
    name = region.name
    zeros = [0] * region.length()
    empty = lambda: {name: _coverage(region, zeros)}
    return tgc.CoverageInfo(
        readcounts={name: sum(coverage)},
        allcoverages={name: _coverage(region, coverage)},
        readstarts=empty(), readends=empty(),
        multaminocoverages=empty(), multaccoverages=empty(), multtrnacoverages=empty(),
        uniquecoverages={name: _coverage(region, coverage)},
        uniquegenomecoverages=empty(), multigenomecoverages=empty(),
        readmismatches={name: _coverage(region, mismatches)},
        adeninemismatches={name: _coverage(region, bases or zeros)},
        thyminemismatches=empty(), cytosinemismatches=empty(), guanosinemismatches=empty(),
        readskips={name: _coverage(region, deletions or zeros)},
        trimcoverage=empty(), trimmismatches=empty(),
    )


class _Samples:
    def __init__(self, samples):
        self._samples = list(samples)

    def getsamples(self):
        return list(self._samples)


# --- significant_mismatch_positions ------------------------------------------------

def test_positions_are_sequence_indices_never_beyond_the_feature():
    """The defect: indices came from the alignment, which is wider than the sequence."""
    length = 5
    cov = {"s1": [100] * length}
    mis = {"s1": [0, 90, 0, 0, 80]}

    hits = list(tgc.significant_mismatch_positions(cov, mis))

    assert [pos for pos, _ in hits] == [1, 4]
    assert all(pos < length for pos, _ in hits)


def test_rate_uses_raw_counts_with_the_trax_pseudocount():
    """tRAX's getgenomicmismatches.py:181 is `mismatch / (10 + coverage)` on raw counts."""
    cov = {"s1": [90]}
    mis = {"s1": [6]}

    (pos, percent), = tgc.significant_mismatch_positions(cov, mis)

    assert pos == 0
    assert percent == pytest.approx(6 / 100)


def test_rate_is_the_maximum_across_samples():
    cov = {"a": [1000], "b": [10]}
    mis = {"a": [1], "b": [9]}

    (_, percent), = tgc.significant_mismatch_positions(cov, mis)

    assert percent == pytest.approx(9 / 20)


def test_threshold_is_exclusive_and_defaults_to_the_trax_value():
    assert tgc.SIGMISMATCH_BED_THRESHOLD == 0.05
    assert tgc.SIGMISMATCH_BED_PSEUDOCOUNT == 10

    # exactly at the threshold: 5 / (10 + 90) == 0.05, must not be reported
    assert list(tgc.significant_mismatch_positions({"s": [90]}, {"s": [5]})) == []
    assert [p for p, _ in tgc.significant_mismatch_positions({"s": [90]}, {"s": [6]})] == [0]


def test_pseudocount_is_configurable_for_graph_time_callers():
    cov = {"s": [90]}
    mis = {"s": [5]}

    # 5 / (10 + 90) == 0.05 exactly, which the exclusive threshold rejects...
    assert list(tgc.significant_mismatch_positions(cov, mis)) == []
    # ...while 5 / 90 clears it once the damping is removed.
    hits = list(tgc.significant_mismatch_positions(cov, mis, pseudocount=0))
    assert [p for p, _ in hits] == [0]


# --- BED emission through transcriptcoverage ---------------------------------------

def _run_transcriptcoverage(alignment, coverage, mismatches, positionnums=None):
    """Drive transcriptcoverage over one feature and return (coverage_rows, bed, table)."""
    seqlen = len([c for c in alignment if c not in tgc.gapchars])
    region = toolsTG.GenomeRange("db", "tRNA-Test-AAA-1", 20, 20 + seqlen, "+", name="tRNA-Test-AAA-1")
    info = _coverage_info(region, coverage=coverage, mismatches=mismatches)

    stk = toolsTG.RnaAlignment({region.name: alignment}, "." * len(alignment))
    if positionnums is None:
        positionnums = [str(i + 1) for i in range(len(alignment))]

    covout, bedout, tableout = io.StringIO(), io.StringIO(), io.StringIO()
    tgc.transcriptcoverage(
        {"s1": info}, covout, [region], _Samples(["s1"]), {"s1": 1.0},
        stk, positionnums, sigmismatch=bedout, sigmismatchtable=tableout,
    )
    return covout.getvalue(), bedout.getvalue(), tableout.getvalue()


def test_bed_coordinates_stay_inside_the_feature_when_the_alignment_has_gaps():
    """The real-data failure: alignment columns were used as genomic offsets."""
    # 6 sequence bases spread across a 10-column alignment; the last base is column 9.
    alignment = "AC--GT--AC"
    coverage = [100] * 6
    mismatches = [0, 0, 0, 0, 0, 90]

    _, bed, _ = _run_transcriptcoverage(alignment, coverage, mismatches)

    fields = bed.strip().split("\t")
    assert fields[0] == "tRNA-Test-AAA-1"
    # sequence index 5 -> start 20 + 5; the buggy version emitted column 9 -> start 29
    assert (int(fields[1]), int(fields[2])) == (25, 26)
    assert fields[3] == "tRNA-Test-AAA-1_5pos"
    assert int(fields[1]) >= 20 and int(fields[2]) <= 26


def test_bed_never_emits_past_the_feature_end():
    alignment = "A-C-G-T-A-C-G-T"
    seqlen = 8
    coverage = [100] * seqlen
    mismatches = [90] * seqlen

    region_end = 20 + seqlen
    _, bed, _ = _run_transcriptcoverage(alignment, coverage, mismatches)

    ends = [int(line.split("\t")[2]) for line in bed.strip().splitlines()]
    assert len(ends) == seqlen
    assert max(ends) <= region_end


def test_bed_emission_does_not_swallow_errors():
    """The old block wrapped the emit in a bare `except Exception: pass`."""
    import inspect
    source = inspect.getsource(tgc.transcriptcoverage)
    assert "except Exception" not in source


# --- the tRAX-shaped sigmismatch table ---------------------------------------------

def test_table_filters_on_normalized_coverage_at_the_trax_threshold():
    """newcoverageplots.R:581 -- mismatchedbases / (coverage + 30) > .1."""
    assert tgc.SIGMISMATCH_TABLE_THRESHOLD == 0.1
    assert tgc.SIGMISMATCH_TABLE_PSEUDOCOUNT == 30

    alignment = "ACGT"
    coverage = [70, 70, 70, 70]
    # 7/(70+30) = .07 (below); 11/(70+30) = .11 (above)
    mismatches = [7, 11, 7, 11]

    _, _, table = _run_transcriptcoverage(alignment, coverage, mismatches)

    rows = table.strip().splitlines()
    assert len(rows) == 2
    assert all('"tRNA-Test-AAA-1" "s1"' in row for row in rows)


def test_table_rows_carry_r_write_table_quoting_and_a_row_index():
    alignment = "ACGT"
    _, _, table = _run_transcriptcoverage(alignment, [70] * 4, [0, 50, 0, 0])

    row, = table.strip().splitlines()
    fields = row.split(" ")
    # rowname, Feature, Sample, position, then 19 columns total after the rowname
    assert len(fields) == 20
    assert fields[0] == '"2"', "row index is the 1-based position in the coverage table"
    assert fields[1] == '"tRNA-Test-AAA-1"'
    assert fields[2] == '"s1"'
    assert fields[3] == '"2"'
    assert fields[12] == '"C"', "actualbase is quoted, every other trailing column is numeric"


def test_table_row_index_continues_across_chunks():
    """main() calls transcriptcoverage once per 50-feature chunk; indices must not restart."""
    alignment = "ACGT"
    seqlen = 4
    region = toolsTG.GenomeRange("db", "tRNA-Test-AAA-1", 20, 20 + seqlen, "+", name="tRNA-Test-AAA-1")
    info = _coverage_info(region, coverage=[70] * seqlen, mismatches=[0, 50, 0, 0])
    stk = toolsTG.RnaAlignment({region.name: alignment}, "." * seqlen)
    nums = [str(i + 1) for i in range(seqlen)]

    covout, bedout, tableout = io.StringIO(), io.StringIO(), io.StringIO()
    written = tgc.transcriptcoverage({"s1": info}, covout, [region], _Samples(["s1"]), {"s1": 1.0},
                                     stk, nums, sigmismatch=bedout, sigmismatchtable=tableout)
    assert written == seqlen

    tgc.transcriptcoverage({"s1": info}, covout, [region], _Samples(["s1"]), {"s1": 1.0},
                           stk, nums, sigmismatch=bedout, sigmismatchtable=tableout,
                           rowoffset=written)

    first, second = tableout.getvalue().strip().splitlines()
    assert first.split(" ")[0] == '"2"'
    assert second.split(" ")[0] == '"6"'


def test_no_table_or_bed_written_when_not_requested():
    alignment = "ACGT"
    region = toolsTG.GenomeRange("db", "tRNA-Test-AAA-1", 20, 24, "+", name="tRNA-Test-AAA-1")
    info = _coverage_info(region, coverage=[70] * 4, mismatches=[50] * 4)
    stk = toolsTG.RnaAlignment({region.name: alignment}, "." * 4)

    covout = io.StringIO()
    tgc.transcriptcoverage({"s1": info}, covout, [region], _Samples(["s1"]), {"s1": 1.0},
                           stk, ["1", "2", "3", "4"])

    assert len(covout.getvalue().strip().splitlines()) == 4
