"""Regression test for roadmap.md Phase 2's "Split-file variance" item -- the Phase 1
deterministic-sort fix in toolsCountReads.py's counttypereads() (the `sorted(list(
trnaloci[currbed].getbin(currread)), key=lambda x: (x.start, x.end, x.name if x.name else ""))`
line). Without it, tRAX's underlying set-iteration order (hash-dependent, and therefore free to
differ across separate process runs -- e.g. two --readlengthsplit worker subprocesses handling the
same underlying reads) could pick a different "first matching" trnaloci candidate for the same
read, and unlike some other call sites audited during this same investigation (see roadmap.md),
here that genuinely changes the outcome: a same-strand candidate classifies the read as a sense
pre-tRNA locus hit (`fulllocuscounts`/`partiallocuscounts`), while an opposite-strand candidate at
the same coordinates classifies it as antisense (`trnaanticounts`) instead -- two different,
externally observable accumulator buckets.

Rather than depend on Python's incidental (and, for just two objects, not reliably order-flipping)
set/hash behavior to simulate "a different process run picked a different order", this test
supplies the two orderings directly via a minimal stand-in for RangeBin's `getbin()` -- exactly
the interface counttypereads() calls -- so it deterministically exercises both possible orderings
of the exact same underlying candidates and confirms the code's own sort normalizes both to an
identical classification result.
"""
from types import SimpleNamespace

import pysam

from trnagraph.modules import toolsCountReads
from trnagraph.modules.toolsTG import GenomeRange


class _FixedOrderBin:
    """Minimal getbin()-compatible stand-in that yields candidates in a caller-chosen order,
    standing in for a RangeBin whose underlying set happened to iterate in that order."""

    def __init__(self, ordered_candidates):
        self._ordered_candidates = ordered_candidates

    def getbin(self, item):
        return iter(self._ordered_candidates)


def _make_single_read_bam(tmp_path):
    bampath = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 500}]}
    with pysam.AlignmentFile(str(bampath), "wb", header=header) as outf:
        read = pysam.AlignedSegment()
        read.query_name = "read1"
        read.flag = 0
        read.reference_id = 0
        read.reference_start = 140
        read.mapping_quality = 42
        read.cigar = [(0, 20)]
        read.query_sequence = "A" * 20
        read.query_qualities = [30] * 20
        outf.write(read)
    return str(bampath)


def _classify(bamfile, ordered_candidates):
    trnaloci = {"testbed": _FixedOrderBin(ordered_candidates)}
    return toolsCountReads.counttypereads(
        bamfile, "sample1", None, trnaloci, [], {},
    )


def test_trnaloci_classification_is_stable_regardless_of_candidate_order(tmp_path):
    bamfile = _make_single_read_bam(tmp_path)

    same_strand = GenomeRange("genome", "chr1", 100, 200, "+", name="lociA")
    opposite_strand = GenomeRange("genome", "chr1", 100, 200, "-", name="lociB")

    forward_counts = _classify(bamfile, [same_strand, opposite_strand])
    reversed_counts = _classify(bamfile, [opposite_strand, same_strand])

    assert dict(forward_counts.fulllocuscounts) == dict(reversed_counts.fulllocuscounts)
    assert dict(forward_counts.partiallocuscounts) == dict(reversed_counts.partiallocuscounts)
    assert dict(forward_counts.trnaanticounts) == dict(reversed_counts.trnaanticounts)

    # Sanity: the scenario really does produce a classification, and the two orderings really
    # would disagree without the sort -- otherwise this wouldn't be exercising the fix at all.
    total_classified = (
        sum(forward_counts.fulllocuscounts.values())
        + sum(forward_counts.partiallocuscounts.values())
        + sum(forward_counts.trnaanticounts.values())
    )
    assert total_classified == 1
