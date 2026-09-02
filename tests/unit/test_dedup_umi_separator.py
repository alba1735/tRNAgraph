"""Tests for UMI-separator detection on mapped BAMs (`preprocess map --dedup`).

tRNAgraph can produce read names in two different shapes, and `umi_tools dedup` has to be
told which one it is looking at:

  * `preprocess trim -u N`          -- fastp appends the UMI to the read name, and tRNAgraph
                                       pins the delimiter to `_` (fastp's own default is `:`).
  * `preprocess trim -u N --umi3`   -- fastp has no tail-anchored UMI option, so tRNAgraph
                                       shells out to `umi_tools extract`, which uses `_`.

Externally-trimmed data is the reason this is detected rather than assumed: a BAM trimmed by
someone else's pipeline (or by an older tRNAgraph that let fastp's `:` default stand) still has
to dedup correctly. `umi_tools dedup` defaults to `--umi-separator=_`, so guessing wrong does
not error -- it mis-parses the UMI, which is exactly the kind of silent wrong answer this
detection exists to prevent.
"""
import pysam

from trnagraph.modules.toolsDedup import detect_umi_separator


def _write_bam(path, read_names):
    """Writes a minimal single-reference BAM carrying the given read names."""
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"SN": "tRNA-Ala-AGC-1", "LN": 100}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for name in read_names:
            record = pysam.AlignedSegment()
            record.query_name = name
            record.query_sequence = "ACGTACGTAC"
            record.flag = 0
            record.reference_id = 0
            record.reference_start = 10
            record.mapping_quality = 255
            record.cigartuples = [(0, 10)]
            record.query_qualities = pysam.qualitystring_to_array("I" * 10)
            bam.write(record)
    return str(path)


def test_detects_underscore_separator(tmp_path):
    """`umi_tools extract` style names -- the shape the real hg38 OTTR BAMs carry."""
    bam = _write_bam(
        tmp_path / "underscore.bam",
        [
            "VH00771:218:AACVJLLHV:1:1102:52921:23794_AAAACAC",
            "VH00771:218:AACVJLLHV:1:1210:46559:20045_AAAACAC",
            "VH00771:218:AACVJLLHV:1:1212:8517:41571_CAAGCGC",
        ],
    )

    assert detect_umi_separator(bam) == "_"


def test_detects_colon_separator_against_colon_delimited_illumina_names(tmp_path):
    """fastp's own default delimiter, which older/external trims leave in place.

    This is the case that cannot be decided by "does the name contain the character": an
    Illumina read name is already colon-delimited, so every name contains colons whether or
    not a UMI was appended. What distinguishes the two is that the appended field is a short
    nucleotide string of consistent length.
    """
    bam = _write_bam(
        tmp_path / "colon.bam",
        [
            "VH00771:218:AACVJLLHV:1:1102:52921:23794:AAAACAC",
            "VH00771:218:AACVJLLHV:1:1210:46559:20045:AAAACAC",
            "VH00771:218:AACVJLLHV:1:1212:8517:41571:CAAGCGC",
        ],
    )

    assert detect_umi_separator(bam) == ":"


def test_plain_illumina_names_report_no_umi(tmp_path):
    """Untrimmed-for-UMI names must come back as None, not as a false colon match.

    This is the case the refusal depends on: `umi_tools dedup` given UMI-less reads collapses
    them by alignment position instead, and in tRNA-seq position-identical reads are
    overwhelmingly real, so a false positive here would silently delete genuine signal.
    """
    bam = _write_bam(
        tmp_path / "no_umi.bam",
        [
            "VH00771:218:AACVJLLHV:1:1102:52921:23794",
            "VH00771:218:AACVJLLHV:1:1210:46559:20045",
            "VH00771:218:AACVJLLHV:1:1212:8517:41571",
        ],
    )

    assert detect_umi_separator(bam) is None
