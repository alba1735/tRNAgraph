"""Tests what `preprocess map --dedup` leaves on disk.

`umi_tools dedup` cannot write in place -- it reads one file and writes another -- so the
wrapper has to decide what the BAM directory looks like afterwards. The decision is that the
deduplicated reads take over the sample's ordinary name, so `analyze build --bamdir` needs no
knowledge that deduplication happened and no downstream path changes. The pre-deduplication BAM
is discarded by default and kept only on request, because mapping is the expensive step and
comparing deduplicated against non-deduplicated output is a real workflow.

umi_tools is stubbed here: whether it removes the right reads is its own behavior, not
tRNAgraph's. What is tested is the file handling around it.
"""
from unittest.mock import patch

import pysam
import pytest

from trnagraph.modules.toolsDedup import dedup_sample


def _write_bam(path, read_names):
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


UMI_READS = [
    "VH00771:218:AACVJLLHV:1:1102:52921:23794_AAAACAC",
    "VH00771:218:AACVJLLHV:1:1210:46559:20045_AAAACAC",
    "VH00771:218:AACVJLLHV:1:1212:8517:41571_CAAGCGC",
]


def _fake_umi_tools(deduped_reads):
    """Stands in for umi_tools: writes a BAM at the path given by --stdout."""
    def run(cmd, *args, **kwargs):
        stdout_path = next(
            part.split("=", 1)[1] for part in cmd if part.startswith("--stdout=")
        )
        _write_bam(stdout_path, deduped_reads)
        return type("Result", (), {"returncode": 0, "stdout": b"", "stderr": b""})()
    return run


def _read_names(path):
    with pysam.AlignmentFile(str(path), "rb") as bam:
        return [record.query_name for record in bam]


def test_deduplicated_reads_take_over_the_sample_bam_name(tmp_path):
    """The default: downstream sees an ordinary <sample>.bam, and nothing else is left behind."""
    bam = _write_bam(tmp_path / "SampleA.bam", UMI_READS)

    with patch("trnagraph.modules.toolsDedup.subprocess.run",
               side_effect=_fake_umi_tools(UMI_READS[:2])):
        dedup_sample(bam)

    assert _read_names(bam) == UMI_READS[:2]
    assert not (tmp_path / "SampleA.prededup.bam").exists()


def test_keep_prededup_retains_the_original_alongside(tmp_path):
    """Opt-in: both files present, so a deduplicated/non-deduplicated comparison needs no remap."""
    bam = _write_bam(tmp_path / "SampleA.bam", UMI_READS)

    with patch("trnagraph.modules.toolsDedup.subprocess.run",
               side_effect=_fake_umi_tools(UMI_READS[:2])):
        dedup_sample(bam, keep_prededup=True)

    assert _read_names(bam) == UMI_READS[:2]
    assert _read_names(tmp_path / "SampleA.prededup.bam") == UMI_READS


def test_failed_dedup_restores_the_original_bam(tmp_path):
    """A umi_tools failure must not cost the mapping.

    The original is moved aside before umi_tools runs, so a failure that went unnoticed would
    leave the sample's BAM missing or truncated while the only good copy had already been
    deleted -- and mapping is by far the most expensive step to repeat.
    """
    bam = _write_bam(tmp_path / "SampleA.bam", UMI_READS)

    def failing_umi_tools(cmd, *args, **kwargs):
        return type("Result", (), {"returncode": 1, "stdout": b"", "stderr": b"boom"})()

    with patch("trnagraph.modules.toolsDedup.subprocess.run", side_effect=failing_umi_tools):
        with pytest.raises(RuntimeError):
            dedup_sample(bam)

    assert _read_names(bam) == UMI_READS
    assert not (tmp_path / "SampleA.prededup.bam").exists()
