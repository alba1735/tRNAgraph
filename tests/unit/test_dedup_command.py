"""Tests the command line tRNAgraph hands to `umi_tools dedup`.

This asserts against subprocess argv rather than observable output because the argv *is* the
contract with an external binary -- there is no in-process boundary between the two, and the
alternative (running umi_tools for real) is an integration test, not a unit one. What the
deduplication itself does to the reads is umi_tools' behavior and is deliberately not asserted
here.

The separator argument is the load-bearing one: umi_tools defaults to `--umi-separator=_`, so a
BAM carrying fastp's colon-delimited names would be parsed with the UMI silently mis-read rather
than failing.
"""
import pathlib
from unittest.mock import patch

import pysam
import pytest

from trnagraph.modules.toolsDedup import dedup_sample


def _write_bam(path, read_names):
    pathlib.Path(path).parent.mkdir(parents=True, exist_ok=True)
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


UNDERSCORE_READS = [
    "VH00771:218:AACVJLLHV:1:1102:52921:23794_AAAACAC",
    "VH00771:218:AACVJLLHV:1:1210:46559:20045_CAAGCGC",
]
COLON_READS = [
    "VH00771:218:AACVJLLHV:1:1102:52921:23794:AAAACAC",
    "VH00771:218:AACVJLLHV:1:1210:46559:20045:CAAGCGC",
]


def _captured_command(tmp_path, reads, **kwargs):
    bam = _write_bam(tmp_path / "SampleA.bam", reads)
    captured = {}

    def fake_run(cmd, *args, **runkwargs):
        captured["cmd"] = cmd
        stdout_path = next(p.split("=", 1)[1] for p in cmd if p.startswith("--stdout="))
        _write_bam(stdout_path, reads[:1])
        return type("Result", (), {"returncode": 0, "stdout": b"", "stderr": b""})()

    with patch("trnagraph.modules.toolsDedup.subprocess.run", side_effect=fake_run):
        dedup_sample(bam, **kwargs)
    return captured["cmd"]


def test_separator_argument_follows_the_detected_read_name_shape(tmp_path):
    underscore_cmd = _captured_command(tmp_path / "u", UNDERSCORE_READS)
    colon_cmd = _captured_command(tmp_path / "c", COLON_READS)

    assert "--umi-separator=_" in underscore_cmd
    assert "--umi-separator=:" in colon_cmd


def test_defaults_to_directional_and_honors_an_explicit_method(tmp_path):
    default_cmd = _captured_command(tmp_path / "d", UNDERSCORE_READS)
    unique_cmd = _captured_command(tmp_path / "u", UNDERSCORE_READS, method="unique")

    assert default_cmd[default_cmd.index("--method") + 1] == "directional"
    assert unique_cmd[unique_cmd.index("--method") + 1] == "unique"


def test_reads_the_moved_aside_original_and_writes_the_sample_name(tmp_path):
    """umi_tools cannot write in place, so the two paths must not be the same file."""
    cmd = _captured_command(tmp_path, UNDERSCORE_READS)

    stdin = next(p.split("=", 1)[1] for p in cmd if p.startswith("--stdin="))
    stdout = next(p.split("=", 1)[1] for p in cmd if p.startswith("--stdout="))

    assert stdin.endswith(".prededup.bam")
    assert stdout.endswith("SampleA.bam")
    assert stdin != stdout
