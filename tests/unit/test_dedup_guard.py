"""Tests that `preprocess map --dedup` refuses to run on reads with no UMI.

`umi_tools dedup` does not fail on UMI-less input -- it falls back to collapsing reads that
share an alignment position. For most assays that is merely wrong; for tRNA-seq it is
destructive, because a mature tRNA is short and deeply covered, so position-identical reads are
overwhelmingly genuine molecules rather than PCR duplicates. Collapsing them would remove real
signal and leave nothing in the output to indicate it had happened.

The guard therefore refuses rather than warning, and refuses before umi_tools is invoked at all.
"""
from unittest.mock import patch

import pysam
import pytest

from trnagraph.modules.toolsDedup import MissingUMIError, dedup_sample


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


def test_refuses_and_never_invokes_umi_tools_when_no_umi_present(tmp_path):
    bam = _write_bam(
        tmp_path / "SampleA.bam",
        [
            "VH00771:218:AACVJLLHV:1:1102:52921:23794",
            "VH00771:218:AACVJLLHV:1:1210:46559:20045",
        ],
    )

    with patch("trnagraph.modules.toolsDedup.subprocess.run") as mock_run:
        with pytest.raises(MissingUMIError) as excinfo:
            dedup_sample(bam)

    mock_run.assert_not_called()
    # The message has to name the offending file, since a run covers many samples at once.
    assert "SampleA.bam" in str(excinfo.value)
