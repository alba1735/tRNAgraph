"""Regression tests for the trim summary's output paths (roadmap.md: "`trim_metadata.tsv` wrong
filename for paired-end samples").

Two defects shared one root cause -- the trimmed-output prefix was derived in two places that
disagreed. `_construct_command()` rewrote a bare manifest prefix to `processed/trimmed/<name>`
locally, but that rewrite never reached `self.samples`, so `_generate_summary()` looked for
fastp's JSON report in the working directory, found nothing, and silently skipped every one of
its three outputs -- exactly the path README.md's own Quick Start walks. Where the prefix DID
contain a directory (the only case that got far enough to write anything), the paired-end row
named `{prefix}_merged_trimmed.fastq.gz`, a file fastp never produces.

These test the emitted `trim_metadata.tsv` rather than the naming helper itself: the file is the
template a user hands to `analyze build -i`, so "every path in it names a file that exists" is
the behavior worth pinning.
"""
import json
import logging
import os
import pathlib
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from trnagraph.modules.toolsTrim import FastpTrimmer


def _write_fastq(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("@read\nACGT\n+\nIIII\n")
    return str(path)


def _make_trimmer(tmp_path, entries):
    """Builds a FastpTrimmer over a manifest of (prefix, paired) entries."""
    lines = []
    for prefix, paired in entries:
        r1 = _write_fastq(tmp_path / "raw" / f"{prefix.replace('/', '_')}_R1.fastq")
        if paired:
            r2 = _write_fastq(tmp_path / "raw" / f"{prefix.replace('/', '_')}_R2.fastq")
            lines.append(f"{prefix}\t{r1}\t{r2}")
        else:
            lines.append(f"{prefix}\t{r1}")
    manifest = tmp_path / "manifest.txt"
    manifest.write_text("\n".join(lines) + "\n")
    args = SimpleNamespace(
        input=str(manifest), adapter1=None, adapter2=None, length=15,
        umilength=0, umi3=False, threads=1, colormap=None,
        log=None, quiet=False, verbose=False,
    )
    return FastpTrimmer(args)


def _write_fastp_report(prefix, paired):
    """Writes the minimal fastp JSON report `_generate_summary()` parses."""
    report = {
        'summary': {
            'before_filtering': {'total_reads': 100},
            'after_filtering': {'total_reads': 90},
        },
        'filtering_result': {'too_short_reads': 5, 'low_quality_reads': 5},
        'adapter_cutting': {'adapter_trimmed_reads': 20},
    }
    if paired:
        report['merged_and_filtered'] = {'total_reads': 80}
    path = pathlib.Path(f"{prefix}.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report))


def _read_metadata(path):
    """Returns {sample: fastq} from a written trim_metadata.tsv."""
    rows = pathlib.Path(path).read_text().strip().split("\n")
    assert rows[0].split("\t") == ["fastq", "sample", "group"]
    return {r.split("\t")[1]: r.split("\t")[0] for r in rows[1:]}


@pytest.fixture(autouse=True)
def _no_plot():
    """The summary's trim_feature_types.pdf is out of scope here (and pulls in matplotlib)."""
    with patch("trnagraph.modules.plotsTrimmingStats.visualizer"):
        yield


def test_bare_prefix_paired_end_names_the_merged_file_fastp_actually_writes(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    trimmer = _make_trimmer(tmp_path, [("SampleA", True)])
    _write_fastp_report("processed/trimmed/SampleA", paired=True)

    trimmer._generate_summary()

    metadata = _read_metadata("processed/trimmed/trim_metadata.tsv")
    assert metadata == {"SampleA": "processed/trimmed/SampleA_merged.fastq.gz"}


def test_directory_prefix_paired_end_drops_the_phantom_trimmed_suffix(tmp_path, monkeypatch):
    """The defect as originally reported: fastp writes `{prefix}_merged.fastq.gz` for merged
    output, but the template named `{prefix}_merged_trimmed.fastq.gz`. Only reachable with a
    directory-qualified prefix -- a bare one never got past the JSON lookup above."""
    monkeypatch.chdir(tmp_path)
    trimmer = _make_trimmer(tmp_path, [("out/dir/SampleA", True)])
    _write_fastp_report("out/dir/SampleA", paired=True)

    trimmer._generate_summary()

    metadata = _read_metadata("out/dir/trim_metadata.tsv")
    assert metadata == {"SampleA": "out/dir/SampleA_merged.fastq.gz"}


@pytest.mark.parametrize("prefix,expected_dir", [
    ("SampleA", "processed/trimmed"),
    ("out/dir/SampleA", "out/dir"),
])
def test_single_end_names_the_trimmed_file(tmp_path, monkeypatch, prefix, expected_dir):
    monkeypatch.chdir(tmp_path)
    trimmer = _make_trimmer(tmp_path, [(prefix, False)])
    _write_fastp_report(f"{expected_dir}/SampleA", paired=False)

    trimmer._generate_summary()

    metadata = _read_metadata(f"{expected_dir}/trim_metadata.tsv")
    assert metadata == {"SampleA": f"{expected_dir}/SampleA_trimmed.fastq.gz"}


def test_every_path_in_the_template_names_a_file_the_command_actually_writes(tmp_path, monkeypatch):
    """The property that ties the two halves together: whatever `_construct_command` tells fastp
    to write is what the template points at, for every prefix/layout combination."""
    monkeypatch.chdir(tmp_path)
    entries = [("BareSE", False), ("BarePE", True), ("d/DirSE", False), ("d/DirPE", True)]
    trimmer = _make_trimmer(tmp_path, entries)

    for output_prefix, files in trimmer.samples.items():
        cmd, primary_output = trimmer._construct_command(output_prefix, files)
        _write_fastp_report(output_prefix, paired=bool(files['r2']))
        # fastp is told to write this exact path, either as --merged_out or as --out1
        assert primary_output in cmd

    trimmer._generate_summary()

    written = _read_metadata("processed/trimmed/trim_metadata.tsv")
    expected = {
        os.path.basename(p): trimmer._primary_output(p, f)
        for p, f in trimmer.samples.items()
    }
    assert written == expected


def test_sample_without_a_fastp_report_gets_no_row(tmp_path, monkeypatch):
    """A sample whose fastp run failed has no trimmed output on disk. Listing it anyway would
    reintroduce the very defect this file exists to prevent -- a template naming a file that
    isn't there -- and `analyze build` fails on the missing fastq. The log already names which
    samples failed, and re-running trimming for them is the recovery."""
    monkeypatch.chdir(tmp_path)
    trimmer = _make_trimmer(tmp_path, [("SampleA", True), ("SampleB", True)])
    _write_fastp_report("processed/trimmed/SampleA", paired=True)  # SampleB failed

    trimmer._generate_summary()

    metadata = _read_metadata("processed/trimmed/trim_metadata.tsv")
    assert metadata == {"SampleA": "processed/trimmed/SampleA_merged.fastq.gz"}


def test_manifest_spanning_two_directories_warns_where_the_template_landed(tmp_path, monkeypatch, caplog):
    """One run produces one template, placed alongside the first sample. That is defensible, but
    silently picking one of several directories is expensive to debug, so say which one."""
    monkeypatch.chdir(tmp_path)
    trimmer = _make_trimmer(tmp_path, [("first/SampleA", False), ("second/SampleB", False)])
    _write_fastp_report("first/SampleA", paired=False)
    _write_fastp_report("second/SampleB", paired=False)

    with caplog.at_level(logging.INFO, logger="trnagraph.modules.toolsTrim"):
        trimmer._generate_summary()

    warnings = [r.message for r in caplog.records if r.levelno == logging.WARNING]
    assert any("first/trim_metadata.tsv" in m and "second" in m for m in warnings)
    # the paths inside it stay correct regardless of where the file sits
    assert _read_metadata("first/trim_metadata.tsv") == {
        "SampleA": "first/SampleA_trimmed.fastq.gz",
        "SampleB": "second/SampleB_trimmed.fastq.gz",
    }


def test_single_directory_manifest_does_not_warn(tmp_path, monkeypatch, caplog):
    monkeypatch.chdir(tmp_path)
    trimmer = _make_trimmer(tmp_path, [("SampleA", False), ("SampleB", False)])
    _write_fastp_report("processed/trimmed/SampleA", paired=False)
    _write_fastp_report("processed/trimmed/SampleB", paired=False)

    with caplog.at_level(logging.INFO, logger="trnagraph.modules.toolsTrim"):
        trimmer._generate_summary()

    assert [r.message for r in caplog.records if r.levelno == logging.WARNING] == []
