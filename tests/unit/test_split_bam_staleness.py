"""A cached split BAM is only reused while it is newer than the BAM it was cut from.

The incident this comes from: on the hg38 c305 dataset, `preprocess map --dedup` wrote the
deduplicated BAMs into processed/bam_dedup/ over the pre-dedup ones, leaving the u60/ and o60/
subdirectories -- cut from the PRE-dedup parents -- in place beside them. The next build found
every split BAM present, logged "Split BAM files found for cutoff 60. Skipping split", and
counted pre-dedup reads into both split variants. u60 + o60 came to 68.7M reads against the full
variant's 6.97M: a tenfold overstatement, invisible in every figure downstream because a count
carries no provenance. The physical giveaway was a 1.44 GB u60 split whose parent was 621 MB --
a length-filtered subset larger than the whole.

Presence was the only thing the cache checked. It now also checks order, and REFUSES rather than
silently re-cutting: re-splitting rewrites gigabytes off an mtime comparison, and a silent fix
would have produced correct numbers while leaving the user believing the splits had been sound.
"""
import logging
import os
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from trnagraph.modules import toolsTG
from trnagraph.modules.adataBuild import (AnnDataBuilder, _refuse_stale_split_bams,
                                          _stale_split_bams)


def _bare_builder():
    builder = object.__new__(AnnDataBuilder)
    builder.logger = logging.getLogger("trnagraph.modules.adataBuild")
    builder._split_dirs_preexisted = {}
    return builder


def _make_manifest(tmp_path, samplenames):
    manifest = tmp_path / "manifest.txt"
    lines = ["fastq\tsample\tgroup"]
    lines += [f"{name}.fastq\t{name}\t{name}" for name in samplenames]
    manifest.write_text("\n".join(lines) + "\n")
    return str(manifest)


def _bamdir(tmp_path, samplenames, cutoff, split_older):
    """A populated bamdir whose splits are written either before or after their parents."""
    bamdir = tmp_path / "bam"
    bamdir.mkdir()
    for tag in (f"u{cutoff}", f"o{cutoff}"):
        (bamdir / tag).mkdir()
    early, late = 1_000_000, 2_000_000
    split_time, parent_time = (early, late) if split_older else (late, early)
    for name in samplenames:
        parent = bamdir / f"{name}.bam"
        parent.write_text("")
        os.utime(parent, (parent_time, parent_time))
        for tag in (f"u{cutoff}", f"o{cutoff}"):
            split = bamdir / tag / f"{name}.bam"
            split.write_text("")
            os.utime(split, (split_time, split_time))
    return bamdir


def test_a_split_older_than_its_parent_is_detected(tmp_path):
    bamdir = _bamdir(tmp_path, ["s0", "s1"], 60, split_older=True)

    stale = _stale_split_bams(str(bamdir), 60, ["s0", "s1"])

    assert {tag for tag, _, _ in stale} == {"u60", "o60"}
    assert {sample for _, sample, _ in stale} == {"s0", "s1"}


def test_a_split_newer_than_its_parent_is_accepted(tmp_path):
    """The ordinary case, and the one that must stay silent: splitting writes the child after
    the parent, so a healthy cache always looks like this."""
    bamdir = _bamdir(tmp_path, ["s0", "s1"], 60, split_older=False)

    assert _stale_split_bams(str(bamdir), 60, ["s0", "s1"]) == []


def test_the_refusal_names_the_offending_files_and_the_flag(tmp_path):
    """An error that says only "something is stale" leaves the user with nothing to act on."""
    bamdir = _bamdir(tmp_path, ["s0"], 60, split_older=True)

    with pytest.raises(toolsTG.InvalidParameterError) as raised:
        _refuse_stale_split_bams(str(bamdir), 60, ["s0"], logging.getLogger("t"))

    message = str(raised.value)
    assert "s0.bam" in message
    assert "--overwritesplits" in message
    assert "OLDER" in message


def test_a_build_refuses_to_reuse_a_stale_split(tmp_path):
    """End to end through _handle_preprocessing: this is the exact shape of the c305 incident,
    and the run must stop instead of counting the pre-dedup reads."""
    samples = ["s0", "s1"]
    bamdir = _bamdir(tmp_path, samples, 60, split_older=True)
    args = SimpleNamespace(
        input=_make_manifest(tmp_path, samples), bamdir=str(bamdir),
        output=str(tmp_path / "exp" / "exp.h5ad"), database=str(tmp_path / "db"),
        overwritebams=False, overwritesplits=False, readlengthsplit=60, threads=1,
    )

    with patch("trnagraph.modules.toolsSplit.BamSplitter.__init__", lambda self, args: None), \
         patch("trnagraph.modules.toolsSplit.BamSplitter.process", lambda self: None), \
         pytest.raises(toolsTG.InvalidParameterError):
        _bare_builder()._handle_preprocessing(args)


def test_overwritesplits_proceeds_without_complaining(tmp_path):
    """The escape hatch has to actually work on the stale case, since that is the only case it
    exists for. The offending files are about to be rewritten, so there is nothing to refuse."""
    samples = ["s0"]
    bamdir = _bamdir(tmp_path, samples, 60, split_older=True)
    args = SimpleNamespace(
        input=_make_manifest(tmp_path, samples), bamdir=str(bamdir),
        output=str(tmp_path / "exp" / "exp.h5ad"), database=str(tmp_path / "db"),
        overwritebams=False, overwritesplits=True, readlengthsplit=60, threads=1,
    )

    called = []
    with patch("trnagraph.modules.toolsSplit.BamSplitter.__init__", lambda self, a: None), \
         patch("trnagraph.modules.toolsSplit.BamSplitter.process",
               lambda self: called.append(True)):
        _bare_builder()._handle_preprocessing(args)

    assert called, "--overwritesplits must re-cut rather than skip"


def test_overwritebams_does_not_re_cut_the_splits(tmp_path):
    """The two flags are separate operations. Conflating them made "re-cut the splits" impossible
    to ask for without also re-mapping from FASTQ -- unreachable once the FASTQs have moved, which
    is exactly the situation a stale split is discovered in."""
    samples = ["s0"]
    bamdir = _bamdir(tmp_path, samples, 60, split_older=True)
    args = SimpleNamespace(
        input=_make_manifest(tmp_path, samples), bamdir=str(bamdir),
        output=str(tmp_path / "exp" / "exp.h5ad"), database=str(tmp_path / "db"),
        overwritebams=True, overwritesplits=False, readlengthsplit=60, threads=1,
    )

    with patch("trnagraph.modules.toolsMap.MapSamples.__init__", lambda self, a: None), \
         patch("trnagraph.modules.toolsMap.MapSamples.main", lambda self: None), \
         patch("trnagraph.modules.toolsSplit.BamSplitter.__init__", lambda self, a: None), \
         patch("trnagraph.modules.toolsSplit.BamSplitter.process", lambda self: None), \
         pytest.raises(toolsTG.InvalidParameterError):
        _bare_builder()._handle_preprocessing(args)


def test_a_missing_split_is_not_reported_as_stale(tmp_path):
    """Nothing to compare is not an error -- that path re-splits anyway."""
    bamdir = _bamdir(tmp_path, ["s0"], 60, split_older=True)
    os.remove(bamdir / "u60" / "s0.bam")

    stale = _stale_split_bams(str(bamdir), 60, ["s0"])

    assert {tag for tag, _, _ in stale} == {"o60"}
