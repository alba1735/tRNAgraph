"""Regression tests for roadmap.md Phase 2's "BAM deletion safety" item. Previously,
AnnDataBuilder._apply_readlength_split_() (and the standalone `add_split()` command) deleted a
u{N}/o{N} split-BAM directory whenever --savesplitbams wasn't passed on *this* run, with no check
for whether this run actually created it -- a real, reproducible data-loss bug: a prior run with
--savesplitbams keeps u60/o60 on disk; a later run without that flag sees the files already exist
(toolsSplit.BamSplitter already skips re-splitting in that case) but deletes them anyway at the
end. The fix records, in _handle_preprocessing() (and inline in add_split()), whether each tag's
directory was already fully populated *before* this run's own splitting step -- via the shared
`_split_bam_dirs_preexisting()` helper -- so the later cleanup only ever deletes what this run
actually created or explicitly regenerated (--overwritebams)."""
import logging
from types import SimpleNamespace
from unittest.mock import patch

from trnagraph.modules.adataBuild import AnnDataBuilder, _split_bam_dirs_preexisting


def _bare_builder():
    """A minimally-initialized AnnDataBuilder, bypassing the full constructor's metadata/AnnData
    setup -- _handle_preprocessing() only needs self.logger and self._split_dirs_preexisted."""
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


def _make_fully_populated_bamdir(tmp_path, samplenames, cutoff):
    bamdir = tmp_path / "bam"
    bamdir.mkdir()
    for name in samplenames:
        (bamdir / f"{name}.bam").write_text("")
    for tag in (f"u{cutoff}", f"o{cutoff}"):
        tagdir = bamdir / tag
        tagdir.mkdir()
        for name in samplenames:
            (tagdir / f"{name}.bam").write_text("")
    return bamdir


def test_handle_preprocessing_marks_fully_populated_split_dirs_as_preexisting_when_not_overwriting(tmp_path):
    samplenames = ["sample0", "sample1"]
    bamdir = _make_fully_populated_bamdir(tmp_path, samplenames, cutoff=60)
    manifest = _make_manifest(tmp_path, samplenames)

    args = SimpleNamespace(
        input=manifest, bamdir=str(bamdir), output=str(tmp_path / "exp" / "exp.h5ad"),
        overwritebams=False, readlengthsplit=60, threads=1,
    )

    builder = _bare_builder()
    builder._handle_preprocessing(args)

    assert builder._split_dirs_preexisted == {"u60": True, "o60": True}
    # Sanity: this check itself must not have touched anything.
    for tag in ("u60", "o60"):
        for name in samplenames:
            assert (bamdir / tag / f"{name}.bam").exists()


def test_handle_preprocessing_does_not_mark_preexisting_when_overwritesplits_forces_a_fresh_split(tmp_path):
    """--overwritesplits means this run regenerates the split BAMs regardless of what was already
    there, so that tag's output is this run's own to clean up afterward -- unlike the untouched
    case above, it should NOT be protected from --savesplitbams-off cleanup."""
    samplenames = ["sample0", "sample1"]
    bamdir = _make_fully_populated_bamdir(tmp_path, samplenames, cutoff=60)
    manifest = _make_manifest(tmp_path, samplenames)

    args = SimpleNamespace(
        input=manifest, bamdir=str(bamdir), output=str(tmp_path / "exp" / "exp.h5ad"),
        database=str(tmp_path / "db"), overwritesplits=True, readlengthsplit=60, threads=1,
    )

    builder = _bare_builder()
    with patch("trnagraph.modules.toolsMap.MapSamples.__init__", lambda self, args: None), \
         patch("trnagraph.modules.toolsMap.MapSamples.main", lambda self: None), \
         patch("trnagraph.modules.toolsSplit.BamSplitter.__init__", lambda self, args: None), \
         patch("trnagraph.modules.toolsSplit.BamSplitter.process", lambda self: None):
        builder._handle_preprocessing(args)

    assert builder._split_dirs_preexisted == {"u60": False, "o60": False}


def test_split_bam_dirs_preexisting_requires_every_sample(tmp_path):
    samplenames = ["sample0", "sample1"]
    bamdir = tmp_path / "bam"
    bamdir.mkdir()
    (bamdir / "u60").mkdir()
    (bamdir / "u60" / "sample0.bam").write_text("")
    # sample1.bam missing from u60 -- not fully populated.
    (bamdir / "o60").mkdir()
    for name in samplenames:
        (bamdir / "o60" / f"{name}.bam").write_text("")

    result = _split_bam_dirs_preexisting(str(bamdir), 60, samplenames)

    assert result == {"u60": False, "o60": True}
