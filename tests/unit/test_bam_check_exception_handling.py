"""Regression test for roadmap.md's "Persistant BAM check warning" bug. Previously,
AnnDataBuilder._handle_preprocessing() wrapped the metadata parse in a bare `except Exception:`
that discarded the real error and fell back to `samples = []` -- which then silently skipped the
per-sample BAM-existence check entirely (the loop over `samples` never executes), leaving only the
bam directory's own existence checked. A partially-populated bam directory would therefore be
silently treated as "BAM files found, skipping mapping". The fix logs the real exception and
re-raises instead of swallowing it and continuing with a degraded check."""
import logging
from types import SimpleNamespace

import pytest

from trnagraph.modules.adataBuild import AnnDataBuilder


def _bare_builder():
    builder = object.__new__(AnnDataBuilder)
    builder.logger = logging.getLogger("trnagraph.modules.adataBuild")
    builder._split_dirs_preexisted = {}
    return builder


def test_handle_preprocessing_raises_on_unparseable_metadata(tmp_path):
    """A metadata file that can't be parsed must abort the pipeline, not silently degrade to a
    directory-only BAM check."""
    bamdir = tmp_path / "bam"
    bamdir.mkdir()
    (bamdir / "sample0.bam").write_text("")  # a partially-populated bamdir, missing sample1.bam

    args = SimpleNamespace(
        input=str(tmp_path / "does_not_exist.tsv"), bamdir=str(bamdir),
        output=str(tmp_path / "exp" / "exp.h5ad"), overwritebams=False,
    )

    builder = _bare_builder()
    with pytest.raises(Exception):
        builder._handle_preprocessing(args)


def test_handle_preprocessing_logs_the_real_exception(tmp_path, caplog):
    """The log must surface the actual parse failure (e.g. the missing filename), not just the
    old generic "Could not parse metadata..." message with no detail."""
    missing_path = str(tmp_path / "does_not_exist.tsv")
    args = SimpleNamespace(
        input=missing_path, bamdir=str(tmp_path / "bam"),
        output=str(tmp_path / "exp" / "exp.h5ad"), overwritebams=False,
    )

    builder = _bare_builder()
    with caplog.at_level(logging.ERROR, logger="trnagraph.modules.adataBuild"):
        with pytest.raises(Exception):
            builder._handle_preprocessing(args)

    messages = "\n".join(r.message for r in caplog.records)
    assert missing_path in messages
