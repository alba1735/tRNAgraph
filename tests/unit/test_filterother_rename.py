"""Regression test for roadmap.md Phase 2's `--dumpother` -> `--filterother` rename (naming-
convention consistency with `--filtermultimapped`/`--filtermismatches`). No behavior change --
just confirms AnalysisPipeline reads `args.filterother` (not the old `args.dumpother`) and still
threads it through to counttypes()'s `bamnofeature` kwarg correctly."""
from types import SimpleNamespace
from unittest.mock import patch

from trnagraph.modules.adataBuild import AnalysisPipeline


def _make_args(tmp_path, filterother):
    return SimpleNamespace(
        database=str(tmp_path / "db"), output=str(tmp_path / "exp" / "exp.h5ad"),
        input=str(tmp_path / "metadata.txt"), gtf=None, bed=[],
        bamdir=str(tmp_path / "bam"), threads=1,
        minnontrnasize=20, maxmismatches=None, mincoverage=None, filtermultimapped=False,
        pairs=None, hubonly=False, hub=False, filterother=filterother, quiet=True,
    )


def test_analysispipeline_reads_filterother_attribute(tmp_path):
    pipeline = AnalysisPipeline(_make_args(tmp_path, filterother=True))
    assert pipeline.filterother is True
    assert not hasattr(pipeline, "dumpother")


def test_counttypes_threads_filterother_into_bamnofeature(tmp_path):
    pipeline = AnalysisPipeline(_make_args(tmp_path, filterother=True))

    captured = {}

    def _fake_main(**kwargs):
        captured.update(kwargs)

    with patch("trnagraph.modules.adataBuild.toolsCountReads.main", _fake_main):
        pipeline.counttypes()

    assert captured["bamnofeature"] is True
