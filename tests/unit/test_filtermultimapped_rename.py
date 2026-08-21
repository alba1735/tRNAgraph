"""Regression test for roadmap.md Phase 2's `--uniqueonly` -> `--filtermultimapped` rename. No
behavior change -- default stays False (matches tRAX's own --uniqueonlycov default), same
mechanism (drops genomically-multi-mapped reads, toolsGetCoverage.py's issinglemapped() check).
Just confirms AnalysisPipeline reads `args.filtermultimapped` (not the old `args.uniqueonly`) and
still threads it through to toolsGetCoverage.main()'s `uniqueonly` kwarg (that internal parameter
name is unchanged -- it's not directly CLI-facing, same scoping precedent as `--filterother`
leaving toolsCountReads.py's internal `bamnofeature` name untouched)."""
from types import SimpleNamespace
from unittest.mock import patch

from trnagraph.modules.adataBuild import AnalysisPipeline


def _make_args(tmp_path, filtermultimapped):
    return SimpleNamespace(
        database=str(tmp_path / "db"), output=str(tmp_path / "exp" / "exp.h5ad"),
        input=str(tmp_path / "metadata.txt"), gtf=None, bed=[],
        bamdir=str(tmp_path / "bam"), threads=1,
        minnontrnasize=20, maxmismatches=None,
        filtermultimapped=filtermultimapped, pairs=None, hubonly=False, hub=False,
        filterother=False, quiet=True,
    )


def test_analysispipeline_reads_filtermultimapped_attribute(tmp_path):
    pipeline = AnalysisPipeline(_make_args(tmp_path, filtermultimapped=True))
    assert pipeline.filtermultimapped is True
    assert not hasattr(pipeline, "uniqueonlycov")


def test_gettrnacoverage_threads_filtermultimapped_into_uniqueonly_kwarg(tmp_path):
    pipeline = AnalysisPipeline(_make_args(tmp_path, filtermultimapped=True))

    captured = {}

    def _fake_main(**kwargs):
        captured.update(kwargs)

    with patch("trnagraph.modules.adataBuild.toolsGetCoverage.main", _fake_main):
        pipeline.gettrnacoverage("euk")

    assert captured["uniqueonly"] is True
