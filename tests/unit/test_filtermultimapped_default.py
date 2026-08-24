"""Coverage always excludes genomically multi-mapped reads, with no flag to change it.

This matches tRAX, which does the same thing unconditionally: `getcoverage.py`
parses a `uniqueonly` setting at line 666 and then passes the literal
`uniqueonly = True` at both transcript-coverage call sites (760 and 766), so its
own `--uniqueonly` flag has no effect and its coverage is always restricted to
MAPQ >= 2 alignments.

tRNAgraph's original port (`09da375`) reproduced that faithfully, hardcoding True
in the same two places. Commit `dca2673` -- a module refactor introducing
adataBuild/adataCluster/adataGraph/adataMerge -- replaced those literals with
`uniqueonly=uniqueonly`, which reads as wiring up a dead parameter but silently
changed what every coverage file contained, because the CLI flag defaulted to
off. On the human hg38 dataset that made tRNAgraph's coverage 1.28x tRAX's,
scaling with how multi-mapped each feature is.

No flag is offered to restore the unfiltered behaviour. It would add no
information -- the AnnData object already carries both unique and total counts
(`nreads_*_unique_raw`/`_norm` alongside `nreads_*_raw`/`_norm`) and the coverage
table keeps its `uniquecoverage`/`multitrnacoverage`/`multianticodoncoverage`/
`multiaminocoverage` breakdown regardless -- while reintroducing a divergence from
tRAX. These tests exist so the hardcoded value is not "tidied up" into a
parameter again.
"""
import inspect
from types import SimpleNamespace
from unittest.mock import patch

from trnagraph.modules.adataBuild import AnalysisPipeline


def test_cli_build_exposes_no_multimapping_flag():
    from trnagraph import cli
    params = inspect.signature(cli.build).parameters
    assert "keepmultimapped" not in params
    assert "filtermultimapped" not in params


def _make_args(tmp_path):
    return SimpleNamespace(
        database=str(tmp_path / "db"), output=str(tmp_path / "exp" / "exp.h5ad"),
        input=str(tmp_path / "metadata.txt"), gtf=None, bed=[],
        bamdir=str(tmp_path / "bam"), threads=1,
        minnontrnasize=20, maxmismatches=None,
        pairs=None, hubonly=False, hub=False,
        filterother=False, quiet=True,
    )


def test_coverage_always_filters_multimapped(tmp_path):
    """Not merely defaulted to True -- unconditional, and not read from args."""
    pipeline = AnalysisPipeline(_make_args(tmp_path))
    captured = {}
    with patch("trnagraph.modules.adataBuild.toolsGetCoverage.main",
               lambda **kw: captured.update(kw)):
        pipeline.gettrnacoverage("euk")
    assert captured["uniqueonly"] is True


def test_pipeline_keeps_no_multimapping_attribute(tmp_path):
    pipeline = AnalysisPipeline(_make_args(tmp_path))
    assert not hasattr(pipeline, "filtermultimapped")
    assert not hasattr(pipeline, "uniqueonlycov")
