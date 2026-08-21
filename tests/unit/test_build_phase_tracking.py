"""Regression tests for AnalysisPipeline.run()'s phase-tracking wiring (roadmap.md Phase 2:
"tqdm" -- post-Stage-3 design follow-up fixing `analyze build`'s progress reporting, which
previously hit 100% early via toolsCountReads.py's per-sample counting milestones and stayed
there through DESeq2 fitting/coverage/VST/writing the h5ad). These tests mock out the heavy
per-phase work (countfeatures/run_deseq2/counttypes/run_unique_deseq2/gettrnacoverage) to isolate
the wiring itself: that AnalysisPipeline.run() advances a toolsTG.PhaseTracker once per phase, in
order, folding in a variant label when one is passed (the split-build case)."""
import logging
from types import SimpleNamespace
from unittest.mock import patch

from trnagraph.modules import toolsTG
from trnagraph.modules.adataBuild import AnalysisPipeline, _analysis_pipeline_phase_names


def _make_args(tmp_path, nosizefactors=False):
    return SimpleNamespace(
        database=str(tmp_path / "db"), output=str(tmp_path / "exp" / "exp.h5ad"),
        input=str(tmp_path / "metadata.txt"), gtf=None, bed=[], nofrag=False,
        nosizefactors=nosizefactors, bamdir=str(tmp_path / "bam"), threads=1,
        minnontrnasize=20, maxmismatches=None, mincoverage=None, filtermultimapped=False,
        pairs=None, hubonly=False, hub=False, filterother=False, quiet=True,
    )


def _patched_pipeline_methods():
    return (
        patch.object(AnalysisPipeline, "makefeaturebed", lambda self: None),
        patch.object(AnalysisPipeline, "countfeatures", lambda self: None),
        patch.object(AnalysisPipeline, "run_deseq2", lambda self: None),
        patch.object(AnalysisPipeline, "counttypes", lambda self: None),
        patch.object(AnalysisPipeline, "run_unique_deseq2", lambda self: None),
        patch.object(AnalysisPipeline, "gettrnacoverage", lambda self, orgtype: None),
        patch.object(AnalysisPipeline, "write_runinfo", lambda self: None),
    )


def test_run_advances_the_default_tracker_once_per_phase_in_order(tmp_path, caplog):
    args = _make_args(tmp_path)
    pipeline = AnalysisPipeline(args)

    patches = _patched_pipeline_methods()
    with patches[0], patches[1], patches[2], patches[3], patches[4], patches[5], patches[6], \
         caplog.at_level(logging.INFO, logger="trnagraph.modules.adataBuild"):
        pipeline.run()

    phase_messages = [r.message for r in caplog.records if r.message.startswith("Build phase")]
    assert phase_messages == [
        "Build phase 1/5 (20%) complete: Counting Reads",
        "Build phase 2/5 (40%) complete: Analyzing counts",
        "Build phase 3/5 (60%) complete: Counting Read Types",
        "Build phase 4/5 (80%) complete: Analyzing unique counts",
        "Build phase 5/5 (100%) complete: Generating Read Coverage plots",
    ]


def test_run_skips_deseq2_phases_when_nosizefactors(tmp_path, caplog):
    args = _make_args(tmp_path, nosizefactors=True)
    pipeline = AnalysisPipeline(args)

    patches = _patched_pipeline_methods()
    with patches[0], patches[1], patches[2], patches[3], patches[4], patches[5], patches[6], \
         caplog.at_level(logging.INFO, logger="trnagraph.modules.adataBuild"):
        pipeline.run()

    phase_messages = [r.message for r in caplog.records if r.message.startswith("Build phase")]
    assert phase_messages == [
        "Build phase 1/3 (33%) complete: Counting Reads",
        "Build phase 2/3 (67%) complete: Counting Read Types",
        "Build phase 3/3 (100%) complete: Generating Read Coverage plots",
    ]


def test_run_folds_a_variant_label_into_a_shared_tracker(tmp_path, caplog):
    """Regression test for the split-build case: a shared tracker (spanning more than just this
    instance's own phases) plus a variant_label should fold the label into each phase message,
    and this instance must NOT construct its own separate default tracker."""
    logger = logging.getLogger("trnagraph.modules.adataBuild")
    shared_tracker = toolsTG.PhaseTracker(
        phases=_analysis_pipeline_phase_names(nosizefactors=False) * 2, logger=logger, desc="Build",
    )
    args = _make_args(tmp_path)
    pipeline = AnalysisPipeline(args, phase_tracker=shared_tracker, variant_label="Under 60")
    assert pipeline.phase_tracker is shared_tracker

    patches = _patched_pipeline_methods()
    with patches[0], patches[1], patches[2], patches[3], patches[4], patches[5], patches[6], \
         caplog.at_level(logging.INFO, logger="trnagraph.modules.adataBuild"):
        pipeline.run()

    phase_messages = [r.message for r in caplog.records if r.message.startswith("Build phase")]
    assert phase_messages == [
        "Build phase 1/10 (10%) complete: [Under 60] Counting Reads",
        "Build phase 2/10 (20%) complete: [Under 60] Analyzing counts",
        "Build phase 3/10 (30%) complete: [Under 60] Counting Read Types",
        "Build phase 4/10 (40%) complete: [Under 60] Analyzing unique counts",
        "Build phase 5/10 (50%) complete: [Under 60] Generating Read Coverage plots",
    ]
    # Only 5 of the shared tracker's 10 phases were consumed by this one instance's run().
    assert shared_tracker._index == 5
