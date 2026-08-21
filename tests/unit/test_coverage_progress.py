"""Regression tests for plotsCoverage.py's progress reporting (roadmap.md Phase 2: "tqdm"
graphing objective). `generate_split()` previously drove its multiprocessing pool with a blocking
`Pool.map()` and no progress feedback -- the one graph type the roadmap explicitly calls out as
producing "hundreds/thousands" of plots, versus other types' "a handful". It now wraps its
`pool.imap_unordered()` in the shared `toolsTG.progress_iterator()` helper (same pattern as
toolsTrim.py/toolsMap.py/toolsCountReads.py) AND ticks an optional `toolsTG.PhaseTracker` once per
completed plot via `advance()`, so a graphing command's outer per-graph-type bar can move in
lockstep with coverage's own per-plot completions instead of only ticking atomically when the
whole "coverage" phase finishes. `generate_split_single` (the actual matplotlib rendering) is
mocked out here -- these tests exercise the progress wiring, not plot correctness."""
import logging
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd

from trnagraph.modules import toolsTG
from trnagraph.modules.plotsCoverage import visualizer


def _fake_generate_split_single(self, covobs):
    '''
    Replaces visualizer.generate_split_single for these tests -- must be a real, named function
    (never a lambda/mock) so `multiprocessing`'s bound-method pickling can reconstruct it inside a
    Pool worker: it serializes a bound method as (getattr, (obj, func.__name__)), keyed off the
    FUNCTION's own __name__, not the class attribute name it was patched onto. A lambda's
    __name__ is literally "<lambda>", which doesn't exist as an attribute on the instance --
    pickling then fails inside the worker process, which crashes and respawns in a loop the main
    process never observes as an error, hanging generate_split() forever instead of raising.
    '''
    return None


_fake_generate_split_single.__name__ = "generate_split_single"


def _make_adata(trnas, group_values):
    n_obs = len(trnas)
    obs = pd.DataFrame({"trna": trnas, "group": group_values}, index=[f"obs{i}" for i in range(n_obs)])
    var = pd.DataFrame({"gap": [False, False], "coverage": ["coverage", "coverage"]}, index=["pos0", "pos1"])
    return ad.AnnData(X=np.zeros((n_obs, 2)), obs=obs, var=var)


def _make_visualizer(tmp_path, n_trnas=3, quiet=False, phase_tracker=None):
    trnas = [f"tRNA{i}" for i in range(n_trnas) for _ in range(2)]
    groups = ["A", "B"] * n_trnas
    adata = _make_adata(trnas, groups)
    return visualizer(
        adata, threads=1, coverage_grp="group", coverage_obs="trna", coverage_type="coverage",
        coverage_gap=[False], coverage_method="mean", colormap=None, output=str(tmp_path) + "/",
        phase_tracker=phase_tracker, quiet=quiet,
    )


def test_generate_split_logs_progress_milestones(tmp_path, caplog):
    viz = _make_visualizer(tmp_path, n_trnas=3)

    with patch.object(visualizer, "generate_split_single", _fake_generate_split_single), \
         caplog.at_level(logging.INFO, logger="trnagraph.modules.plotsCoverage"):
        viz.generate_split()

    messages = "\n".join(r.message for r in caplog.records)
    assert "Coverage plots" in messages
    assert "100%" in messages


def test_generate_split_advances_a_shared_phase_tracker_per_plot(tmp_path):
    logger = logging.getLogger("trnagraph.modules.plotsCoverage")
    tracker = toolsTG.PhaseTracker(phases=["coverage"], logger=logger, desc="Graphing", weights=[3])
    viz = _make_visualizer(tmp_path, n_trnas=3, phase_tracker=tracker)

    with patch.object(visualizer, "generate_split_single", _fake_generate_split_single):
        with tracker.phase():
            viz.generate_split()

    assert tracker._done_weight == 3


def test_generate_split_without_a_phase_tracker_still_works(tmp_path):
    """phase_tracker defaults to None -- plotsCoverage.py must remain usable standalone (e.g. if
    ever called outside adataGraph.py's orchestration) without crashing on a missing tracker."""
    viz = _make_visualizer(tmp_path, n_trnas=2, phase_tracker=None)

    with patch.object(visualizer, "generate_split_single", _fake_generate_split_single):
        viz.generate_split()  # must not raise
