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
from matplotlib.backends.backend_pdf import PdfPages

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
    # gtRNAdb-shaped names: _covobs_list() sorts `trna` values by anticodon then copy number and
    # parses that trailing integer, so a bare "tRNA0" is not a name this code ever sees.
    trnas = [f"tRNA-Ala-AGC-{i + 1}" for i in range(n_trnas) for _ in range(2)]
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


# ---------------------------------------------------------------------------------------------
# The other three coverage sub-steps. Measured on the hg38 object (436 tRNAs x 3 groups), a
# `-g coverage` run spent 6m21s of which 5m42s produced no output at all: generate_split is the
# only step that ever reported anything, so a run that was working normally was indistinguishable
# from one that had wedged. These pin that every step reports, and that the combined pages are
# written as they are rendered rather than all being held until the last one is drawn.


def _fake_combine_page(self, ulist, coverage_fill):
    """Stands in for the real matplotlib render; returns a real (empty) Figure because the
    caller hands it straight to PdfPages.savefig()."""
    import matplotlib.pyplot as plt
    return plt.figure()


_fake_combine_page.__name__ = "generate_combine_page"


def _fake_partition_frame(self):
    """Two tRNAs x two groups of position-by-category coverage, the shape the real one returns."""
    return {
        trna: {group: pd.DataFrame({'coverage': [1.0, 2.0]}) for group in ('A', 'B')}
        for trna in ('tRNA-Ala-AGC-1', 'tRNA-Ala-AGC-2')
    }


_fake_partition_frame.__name__ = "_partition_frame"


def test_generate_combine_reports_progress_and_advances_the_tracker(tmp_path):
    logger = logging.getLogger("trnagraph.modules.plotsCoverage")
    # Three tRNAs is a single page of the 16-per-page layout, rendered once per fill style.
    tracker = toolsTG.PhaseTracker(phases=["coverage"], logger=logger, desc="Graphing", weights=[2])
    viz = _make_visualizer(tmp_path, n_trnas=3, phase_tracker=tracker)
    viz.build_output_dirs()

    with patch.object(visualizer, "generate_combine_page", _fake_combine_page):
        with tracker.phase():
            viz.generate_combine()

    assert tracker._done_weight == 2, "expected one tick per combined page written"


def test_generate_combine_saves_each_page_as_it_is_rendered(tmp_path):
    """Streaming, not accumulation: the old implementation rendered every page into a list before
    saving any of them, so peak memory held all 28 pages of an hg38 run at once and no page could
    be counted until the last was drawn."""
    events = []
    viz = _make_visualizer(tmp_path, n_trnas=40)  # 40 tRNAs -> 3 pages per fill style
    viz.build_output_dirs()

    def _recording_page(self, ulist, coverage_fill):
        events.append('render')
        return _fake_combine_page(self, ulist, coverage_fill)
    _recording_page.__name__ = "generate_combine_page"

    real_savefig = PdfPages.savefig

    def _recording_savefig(self, *args, **kwargs):
        events.append('save')
        return real_savefig(self, *args, **kwargs)

    with patch.object(visualizer, "generate_combine_page", _recording_page), \
         patch.object(PdfPages, "savefig", _recording_savefig):
        viz.generate_combine()

    assert events.count('render') == 6 and events.count('save') == 6, events
    assert events == ['render', 'save'] * 6, (
        f'pages were not streamed; render/save order was {events}'
    )


def test_generate_partition_split_reports_progress_and_advances_the_tracker(tmp_path, caplog):
    logger = logging.getLogger("trnagraph.modules.plotsCoverage")
    # Two tRNAs x two groups = four individual specificity plots.
    tracker = toolsTG.PhaseTracker(phases=["coverage"], logger=logger, desc="Graphing", weights=[4])
    viz = _make_visualizer(tmp_path, n_trnas=2, phase_tracker=tracker)
    viz.build_output_dirs()

    with patch.object(visualizer, "_partition_frame", _fake_partition_frame), \
         patch.object(visualizer, "generate_partition_plot", lambda *a, **k: None), \
         caplog.at_level(logging.INFO, logger="trnagraph.modules.plotsCoverage"):
        with tracker.phase():
            viz.generate_partition_split()

    assert tracker._done_weight == 4, "expected one tick per specificity plot written"
    assert "Specificity plots" in "\n".join(r.message for r in caplog.records)


def _make_grid_visualizer(tmp_path, n_trnas, n_groups=2, phase_tracker=None):
    """A visualizer whose var carries the full read-specificity partition, which
    generate_partition_overview() requires before it will draw anything."""
    trnas = [f"tRNA-Ala-AGC-{i + 1}" for i in range(n_trnas) for _ in range(n_groups)]
    groups = [chr(ord('A') + g) for _ in range(n_trnas) for g in range(n_groups)]
    n_obs = len(trnas)
    obs = pd.DataFrame({"trna": trnas, "group": groups},
                       index=[f"obs{i}" for i in range(n_obs)])
    var = pd.DataFrame(
        {"gap": [False] * len(toolsTG.COVERAGE_PARTITION),
         "coverage": list(toolsTG.COVERAGE_PARTITION)},
        index=list(toolsTG.COVERAGE_PARTITION),
    )
    adata = ad.AnnData(X=np.zeros((n_obs, len(var))), obs=obs, var=var)
    return visualizer(
        adata, threads=1, coverage_grp="group", coverage_obs="trna", coverage_type="coverage",
        coverage_gap=[False], coverage_method="mean", colormap=None, output=str(tmp_path) + "/",
        phase_tracker=phase_tracker, quiet=False,
    )


def test_specificity_grid_saves_each_page_as_it_is_built(tmp_path):
    """The grid built every page into a list before saving any: on the human build that was 107
    pages held at once, and 56 of the step's 91 seconds passed with the progress bar still at
    zero, because nothing could be counted until the last page existed."""
    import matplotlib.pyplot as plt

    events = []
    viz = _make_grid_visualizer(tmp_path, n_trnas=9)  # 9 tRNAs / 4 rows per page -> 3 pages
    viz.build_output_dirs()

    def _recording_page(self, rows, columns, frames, all_groups, palette, labels):
        events.append('build')
        return plt.figure()
    _recording_page.__name__ = "_partition_page"

    real_savefig = PdfPages.savefig

    def _recording_savefig(self, *args, **kwargs):
        events.append('save')
        return real_savefig(self, *args, **kwargs)

    def _frames(self):
        return {trna: {group: pd.DataFrame({'coverage': [1.0, 2.0]}) for group in ('A', 'B')}
                for trna in viz._covobs_list()}
    _frames.__name__ = "_partition_frame"

    with patch.object(visualizer, "_partition_frame", _frames), \
         patch.object(visualizer, "_partition_page", _recording_page), \
         patch.object(PdfPages, "savefig", _recording_savefig):
        viz.generate_partition_overview()

    assert events == ['build', 'save'] * 3, (
        f'grid pages were not streamed; build/save order was {events}'
    )
