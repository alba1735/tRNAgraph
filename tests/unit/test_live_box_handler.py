"""Regression tests for toolsTestSuite.py's _LiveBoxHandler (roadmap.md Phase 2: "tqdm" --
post-Stage-3 design follow-up). The demo pipeline's box was hitting 100% early and staying
there: it drove its bar off the LAST log line matching progress_iterator's bare per-item
milestone format ("N/Total (P%) complete"), regardless of which sub-step emitted it -- so
toolsCountReads.py's per-sample counting milestones (finishing almost immediately relative to
the rest of a build) pinned the box at 100% through DESeq2 fitting/coverage/VST/writing the
h5ad. _LiveBoxHandler now recognizes toolsTG.PhaseTracker's distinct outer "<desc> phase N/Total
(P%) complete: <label>" format (note the literal word "phase" before the fraction).

`phase_only=True` (set upfront by whichever _live_box() call wraps a phase-tracked command, e.g.
"Building AnnData object...") makes bare per-item milestones never drive the bar at all -- this
has to be an upfront flag, not inferred reactively from whichever milestone type arrives first:
a phase-tracked command's FIRST phase can itself wrap an inner per-item loop (e.g. "Counting
Reads"), whose own milestones would otherwise drive the bar up to 100% before that phase's own
(much lower) outer percentage arrives, then visibly drop back down -- a glitch. `phase_only=False`
(the default, used by every other step -- trim/map/etc, which have no phase concept) leaves the
original behavior fully intact."""
import logging

from trnagraph.modules.toolsTestSuite import _LiveBoxHandler


def _emit(handler, message):
    record = logging.LogRecord(name="test", level=logging.INFO, pathname="", lineno=0, msg=message, args=None, exc_info=None)
    handler.handle(record)


def test_bare_item_milestone_drives_the_bar_by_default():
    milestones = []
    handler = _LiveBoxHandler(tail=[], on_milestone=lambda c, t: milestones.append((c, t)))
    handler.setFormatter(logging.Formatter('%(message)s'))

    _emit(handler, "Counting reads: 7/10 (70%) complete")

    assert milestones == [(7, 10)]


def test_phase_milestone_drives_the_bar():
    milestones = []
    handler = _LiveBoxHandler(tail=[], on_milestone=lambda c, t: milestones.append((c, t)))
    handler.setFormatter(logging.Formatter('%(message)s'))

    _emit(handler, "Build phase 2/6 (33%) complete: Analyzing counts")

    assert milestones == [(2, 6)]


def test_phase_only_ignores_bare_item_milestones_even_before_any_phase_signal():
    """The corrected regression fix: with phase_only=True, a per-sample counting milestone (bare
    format) firing during the build's first phase must never drive the bar at all -- not even
    before the first phase-completion line arrives (which would otherwise let the bar climb
    through the inner loop's progress and then jump back down once the real, lower outer
    percentage for that phase arrives)."""
    milestones = []
    handler = _LiveBoxHandler(tail=[], on_milestone=lambda c, t: milestones.append((c, t)), phase_only=True)
    handler.setFormatter(logging.Formatter('%(message)s'))

    _emit(handler, "Counting reads: 10/10 (100%) complete")  # inner milestone, phase 1 in progress
    _emit(handler, "Build phase 1/6 (17%) complete: Counting Reads")
    _emit(handler, "Counting reads: 10/10 (100%) complete")  # inner milestone, phase 2 in progress
    _emit(handler, "Build phase 2/6 (33%) complete: Analyzing counts")

    assert milestones == [(1, 6), (2, 6)]


def test_phase_only_false_matches_the_original_unconditional_bare_milestone_behavior():
    """Every step other than a phase-tracked one (trim/map/etc) passes phase_only=False (the
    default) and must behave exactly as before this change."""
    milestones = []
    handler = _LiveBoxHandler(tail=[], on_milestone=lambda c, t: milestones.append((c, t)))
    handler.setFormatter(logging.Formatter('%(message)s'))

    _emit(handler, "Trimming samples: 3/10 (30%) complete")
    _emit(handler, "Trimming samples: 10/10 (100%) complete")

    assert milestones == [(3, 10), (10, 10)]


def test_every_message_still_reaches_the_tail_regardless_of_milestone_type_or_phase_only():
    tail = []
    handler = _LiveBoxHandler(tail=tail, on_milestone=lambda c, t: None, phase_only=True)
    handler.setFormatter(logging.Formatter('%(message)s'))

    _emit(handler, "Build phase 1/6 (17%) complete: Counting Reads")
    _emit(handler, "Counting reads: 10/10 (100%) complete")
    _emit(handler, "some other plain log line")

    assert tail == [
        "Build phase 1/6 (17%) complete: Counting Reads",
        "Counting reads: 10/10 (100%) complete",
        "some other plain log line",
    ]
