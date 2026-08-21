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

from trnagraph.modules.toolsTestSuite import _LiveBoxHandler, _friendly_variant_title


def _emit(handler, message):
    record = logging.LogRecord(name="test", level=logging.INFO, pathname="", lineno=0, msg=message, args=None, exc_info=None)
    handler.handle(record)


def test_bare_item_milestone_drives_the_bar_by_default():
    milestones = []
    handler = _LiveBoxHandler(tail=[], on_milestone=lambda c, t, label=None: milestones.append((c, t)))
    handler.setFormatter(logging.Formatter('%(message)s'))

    _emit(handler, "Counting reads: 7/10 (70%) complete")

    assert milestones == [(7, 10)]


def test_phase_milestone_drives_the_bar():
    milestones = []
    handler = _LiveBoxHandler(tail=[], on_milestone=lambda c, t, label=None: milestones.append((c, t)))
    handler.setFormatter(logging.Formatter('%(message)s'))

    _emit(handler, "Build phase 2/6 (33%) complete: Analyzing counts")

    assert milestones == [(2, 6)]


def test_phase_advance_milestone_without_the_word_complete_also_drives_the_bar():
    """toolsTG.PhaseTracker.advance()'s intermediate milestones (e.g. one per coverage plot
    within a large weighted phase) omit the word "complete" -- only the phase's own final
    completion line includes it -- so the box must recognize both forms."""
    milestones = []
    handler = _LiveBoxHandler(tail=[], on_milestone=lambda c, t, label=None: milestones.append((c, t)))
    handler.setFormatter(logging.Formatter('%(message)s'))

    _emit(handler, "Graphing phase 3/10 (42%): coverage")

    assert milestones == [(3, 10)]


def test_phase_only_ignores_bare_item_milestones_even_before_any_phase_signal():
    """The corrected regression fix: with phase_only=True, a per-sample counting milestone (bare
    format) firing during the build's first phase must never drive the bar at all -- not even
    before the first phase-completion line arrives (which would otherwise let the bar climb
    through the inner loop's progress and then jump back down once the real, lower outer
    percentage for that phase arrives)."""
    milestones = []
    handler = _LiveBoxHandler(tail=[], on_milestone=lambda c, t, label=None: milestones.append((c, t)), phase_only=True)
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
    handler = _LiveBoxHandler(tail=[], on_milestone=lambda c, t, label=None: milestones.append((c, t)))
    handler.setFormatter(logging.Formatter('%(message)s'))

    _emit(handler, "Trimming samples: 3/10 (30%) complete")
    _emit(handler, "Trimming samples: 10/10 (100%) complete")

    assert milestones == [(3, 10), (10, 10)]


def test_every_message_still_reaches_the_tail_regardless_of_milestone_type_or_phase_only():
    tail = []
    handler = _LiveBoxHandler(tail=tail, on_milestone=lambda c, t, label=None: None, phase_only=True)
    handler.setFormatter(logging.Formatter('%(message)s'))

    _emit(handler, "Build phase 1/6 (17%) complete: Counting Reads")
    _emit(handler, "Counting reads: 10/10 (100%) complete")
    _emit(handler, "some other plain log line")

    assert tail == [
        "Build phase 1/6 (17%) complete: Counting Reads",
        "Counting reads: 10/10 (100%) complete",
        "some other plain log line",
    ]


def test_phase_milestone_passes_the_label_through_to_on_milestone():
    """The box's title needs the phase's label (e.g. to notice a split-variant bracket prefix and
    switch the displayed title) -- _LiveBoxHandler must extract it from the trailing ": <label>"
    portion of a phase-tracker message and pass it through, not just the completed/total counts."""
    seen = []
    handler = _LiveBoxHandler(tail=[], on_milestone=lambda c, t, label=None: seen.append(label))
    handler.setFormatter(logging.Formatter('%(message)s'))

    _emit(handler, "Build phase 8/18 (44%) complete: [Under 60] Counting Reads")

    assert seen == ["[Under 60] Counting Reads"]


def test_bare_item_milestone_passes_no_label():
    seen = []
    handler = _LiveBoxHandler(tail=[], on_milestone=lambda c, t, label=None: seen.append(label))
    handler.setFormatter(logging.Formatter('%(message)s'))

    _emit(handler, "Counting reads: 7/10 (70%) complete")

    assert seen == [None]


def test_friendly_variant_title_maps_under_to_fragments():
    assert _friendly_variant_title("[Under 60] Counting Reads") == "Building under 60bp split (Fragments)..."


def test_friendly_variant_title_maps_over_to_full_length():
    assert _friendly_variant_title("[Over 60] Generating Read Coverage plots") == "Building over 60bp split (Full-length)..."


def test_friendly_variant_title_returns_none_for_a_non_variant_label():
    """Non-variant phase labels (the main/full build) must not override the box's title -- it
    should stay whatever _live_box() was originally given (e.g. "Building AnnData object...")."""
    assert _friendly_variant_title("Counting Reads") is None
    assert _friendly_variant_title("Building AnnData object") is None
