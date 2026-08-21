"""Regression tests for toolsTG.PhaseTracker (roadmap.md Phase 2: "tqdm" -- follow-up design
discussion after Stage 3 shipped). `analyze build`'s `tools test` box was hitting 100% early and
staying there: toolsTestSuite.py's _LiveBoxHandler drives its bar off the LAST log line matching
progress_iterator's per-item milestone format, and Stage 3 wired toolsCountReads.py to emit that
format for per-sample counting, which finishes almost immediately relative to the rest of a build
(DESeq2 fitting, coverage generation, VST, writing the h5ad). PhaseTracker adds a coarser, outer
"phase" layer on top of a fixed named sequence of build steps -- each phase weighted (default
equal weight), advancing only once per phase (not per inner item), and logged in a message format
deliberately distinct from progress_iterator's bare per-item one so toolsTestSuite.py can tell the
two apart and never let an inner per-sample milestone override the outer phase-level percentage."""
import logging

import pytest

from trnagraph.modules import toolsTG


def test_equal_weight_phases_report_correct_cumulative_percentage(caplog):
    logger = logging.getLogger("test_phase_tracker.equal_weight")
    tracker = toolsTG.PhaseTracker(phases=["Counting Reads", "Analyzing counts", "Writing h5ad"], logger=logger, desc="Build")

    with caplog.at_level(logging.INFO, logger=logger.name):
        with tracker.phase():
            pass
        with tracker.phase():
            pass
        with tracker.phase():
            pass

    messages = [r.message for r in caplog.records]
    assert messages == [
        "Build phase 1/3 (33%) complete: Counting Reads",
        "Build phase 2/3 (67%) complete: Analyzing counts",
        "Build phase 3/3 (100%) complete: Writing h5ad",
    ]


def test_weighted_phases_report_percentage_proportional_to_weight(caplog):
    logger = logging.getLogger("test_phase_tracker.weighted")
    tracker = toolsTG.PhaseTracker(
        phases=["coverage", "pca"], logger=logger, desc="Graphing", weights=[90, 10],
    )

    with caplog.at_level(logging.INFO, logger=logger.name):
        with tracker.phase():
            pass
        with tracker.phase():
            pass

    messages = [r.message for r in caplog.records]
    assert messages == [
        "Graphing phase 1/2 (90%) complete: coverage",
        "Graphing phase 2/2 (100%) complete: pca",
    ]


def test_variant_label_is_folded_into_the_phase_message(caplog):
    logger = logging.getLogger("test_phase_tracker.variant")
    tracker = toolsTG.PhaseTracker(phases=["Counting Reads"], logger=logger, desc="Build")

    with caplog.at_level(logging.INFO, logger=logger.name):
        with tracker.phase(variant="Under60"):
            pass

    assert caplog.records[0].message == "Build phase 1/1 (100%) complete: [Under60] Counting Reads"


def test_calling_phase_more_times_than_declared_raises(caplog):
    logger = logging.getLogger("test_phase_tracker.overflow")
    tracker = toolsTG.PhaseTracker(phases=["only one phase"], logger=logger, desc="Build")

    with tracker.phase():
        pass

    with pytest.raises(IndexError):
        with tracker.phase():
            pass


def test_exception_inside_a_phase_does_not_log_a_false_completion(caplog):
    logger = logging.getLogger("test_phase_tracker.exception")
    tracker = toolsTG.PhaseTracker(phases=["Counting Reads", "Writing h5ad"], logger=logger, desc="Build")

    with caplog.at_level(logging.INFO, logger=logger.name):
        with pytest.raises(RuntimeError):
            with tracker.phase():
                raise RuntimeError("boom")

    assert caplog.records == []


def test_advance_ticks_the_outer_percentage_in_lockstep_with_inner_items():
    """Regression coverage for the graphing tqdm objective's design: a phase whose weight equals
    its own inner item count (e.g. a graph-type weighted by its expected plot count) should tick
    the outer percentage per completed item via advance(), rather than only completing atomically
    when the `with tracker.phase():` block exits."""
    logger = logging.getLogger("test_phase_tracker.advance")
    tracker = toolsTG.PhaseTracker(phases=["coverage", "pca"], logger=logger, desc="Graphing", weights=[4, 1])

    with tracker.phase():
        tracker.advance(1)
        assert tracker._done_weight == 1
        tracker.advance(1)
        tracker.advance(1)
        tracker.advance(1)
        assert tracker._done_weight == 4

    with tracker.phase():
        pass

    assert tracker._done_weight == 5


def test_advance_partial_progress_is_topped_up_atomically_when_phase_exits(caplog):
    """A phase that only partially ticks via advance() (e.g. it errors out after 2 of 4 expected
    plots, but still completes the `with` block normally) should still have its full weight
    counted once the block exits -- advance() is an optional finer-grained signal, not a
    replacement for the phase's own weight accounting."""
    logger = logging.getLogger("test_phase_tracker.advance_partial")
    tracker = toolsTG.PhaseTracker(phases=["coverage", "pca"], logger=logger, desc="Graphing", weights=[4, 1])

    with caplog.at_level(logging.INFO, logger=logger.name):
        with tracker.phase():
            tracker.advance(1)
            tracker.advance(1)
        with tracker.phase():
            pass

    completion_messages = [r.message for r in caplog.records if "complete" in r.message]
    assert completion_messages == [
        "Graphing phase 1/2 (80%) complete: coverage",
        "Graphing phase 2/2 (100%) complete: pca",
    ]


def test_advance_logs_a_throttled_milestone_at_ten_percent_steps_within_a_large_phase(caplog):
    """A phase weighted by a large plot count (e.g. coverage's hundreds/thousands) must not log
    on every single advance() call -- mirrors progress_iterator's own ~10%-of-total throttling so
    a long-running phase still gets periodic outer-percentage feedback without flooding the log."""
    logger = logging.getLogger("test_phase_tracker.advance_throttled")
    tracker = toolsTG.PhaseTracker(phases=["coverage"], logger=logger, desc="Graphing", weights=[20])

    with caplog.at_level(logging.INFO, logger=logger.name):
        with tracker.phase():
            for _ in range(20):
                tracker.advance(1)

    progress_messages = [r.message for r in caplog.records if "complete" not in r.message]
    assert len(progress_messages) == 10
    assert progress_messages[0] == "Graphing phase 1/1 (10%): coverage"
    assert progress_messages[-1] == "Graphing phase 1/1 (100%): coverage"


def test_advance_outside_a_phase_raises():
    logger = logging.getLogger("test_phase_tracker.advance_outside")
    tracker = toolsTG.PhaseTracker(phases=["coverage"], logger=logger, desc="Graphing")

    with pytest.raises(RuntimeError):
        tracker.advance(1)


def test_message_format_is_distinguishable_from_progress_iterator_bare_milestone():
    """toolsTestSuite.py's _LiveBoxHandler needs to tell PhaseTracker's outer messages apart from
    progress_iterator's inner per-item ones -- the literal word "phase" immediately before the
    fraction is the distinguishing marker; progress_iterator's own bare format never contains it."""
    logger = logging.getLogger("test_phase_tracker.format")
    tracker = toolsTG.PhaseTracker(phases=["Counting Reads"], logger=logger, desc="Build")

    import re
    bare_item_re = re.compile(r'(\d+)/(\d+) \(\d+%\) complete')
    phase_re = re.compile(r'\bphase (\d+)/(\d+) \((\d+)%\) complete\b')

    messages = []
    class _Capture(logging.Handler):
        def emit(self, record):
            messages.append(record.getMessage())
    logger.addHandler(_Capture())
    logger.setLevel(logging.INFO)

    with tracker.phase():
        pass

    assert len(messages) == 1
    assert phase_re.search(messages[0])
    # It also happens to satisfy the old bare regex (a substring match) -- that's fine and
    # expected, since toolsTestSuite.py checks the more specific phase_re FIRST.
    assert bare_item_re.search(messages[0])
