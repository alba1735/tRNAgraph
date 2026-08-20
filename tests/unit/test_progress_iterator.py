"""Regression tests for toolsTG.progress_iterator() (roadmap.md Phase 2: "tqdm" + the
file-logging architecture) -- the shared progress-reporting helper for long-running per-item
loops (trimming/mapping/counting/graphing), designed not to "poison" the persisted .log/ file or
a non-interactive/background run: a live rich bar+tail only on a real interactive terminal,
percentage-milestone log lines everywhere else, nothing at all under --quiet. File persistence
(propagation up to the centrally-configured 'trnagraph' logger's FileHandler) is unconditional
and independent of which console mode is active -- only the console-tagged handler is ever
temporarily detached, never the FileHandler, never propagation itself."""
import logging

import pytest

from trnagraph.modules import toolsTG


@pytest.fixture(autouse=True)
def _reset_trnagraph_logger():
    """Several tests attach mock handlers to the shared 'trnagraph' logger to exercise
    progress_iterator's real console-handler-detection mechanism -- reset it after each test."""
    yield
    logger = logging.getLogger('trnagraph')
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)


def test_quiet_still_emits_milestones_for_file_persistence(caplog):
    """
    --quiet must not skip milestone logging: file persistence (the whole point of the .log/
    architecture) is unconditional. In production, configure_logging() simply omits the console
    handler under --quiet, so these calls naturally never reach a real terminal -- but they must
    still reach the FileHandler, which this test's caplog handler stands in for.
    """
    logger = logging.getLogger("test_progress_iterator.quiet")

    with caplog.at_level(logging.INFO, logger=logger.name):
        result = list(toolsTG.progress_iterator(
            range(5), total=5, desc="Testing", logger=logger,
            quiet=True, isatty_fn=lambda: False,
        ))

    assert result == [0, 1, 2, 3, 4]
    assert len(caplog.records) > 0


def test_quiet_never_starts_the_interactive_display_even_on_a_real_tty():
    """A live rich Live session writes straight to the real terminal via stderr, bypassing the
    logging system (and thus --quiet) entirely -- quiet must prevent it from ever starting, even
    when isatty_fn reports a real terminal."""
    logger = logging.getLogger("test_progress_iterator.quiet_tty")

    result = list(toolsTG.progress_iterator(
        range(3), total=3, desc="Testing", logger=logger,
        quiet=True, isatty_fn=lambda: True,
    ))

    tail_handlers = [h for h in logger.handlers if isinstance(h, toolsTG._TailCaptureHandler)]
    assert result == [0, 1, 2]
    assert tail_handlers == []


def test_non_tty_emits_ten_percent_milestone_messages(caplog):
    logger = logging.getLogger("test_progress_iterator.non_tty_milestones")

    with caplog.at_level(logging.INFO, logger=logger.name):
        result = list(toolsTG.progress_iterator(
            range(20), total=20, desc="Trimming samples", logger=logger,
            quiet=False, isatty_fn=lambda: False,
        ))

    assert result == list(range(20))
    # 10% of 20 = every 2 items -> 10 milestone messages (10%, 20%, ..., 100%)
    assert len(caplog.records) == 10
    assert "10%" in caplog.records[0].message
    assert "100%" in caplog.records[-1].message
    assert "Trimming samples" in caplog.records[0].message


def test_small_total_logs_every_item(caplog):
    logger = logging.getLogger("test_progress_iterator.small_total")

    with caplog.at_level(logging.INFO, logger=logger.name):
        result = list(toolsTG.progress_iterator(
            range(3), total=3, desc="Trimming samples", logger=logger,
            quiet=False, isatty_fn=lambda: False,
        ))

    assert result == [0, 1, 2]
    assert len(caplog.records) == 3


def test_tty_without_quiet_yields_all_items():
    """
    Real rich Progress/Live objects are used here (not mocked): rich degrades gracefully on a
    non-tty output stream (which pytest's captured stdout/stderr already is) -- it just skips
    the animation, so exercising the real objects is both simpler and more honest than faking
    them, given rich.progress.Progress constructs its own internal Live regardless (a real
    dependency, not something progress_iterator can inject around).
    """
    logger = logging.getLogger("test_progress_iterator.tty")

    result = list(toolsTG.progress_iterator(
        range(5), total=5, desc="Trimming samples", logger=logger,
        quiet=False, isatty_fn=lambda: True,
    ))

    assert result == [0, 1, 2, 3, 4]


def test_tty_branch_captures_logged_messages_into_the_tail_panel():
    logger = logging.getLogger("test_progress_iterator.tty_tail")
    logger.setLevel(logging.INFO)

    gen = toolsTG.progress_iterator(
        range(2), total=2, desc="Trimming samples", logger=logger,
        quiet=False, isatty_fn=lambda: True,
    )
    next(gen)
    tail_handlers = [h for h in logger.handlers if isinstance(h, toolsTG._TailCaptureHandler)]
    assert len(tail_handlers) == 1

    logger.info("Finished: sample-0")

    assert "Finished: sample-0" in tail_handlers[0].lines
    list(gen)  # drain the rest


def test_tty_branch_removes_tail_handler_and_leaves_propagate_untouched():
    logger = logging.getLogger("test_progress_iterator.tty_cleanup")
    logger.propagate = True
    handlers_before = list(logger.handlers)

    list(toolsTG.progress_iterator(
        range(3), total=3, desc="Trimming samples", logger=logger,
        quiet=False, isatty_fn=lambda: True,
    ))

    assert logger.handlers == handlers_before
    assert logger.propagate is True


def test_tty_branch_persists_to_file_handler_while_skipping_console_handler():
    """
    The actual production mechanism (cli.py's configure_logging()): the shared 'trnagraph'
    logger carries a FileHandler (always) and, unless --quiet, a StreamHandler tagged
    `_is_console_handler`. During a live session, progress_iterator must keep messages flowing
    to the FileHandler (so the persisted .log/ file still gets a full record even while the
    console is showing the rich display instead) while skipping only the console-tagged one.
    """
    trnagraph_logger = logging.getLogger('trnagraph')
    file_records, console_records = [], []

    class _RecordingHandler(logging.Handler):
        def __init__(self, sink):
            super().__init__()
            self.sink = sink

        def emit(self, record):
            self.sink.append(record.getMessage())

    file_handler = _RecordingHandler(file_records)
    console_handler = _RecordingHandler(console_records)
    console_handler._is_console_handler = True
    trnagraph_logger.addHandler(file_handler)
    trnagraph_logger.addHandler(console_handler)
    trnagraph_logger.propagate = False  # matches configure_logging()'s own setting

    child_logger = logging.getLogger('trnagraph.modules.test_fake_tool')
    child_logger.setLevel(logging.INFO)

    gen = toolsTG.progress_iterator(
        range(1), total=1, desc="Testing", logger=child_logger,
        quiet=False, isatty_fn=lambda: True,
    )
    next(gen)
    child_logger.info("Finished: sample-0")
    list(gen)  # drain, restoring the console handler

    assert "Finished: sample-0" in file_records
    assert "Finished: sample-0" not in console_records
    # Restored after the live session ends.
    assert console_handler in trnagraph_logger.handlers
