"""Regression tests for cli.py's configure_logging() (roadmap.md Phase 2: "Trimmer log/quiet"
+ "Logging") -- the shared, application-level place --log/--quiet handlers get attached to the
'trnagraph' logger. Per Python's logging documentation, library/module code (trnagraph.modules.*)
should never configure its own handlers, only call logging.getLogger(__name__) and log; this is
what everything propagates up to. Centralizing it here (instead of each module reimplementing
its own handler setup) is what lets the roadmap's incremental print()->logging conversion of
other files get --log/--quiet support for free."""
import logging

import pytest

from trnagraph.cli import configure_logging


@pytest.fixture(autouse=True)
def _reset_trnagraph_logger():
    """Configure_logging mutates the shared, process-wide 'trnagraph' logger -- reset it after
    each test so tests don't leak handlers into each other."""
    yield
    logger = logging.getLogger('trnagraph')
    for handler in logger.handlers[:]:
        handler.close()
        logger.removeHandler(handler)


def test_configure_logging_with_log_file_attaches_only_a_file_handler(tmp_path):
    log_path = tmp_path / "run.log"
    logger = configure_logging(log_file=str(log_path))

    assert len(logger.handlers) == 1
    assert isinstance(logger.handlers[0], logging.FileHandler)


def test_configure_logging_quiet_attaches_no_handlers():
    logger = configure_logging(log_file=None, quiet=True)

    assert logger.handlers == []


def test_configure_logging_default_attaches_a_stream_handler():
    logger = configure_logging(log_file=None, quiet=False)

    assert len(logger.handlers) == 1
    assert isinstance(logger.handlers[0], logging.StreamHandler)


def test_configure_logging_does_not_accumulate_handlers_across_calls(tmp_path):
    """A second CLI invocation in the same process (e.g. across pytest tests, or repeated
    programmatic use) must not pile up handlers on the shared logger name."""
    configure_logging(log_file=str(tmp_path / "first.log"))
    logger = configure_logging(log_file=str(tmp_path / "second.log"))

    assert len(logger.handlers) == 1


def test_module_logger_messages_propagate_into_the_configured_log_file(tmp_path):
    """
    The whole point of the design: a module deep in trnagraph.modules (e.g. toolsTrim.py) just
    calls logging.getLogger(__name__) and logs -- it never touches handlers itself. Confirms
    that actually works end-to-end via propagation up to the 'trnagraph' logger configured here.
    """
    log_path = tmp_path / "run.log"
    configure_logging(log_file=str(log_path))

    module_logger = logging.getLogger('trnagraph.modules.some_future_module')
    module_logger.info('a message from deep inside a converted module')

    log_contents = log_path.read_text()
    assert 'a message from deep inside a converted module' in log_contents
