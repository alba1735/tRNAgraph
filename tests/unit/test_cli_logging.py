"""Regression tests for cli.py's logging architecture (roadmap.md Phase 2: "Logging" +
"File-logging architecture") -- configure_logging() is the shared, application-level place
console/file handlers get attached to the 'trnagraph' logger (library/module code never
configures its own, only calls logging.getLogger(__name__) and logs; this is what everything
propagates up to), and handle_output() is the per-command context manager that always persists
a timestamped log under ./.log/, moves it to the command's real output on success, and warns
(without moving) on failure."""
import logging

import pytest

from trnagraph.cli import configure_logging, handle_output


@pytest.fixture(autouse=True)
def _reset_trnagraph_logger():
    """These mutate the shared, process-wide 'trnagraph' logger -- reset it after each test so
    tests don't leak handlers into each other."""
    yield
    logger = logging.getLogger('trnagraph')
    for handler in logger.handlers[:]:
        handler.close()
        logger.removeHandler(handler)


def test_configure_logging_always_attaches_a_file_handler(tmp_path):
    log_path = tmp_path / "run.log"
    logger = configure_logging(str(log_path), quiet=True)

    assert len(logger.handlers) == 1
    assert isinstance(logger.handlers[0], logging.FileHandler)


def test_configure_logging_not_quiet_also_attaches_a_tagged_console_handler(tmp_path):
    log_path = tmp_path / "run.log"
    logger = configure_logging(str(log_path), quiet=False)

    assert len(logger.handlers) == 2
    console_handlers = [h for h in logger.handlers if getattr(h, '_is_console_handler', False)]
    assert len(console_handlers) == 1
    assert isinstance(console_handlers[0], logging.StreamHandler)


def test_configure_logging_does_not_accumulate_handlers_across_calls(tmp_path):
    """A second CLI invocation in the same process (e.g. across pytest tests, or repeated
    programmatic use) must not pile up handlers on the shared logger name."""
    configure_logging(str(tmp_path / "first.log"), quiet=False)
    logger = configure_logging(str(tmp_path / "second.log"), quiet=False)

    assert len(logger.handlers) == 2


def test_module_logger_messages_propagate_into_the_configured_log_file(tmp_path):
    """
    The whole point of the design: a module deep in trnagraph.modules (e.g. toolsTrim.py) just
    calls logging.getLogger(__name__) and logs -- it never touches handlers itself. Confirms
    that actually works end-to-end via propagation up to the 'trnagraph' logger configured here.
    """
    log_path = tmp_path / "run.log"
    configure_logging(str(log_path), quiet=True)

    module_logger = logging.getLogger('trnagraph.modules.some_future_module')
    module_logger.info('a message from deep inside a converted module')

    log_contents = log_path.read_text()
    assert 'a message from deep inside a converted module' in log_contents


def test_handle_output_persists_a_timestamped_log_under_dot_log(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    with handle_output(quiet=True, tool="testtool"):
        pass

    log_files = list((tmp_path / ".log").glob("*_testtool.log"))
    assert len(log_files) == 1


def test_handle_output_moves_the_log_to_destination_on_success(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    destination = tmp_path / "results"

    with handle_output(quiet=True, tool="testtool", destination=str(destination)):
        pass

    assert not list((tmp_path / ".log").glob("*_testtool.log"))
    assert len(list(destination.glob("*_testtool.log"))) == 1


def test_handle_output_leaves_log_in_dot_log_and_warns_on_failure(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    destination = tmp_path / "results"

    with pytest.raises(ValueError):
        with handle_output(quiet=True, tool="testtool", destination=str(destination)):
            raise ValueError("boom")

    # Log stays in .log/ -- the destination may not have been correctly produced by a run that crashed.
    assert len(list((tmp_path / ".log").glob("*_testtool.log"))) == 1
    assert not destination.exists()
    stderr = capsys.readouterr().err
    assert "testtool" in stderr
    assert ".log" in stderr


def test_handle_output_logs_the_exception_to_the_log_file_on_failure(tmp_path, monkeypatch):
    """A real server failure (trnagraph update's _check_clean_working_tree raising ValueError)
    showed the log file containing only the startup banner lines -- the actual exception/
    traceback never reached it, only the terminal (rendered by Typer/Click's own top-level
    handler, outside this function, by which point handle_output's `finally` had already closed
    the FileHandler). The exception must be logged -- with its traceback -- before that happens."""
    monkeypatch.chdir(tmp_path)

    with pytest.raises(ValueError):
        with handle_output(quiet=True, tool="testtool"):
            raise ValueError("boom: something specific went wrong")

    log_files = list((tmp_path / ".log").glob("*_testtool.log"))
    assert len(log_files) == 1
    content = log_files[0].read_text()
    assert "ValueError" in content
    assert "boom: something specific went wrong" in content
    assert "Traceback" in content


def test_handle_output_still_writes_file_log_under_quiet(tmp_path, monkeypatch):
    """--quiet only suppresses the console -- the persisted file log is unconditional."""
    monkeypatch.chdir(tmp_path)

    with handle_output(quiet=True, tool="testtool"):
        print("this should reach the log file even though quiet")

    log_files = list((tmp_path / ".log").glob("*_testtool.log"))
    assert len(log_files) == 1
    assert "this should reach the log file even though quiet" in log_files[0].read_text()


def test_handle_output_filename_includes_name_suffix(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    with handle_output(quiet=True, tool="addsplit", name_suffix="vibrChol1"):
        pass

    log_files = list((tmp_path / ".log").glob("*_addsplit_vibrChol1.log"))
    assert len(log_files) == 1
