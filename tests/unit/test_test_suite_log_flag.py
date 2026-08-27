"""`tools test --log` is a switch, not a path.

The flag was declared as `Optional[str]` with the help text "Log output to file", but its only
consumer (`toolsTestSuite._live_box`) reads nothing except its truthiness -- it decides whether
to render the live progress panel. Nothing ever opened the path, so `--log run.txt` produced no
file, no warning and no error.

Declaring it a boolean makes the flag honest and makes the old spelling fail loudly (an
unexpected extra argument) rather than silently doing nothing. The suite's own output still goes
where it always did: a fixed `toolsTestSuite.log`, plus a timestamped `.log/` entry per
`trnagraph` invocation it makes.
"""
import inspect

from trnagraph import cli


def test_log_is_a_boolean_switch_not_a_path():
    params = inspect.signature(cli.test).parameters

    assert params["log"].annotation is bool
    assert params["log"].default.default is False


def test_log_help_describes_disabling_the_live_panel():
    """The old help promised a file that was never written."""
    params = inspect.signature(cli.test).parameters
    help_text = params["log"].default.help.lower()

    assert "live" in help_text or "panel" in help_text
    assert "log output to file" not in help_text
