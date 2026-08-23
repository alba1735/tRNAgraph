"""Regression test for roadmap.md's "screen terminal bug". Root cause turned out not to be a
tRNAgraph bug at all: `conda run` defaults to fully capturing/buffering stdout+stderr until the
wrapped process exits (its own --help calls --no-capture-output what "enables interactive mode"),
which run_tg.bash's `conda run -n $ENV_NAME trnagraph ...` invocations never passed -- making
every step look silent for its entire duration, in any terminal, screen included. Since the same
footgun could recur for anyone wrapping trnagraph in an output-capturing command, this adds a
one-time startup hint: if stdout/stderr are non-interactive despite $STY/$TMUX suggesting a
human-attended screen/tmux session, say so."""
from trnagraph.modules.env_check import warn_if_output_capture_suspected


def test_no_warning_when_genuinely_interactive():
    msg = warn_if_output_capture_suspected(isatty_fn=lambda: True, environ={'TMUX': '/tmp/tmux-1000/default,123,0'})
    assert msg is None


def test_no_warning_when_non_interactive_and_no_multiplexer_env():
    """A genuine background/cron/CI invocation (no $STY/$TMUX) shouldn't get this hint -- it's
    correctly non-interactive, nothing suspicious about it."""
    msg = warn_if_output_capture_suspected(isatty_fn=lambda: False, environ={})
    assert msg is None


def test_warns_when_non_interactive_inside_tmux():
    msg = warn_if_output_capture_suspected(isatty_fn=lambda: False, environ={'TMUX': '/tmp/tmux-1000/default,123,0'})
    assert msg is not None
    assert '--no-capture-output' in msg


def test_warns_when_non_interactive_inside_screen():
    msg = warn_if_output_capture_suspected(isatty_fn=lambda: False, environ={'STY': '12345.pts-0.host'})
    assert msg is not None
    assert '--no-capture-output' in msg
