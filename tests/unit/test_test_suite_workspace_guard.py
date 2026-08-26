"""`tools test --all` wipes its workspace; this guards what it is allowed to wipe.

The cleanup is an unconditional `find . -maxdepth 1 -mindepth 1 -not -name
"toolsTestSuite.log" -exec rm -rf {} +` run from inside the resolved workspace. Nothing
checked *what* was in there, so pointing `-d` at a directory holding unrelated work
destroyed it without challenge.

The suite's own top-level footprint is a small fixed set, so the check is an allowlist:
anything else present means this is not a directory the suite owns, and the wipe refuses
rather than proceeding. There is deliberately no --force override -- a flag that skips a
recursive delete check gets reused reflexively once it is in someone's shell history.
"""
import pytest

from trnagraph.modules.toolsTestSuite import (
    WORKSPACE_ENTRIES,
    WorkspaceNotOwnedError,
    unexpected_workspace_entries,
    verify_workspace_is_ours,
)


def _populate(path, names):
    for name in names:
        target = path / name
        if name.endswith(".log"):
            target.write_text("log\n")
        else:
            target.mkdir()
    return path


def test_a_directory_holding_only_suite_output_is_accepted(tmp_path):
    _populate(tmp_path, sorted(WORKSPACE_ENTRIES))

    assert unexpected_workspace_entries(str(tmp_path)) == []
    verify_workspace_is_ours(str(tmp_path))  # must not raise


def test_an_empty_or_missing_directory_is_accepted(tmp_path):
    """A first run creates the workspace; there is nothing to protect yet."""
    assert unexpected_workspace_entries(str(tmp_path)) == []
    assert unexpected_workspace_entries(str(tmp_path / "not-created-yet")) == []


def test_a_partial_workspace_is_accepted(tmp_path):
    """An interrupted run leaves only some of the entries; that is still ours."""
    _populate(tmp_path, ["config", "raw"])

    assert unexpected_workspace_entries(str(tmp_path)) == []


def test_unrelated_work_is_reported(tmp_path):
    _populate(tmp_path, ["config", "raw", "my_thesis", "notes.log"])

    assert unexpected_workspace_entries(str(tmp_path)) == ["my_thesis", "notes.log"]


def test_the_wipe_refuses_and_names_what_it_found(tmp_path):
    _populate(tmp_path, ["config", "important_data"])

    with pytest.raises(WorkspaceNotOwnedError) as excinfo:
        verify_workspace_is_ours(str(tmp_path))

    message = str(excinfo.value)
    assert "important_data" in message, "the user must be told what blocked the wipe"
    assert str(tmp_path) in message, "and which directory it was looking at"
    assert "config" not in message.replace(str(tmp_path), ""), "only the unexpected entries"


def test_hidden_entries_are_checked_too(tmp_path):
    """`find -mindepth 1` deletes dotfiles as well, so the guard must see them."""
    _populate(tmp_path, [".log"])
    (tmp_path / ".ssh").mkdir()

    assert unexpected_workspace_entries(str(tmp_path)) == [".ssh"]


def test_results_are_sorted_for_a_stable_message(tmp_path):
    _populate(tmp_path, ["zeta", "alpha", "mid"])

    assert unexpected_workspace_entries(str(tmp_path)) == ["alpha", "mid", "zeta"]


def test_allowlist_covers_what_a_real_run_leaves_behind():
    """Measured against a completed `tools test --all` workspace."""
    observed = {".log", "config", "processed", "raw", "references",
                "toolsTestSuite.log", "vibrChol1"}

    assert observed <= WORKSPACE_ENTRIES, f"missing: {observed - WORKSPACE_ENTRIES}"


def test_no_force_style_override_exists():
    """A flag that skips the check would be reused reflexively; there must not be one."""
    import inspect
    from trnagraph.modules import toolsTestSuite

    source = inspect.getsource(toolsTestSuite.verify_workspace_is_ours)
    assert "force" not in source.lower()
