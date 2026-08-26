"""Regression tests for the source-root detection behind `update` and the version channel.

`env_check.get_project_root()` used to walk three directories up from its own `__file__` and
return the result unconditionally. In a source checkout that lands on the repo root; from a
non-editable `pip install .` it lands on `<venv>/lib/pythonX.Y`, which is not a checkout at all.

That is not merely a wrong path. `UpdateManager` passes it as `cwd` to every git command, and git
searches *upward* from cwd -- so with the extremely common `python -m venv .venv` layout inside a
project, `update`'s `_check_clean_working_tree()` preflight inspects the **enclosing repository**,
reports it clean, and `update` then proceeds to fetch/checkout/pull against a repository that has
nothing to do with tRNAgraph. `get_version_channel()` mis-reports that repo's branch for the same
reason.

The fix makes the walk a *candidate* that must be positively identified as this project's source
tree before it is returned, and gives every caller an explicit answer for "there isn't one".
"""
import logging
import subprocess
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from trnagraph.modules import env_check, toolsUpdate


def _make_source_root(path, project_name="tRNAgraph"):
    """A directory shaped like this project's own checkout."""
    (path / "src" / "trnagraph").mkdir(parents=True)
    (path / "src" / "trnagraph" / "__init__.py").write_text('__version__ = "0.0.0"\n')
    (path / "pyproject.toml").write_text(f'[project]\nname = "{project_name}"\n')
    return path


def test_identifies_this_projects_source_checkout(tmp_path):
    assert env_check._is_trnagraph_source_root(str(_make_source_root(tmp_path))) is True


def test_rejects_an_unrelated_repository(tmp_path):
    """The actual hazard: a venv created inside someone else's project. The candidate path is a
    real directory, and git will happily answer questions about the repo above it -- but it is
    not tRNAgraph's source tree and must never be handed to `update`."""
    other = tmp_path / "someones-project"
    (other / "src").mkdir(parents=True)
    (other / "pyproject.toml").write_text('[project]\nname = "someones-project"\n')
    subprocess.run(["git", "init", "-q", str(other)], check=True, capture_output=True)

    assert env_check._is_trnagraph_source_root(str(other)) is False


def test_rejects_a_directory_with_no_pyproject(tmp_path):
    (tmp_path / "src" / "trnagraph").mkdir(parents=True)
    (tmp_path / "src" / "trnagraph" / "__init__.py").write_text("")

    assert env_check._is_trnagraph_source_root(str(tmp_path)) is False


def test_rejects_a_pyproject_without_the_package_source(tmp_path):
    (tmp_path / "pyproject.toml").write_text('[project]\nname = "tRNAgraph"\n')

    assert env_check._is_trnagraph_source_root(str(tmp_path)) is False


def test_rejects_a_malformed_pyproject(tmp_path):
    _make_source_root(tmp_path)
    (tmp_path / "pyproject.toml").write_text("this is not toml [[[")

    assert env_check._is_trnagraph_source_root(str(tmp_path)) is False


def test_get_project_root_is_none_without_a_source_checkout(tmp_path):
    with patch.object(env_check, "_is_trnagraph_source_root", return_value=False):
        assert env_check.get_project_root() is None


def test_get_requirements_path_is_none_without_a_source_checkout():
    with patch.object(env_check, "get_project_root", return_value=None):
        assert env_check.get_requirements_path() is None


def test_version_channel_reports_installed_without_a_source_checkout():
    """A wheel install has no branch, tag or hash to report. 'stable' would be an outright lie --
    a wheel built from a dirty dev branch would claim it."""
    assert env_check.get_version_channel(None) == "installed"


def test_validate_environment_is_a_no_op_without_a_source_checkout(capsys):
    """requirements.yaml is not part of the distribution, so there is nothing to validate
    against -- but the CLI must still run."""
    with patch.object(env_check, "get_project_root", return_value=None):
        env_check.validate_environment()

    assert "Validating environment" not in capsys.readouterr().out


def test_update_check_touches_no_network_without_a_source_checkout():
    with patch.object(env_check, "get_project_root", return_value=None), \
         patch.object(env_check, "_update_check_cache_path", return_value="/nonexistent/cache.json"), \
         patch("subprocess.run") as mock_run:
        env_check.check_for_updates()

    mock_run.assert_not_called()


def test_update_refuses_to_construct_without_a_source_checkout():
    """The load-bearing test. `update` must fail before any git command runs, because the first
    one it would run (`git status` in `_check_clean_working_tree`) resolves against whatever
    repository happens to sit above the install location."""
    with patch.object(toolsUpdate.env_check, "get_project_root", return_value=None), \
         patch("subprocess.run") as mock_run, \
         patch("subprocess.Popen") as mock_popen:
        with pytest.raises(ValueError, match="source checkout"):
            toolsUpdate.UpdateManager(SimpleNamespace(branch=None, tag=None, quiet=False))

    mock_run.assert_not_called()
    mock_popen.assert_not_called()


def test_update_still_constructs_from_a_source_checkout(tmp_path):
    root = str(_make_source_root(tmp_path))
    with patch.object(toolsUpdate.env_check, "get_project_root", return_value=root):
        manager = toolsUpdate.UpdateManager(SimpleNamespace(branch=None, tag=None, quiet=False))

    assert manager.project_root == root


def test_makedb_provenance_lookup_skips_git_without_a_source_checkout():
    """`makedb` stamps a git hash into its database provenance. It derived the repo root from
    its own `script_dir` walk -- a third copy of the same guess -- so it now shares the validated
    helper. Without a checkout there is no history to read and git must not be invoked."""
    from trnagraph.modules import toolsTDatabase

    builder = toolsTDatabase.tRNADatabaseBuilder.__new__(toolsTDatabase.tRNADatabaseBuilder)
    with patch.object(toolsTDatabase.env_check, "get_project_root", return_value=None), \
         patch("subprocess.check_output") as mock_check_output:
        _version_str, hash_str = builder.get_git_hash()

    assert hash_str == "Unknown"
    mock_check_output.assert_not_called()
