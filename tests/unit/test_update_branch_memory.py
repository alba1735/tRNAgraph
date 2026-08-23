"""Regression test for roadmap.md's update-tool branch-memory item. `trnagraph update` (with no
--branch) always defaulted to 'main' unconditionally, which meant a checkout deliberately moved to
'dev' via `--branch dev` would silently snap back to 'main' on the next plain `update` call --
main being deliberately behind dev during this project's stabilization is exactly what triggered
the real server failure in testing_docs/error_log_update_server.txt. Fix: default to whatever
branch is currently checked out (git IS the memory) instead of a hardcoded 'main'; --branch still
actively switches branches as before, unchanged."""
import logging
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from trnagraph.modules.toolsUpdate import UpdateManager


def _make_manager(branch=None, current_branch='dev'):
    manager = UpdateManager.__new__(UpdateManager)
    manager.args = SimpleNamespace(branch=branch, tag=None, quiet=False)
    manager.logger = logging.getLogger('trnagraph.modules.toolsUpdate')
    manager.project_root = '/fake/project/root'
    manager._run = MagicMock(return_value=MagicMock(stdout=current_branch, returncode=0))
    return manager


def test_defaults_to_the_currently_checked_out_branch_when_no_branch_flag_given():
    manager = _make_manager(branch=None, current_branch='dev')

    result = manager._resolve_target_branch()

    assert result == 'dev'


def test_explicit_branch_flag_is_used_verbatim_without_checking_git():
    manager = _make_manager(branch='some-feature-branch')

    result = manager._resolve_target_branch()

    assert result == 'some-feature-branch'
    manager._run.assert_not_called()


def test_raises_a_clear_error_on_detached_head():
    manager = _make_manager(branch=None, current_branch='HEAD')

    with pytest.raises(ValueError, match='detached HEAD'):
        manager._resolve_target_branch()
