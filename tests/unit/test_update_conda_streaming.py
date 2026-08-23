"""Regression test for roadmap.md's update-tool "hangs on Syncing conda environment" complaint.
UpdateManager._run() buffers a subprocess's entire stdout/stderr via subprocess.run(...,
capture_output=True) and only logs it once the process exits -- for a slow conda/mamba solve
(a well-known, not-really-fixable conda characteristic), that means literally nothing is visible
for however long it takes, indistinguishable from an actual hang. Fix: stream the conda/mamba
env update step's output live instead."""
import logging
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from trnagraph.modules.toolsUpdate import UpdateManager


def _make_manager():
    manager = UpdateManager.__new__(UpdateManager)
    manager.args = SimpleNamespace(branch=None, tag=None, quiet=False)
    manager.logger = logging.getLogger('trnagraph.modules.toolsUpdate')
    manager.project_root = '/fake/project/root'
    return manager


def _fake_popen(lines, returncode=0):
    proc = MagicMock()
    proc.stdout = iter(lines)
    proc.wait.return_value = returncode
    proc.returncode = returncode
    return proc


def test_run_streaming_logs_each_line_as_it_arrives(caplog):
    manager = _make_manager()
    lines = ['Collecting package metadata...\n', 'Solving environment...\n', 'Done\n']

    with patch('subprocess.Popen', return_value=_fake_popen(lines)) as mock_popen, \
         caplog.at_level(logging.INFO, logger='trnagraph.modules.toolsUpdate'):
        manager._run_streaming(['mamba', 'env', 'update', '-f', 'requirements.yaml', '--prune'])

    messages = [r.message for r in caplog.records]
    assert 'Collecting package metadata...' in messages
    assert 'Solving environment...' in messages
    assert 'Done' in messages
    mock_popen.assert_called_once()


def test_run_streaming_raises_on_nonzero_returncode():
    manager = _make_manager()

    with patch('subprocess.Popen', return_value=_fake_popen(['error output\n'], returncode=1)):
        with pytest.raises(ValueError, match='error output'):
            manager._run_streaming(['mamba', 'env', 'update'])


def test_conda_env_update_uses_streaming_not_buffered_run():
    manager = _make_manager()
    manager._run_streaming = MagicMock()
    manager._run = MagicMock()

    with patch('shutil.which', return_value='/usr/bin/mamba'), \
         patch('os.path.exists', return_value=True), \
         patch('trnagraph.modules.toolsUpdate.env_check.get_requirements_path', return_value='/fake/requirements.yaml'):
        manager._conda_env_update('mamba')

    manager._run_streaming.assert_called_once()
    manager._run.assert_not_called()
