"""Regression test for env_check.py's binary version checking. A server report of `umi_tools`
failing with "Command 'umi_tools' found but version could not be determined" (while the same
version worked fine locally) pointed at a real gap: get_binary_version() had a fixed 5s
subprocess timeout and collapsed every failure mode (timeout, no regex match, any other
exception) into the same generic, undiagnosable message. A Python-based CLI tool has to boot a
full interpreter and import numpy/pandas/pysam just to print --version, which can push past a
short timeout under heavy concurrent load (e.g. a pipeline already running on the same machine) --
that failure mode looks identical to "the output format doesn't match" unless reported
separately."""
import subprocess
from unittest.mock import patch

from trnagraph.modules.env_check import get_binary_version, check_binary_package


def test_returns_the_version_on_a_normal_match():
    with patch('shutil.which', return_value='/usr/bin/fake'), \
         patch('subprocess.run', return_value=subprocess.CompletedProcess([], 0, stdout='tool version: 1.2.3\n', stderr='')):
        version, reason = get_binary_version('fake', '--version', r': ([\d\.]+)')

    assert version == '1.2.3'
    assert reason is None


def test_reports_timeout_distinctly_from_a_generic_failure():
    with patch('shutil.which', return_value='/usr/bin/fake'), \
         patch('subprocess.run', side_effect=subprocess.TimeoutExpired(cmd='fake', timeout=20)):
        version, reason = get_binary_version('fake', '--version', r': ([\d\.]+)', timeout=20)

    assert version is None
    assert 'timed out' in reason
    assert '20' in reason


def test_reports_no_regex_match_distinctly():
    with patch('shutil.which', return_value='/usr/bin/fake'), \
         patch('subprocess.run', return_value=subprocess.CompletedProcess([], 0, stdout='unexpected output\n', stderr='')):
        version, reason = get_binary_version('fake', '--version', r': ([\d\.]+)')

    assert version is None
    assert reason is not None
    assert 'timed out' not in reason


def test_check_binary_package_surfaces_the_failure_reason():
    with patch('shutil.which', return_value='/usr/bin/umi_tools'), \
         patch('subprocess.run', side_effect=subprocess.TimeoutExpired(cmd='umi_tools', timeout=20)):
        ok, message = check_binary_package('umi_tools', ('>=', '1.1.6'))

    assert ok is False
    assert 'timed out' in message


def test_default_timeout_is_generous_enough_for_a_python_based_tool():
    """20s, not the old 5s -- a Python CLI tool booting an interpreter + numpy/pandas/pysam
    imports under load needs real headroom, not a timeout tuned for a fast C binary."""
    import inspect
    assert inspect.signature(get_binary_version).parameters['timeout'].default >= 15
