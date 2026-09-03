"""A failed step makes `tools test` exit non-zero.

main() wrapped the whole pipeline in `except Exception` and only LOGGED, then returned
normally. So a run that aborted partway printed an error line, followed by "Done!", and exited
0 -- indistinguishable from success to anything reading the exit code, and easy to miss in a
long log.

That is exactly how it behaved when the build broke: the suite reported

    An error occurred during execution: Command 'trnagraph analyze build ...' returned
    non-zero exit status 1.
    Done!

and exited 0. A demo pipeline whose whole job is to tell you the pipeline works cannot report
success when it did not finish, and it cannot be used as a CI gate while it does.

The error stays a clean one-line message rather than a traceback -- that part was right, and is
the same treatment WorkspaceNotOwnedError gets.
"""
from types import SimpleNamespace

import pytest

from trnagraph.modules import toolsTestSuite


def _pipeline(failing_step=None):
    import logging
    suite = toolsTestSuite.demoPipeline.__new__(toolsTestSuite.demoPipeline)
    suite.args = SimpleNamespace(
        all=True, metadata=False, fastq=False, trna=False, genome=False, trim=False,
        makedb=False, map=False, build=False, hubonly=False, split_build=False,
        cluster=False, merge=False, graph=False, split_graph=False, skip_download=True, cleanrun=False)
    suite.logger = logging.getLogger('test-exit-code')

    def step(name):
        def run(*args, **kwargs):
            if name == failing_step:
                raise RuntimeError(f'{name} failed')
        return run

    for name in ('download_metadata', 'download_fastq', 'download_trna', 'download_genome',
                 'trim_fastq', 'create_index', 'map_reads', 'build_db', 'cluster_db',
                 'merge_db', 'graph_db', 'graph_split_db', 'generate_hub'):
        setattr(suite, name, step(name))
    return suite


def test_a_failing_step_raises_rather_than_returning():
    suite = _pipeline(failing_step='build_db')

    with pytest.raises(Exception):
        toolsTestSuite.demoPipeline.main(suite)


def test_the_failure_names_the_step_that_broke():
    suite = _pipeline(failing_step='cluster_db')

    with pytest.raises(Exception) as excinfo:
        toolsTestSuite.demoPipeline.main(suite)

    assert 'cluster_db failed' in str(excinfo.value)


def test_a_clean_run_still_returns_normally():
    toolsTestSuite.demoPipeline.main(_pipeline())
