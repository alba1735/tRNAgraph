"""A full `tools test` run builds once, with the read-length split.

It used to build twice -- once plain, once split -- into one directory, which left the plain
run's flat `results/` beside the split run's `results/complete/` holding byte-identical copies.
Separating them into two directories fixed the collision but kept the duplication: every file in
`vibrChol1_nosplit/results/` matched `vibrChol1/results/complete/` except the runinfo line, and
nothing downstream read the plain tree.

So the plain build is gone from a full run. `results/complete/` IS the non-split output, and the
flat-vs-nested path rule it incidentally exercised is pinned directly in
test_variant_dir_names.py instead, where it costs a millisecond rather than a whole build.
"""
import contextlib
import logging
from types import SimpleNamespace
from unittest.mock import patch

from trnagraph.modules import toolsTestSuite


_STEPS = ["download_metadata", "download_fastq", "download_trna", "download_genome",
          "trim_fastq", "create_index", "map_reads", "build_db", "cluster_db",
          "graph_db", "graph_split_db", "_cleanup_workspace"]


def _pipeline(**flags):
    defaults = dict(all=False, metadata=False, fastq=False, trna=False, genome=False,
                    trim=False, makedb=False, map=False, build=False, hubonly=False,
                    split_build=False, cluster=False, merge=False, graph=False,
                    split_graph=False, skip_download=True, cleanrun=False)
    defaults.update(flags)
    suite = toolsTestSuite.demoPipeline.__new__(toolsTestSuite.demoPipeline)
    suite.args = SimpleNamespace(**defaults)
    suite.logger = logging.getLogger('test_testsuite_build_sequence')
    return suite


def _steps_run(suite):
    called = []
    with contextlib.ExitStack() as stack:
        for name in _STEPS:
            stack.enter_context(patch.object(toolsTestSuite.demoPipeline, name,
                                             lambda self, _n=name: called.append(_n)))
        suite.main()
    return called


def test_a_full_run_builds_once():
    """Twice produced two trees with the same numbers in them."""
    assert _steps_run(_pipeline(all=True)).count('build_db') == 1


def test_the_full_run_builds_with_the_split():
    """cluster and graph read the object this writes, so it has to carry every variant."""
    suite = _pipeline(all=True)
    _steps_run(suite)

    assert suite.args.split_build is True


def test_build_alone_still_builds_without_the_split():
    """The flag keeps its meaning for anyone driving one step by hand."""
    suite = _pipeline(build=True)
    called = _steps_run(suite)

    assert called.count('build_db') == 1
    assert suite.args.split_build is False
