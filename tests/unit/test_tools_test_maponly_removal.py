"""Regression tests for roadmap.md Phase 2's `-maponly` removal from `tools test`. Confirmed a
pure convenience shortcut for stopping the demo pipeline after mapping, with no unique
test-harness value beyond what the already-separate --trim/--makedb/--map flags provide: passing
those three specific flags together already runs exactly download->trim->makedb->map and stops,
with no --maponly needed. These tests exercise demoPipeline.main()'s step-gating logic directly
(mocking out the real step methods, which do real network/subprocess work) rather than the CLI
parsing layer."""
from types import SimpleNamespace
from unittest.mock import patch

from trnagraph.modules.toolsTestSuite import demoPipeline

_STEP_METHODS = [
    "download_metadata", "download_fastq", "download_trna", "download_genome",
    "trim_fastq", "create_index", "map_reads", "build_db", "cluster_db",
    "graph_db", "graph_split_db", "_cleanup_workspace",
]


def _make_pipeline(**flag_overrides):
    pipeline = demoPipeline.__new__(demoPipeline)
    defaults = dict(
        metadata=False, fastq=False, trna=False, genome=False, trim=False, makedb=False,
        map=False, hubonly=False, build=False, split_build=False, cluster=False, merge=False,
        graph=False, split_graph=False, all=False, skip_download=True, cleanrun=False,
    )
    defaults.update(flag_overrides)
    pipeline.args = SimpleNamespace(**defaults)
    import logging
    pipeline.logger = logging.getLogger("test_tools_test_maponly_removal")
    return pipeline


def _run_and_record_called_steps(pipeline):
    called = []
    patches = [patch.object(demoPipeline, name, lambda self, n=name: called.append(n)) for name in _STEP_METHODS]
    with patches[0], patches[1], patches[2], patches[3], patches[4], patches[5], \
         patches[6], patches[7], patches[8], patches[9], patches[10], patches[11]:
        pipeline.main()
    return called


def test_maponly_is_not_a_recognized_attribute_on_the_args_namespace():
    """The flag itself is gone -- demoPipeline.main() must not reference self.args.maponly at
    all. main() catches and logs exceptions internally rather than raising, so the real
    assertion here is that steps actually ran, not just that nothing escaped."""
    pipeline = _make_pipeline()
    assert not hasattr(pipeline.args, "maponly")
    called = _run_and_record_called_steps(pipeline)
    assert called, "main() must not silently fail via a caught AttributeError on .maponly"


def test_no_flags_runs_every_step_including_build_cluster_and_graph():
    """run_all (no specific flags given) must still run the full pipeline through graph/split
    steps -- this used to be gated by `and not self.args.maponly`, which only ever mattered when
    --maponly was passed; with the flag gone, run_all should behave as if maponly was never set."""
    pipeline = _make_pipeline()
    called = _run_and_record_called_steps(pipeline)

    for step in ["trim_fastq", "create_index", "map_reads", "build_db", "cluster_db", "graph_db", "graph_split_db"]:
        assert step in called, f"{step} should run under run_all"
    assert called.count("build_db") == 2  # once plain, once for split_build


def test_map_trim_makedb_flags_alone_stop_after_mapping():
    """The documented replacement for --maponly: passing --trim --makedb --map together already
    runs exactly download->trim->makedb->map and stops, with no dedicated flag needed."""
    pipeline = _make_pipeline(trim=True, makedb=True, map=True)
    called = _run_and_record_called_steps(pipeline)

    assert called == ["trim_fastq", "create_index", "map_reads"]
