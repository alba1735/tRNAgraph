"""Tests that `preprocess map`'s deduplication flags reach the mapping stage.

Deduplication is opt-in and must stay that way: running it by default would silently change
every existing user's counts, and on a dataset without UMIs it would be destructive. These pin
that the flags exist, that they carry through to the stage that acts on them, and -- most
importantly -- that nothing deduplicates unless asked.
"""
import inspect
import os
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from trnagraph import cli
from trnagraph.modules import toolsMap


def test_map_exposes_the_dedup_flags():
    params = inspect.signature(cli.map_cmd).parameters

    assert "dedup" in params
    assert "keep_prededup" in params
    assert "dedup_method" in params


def test_dedup_is_off_by_default():
    """Opt-in, because deduplication changes counts and cannot be undone without a remap."""
    params = inspect.signature(cli.map_cmd).parameters

    assert params["dedup"].default.default is False
    assert params["keep_prededup"].default.default is False
    assert params["dedup_method"].default.default == "directional"


def _make_mapper(tmp_path, **overrides):
    args = SimpleNamespace(
        mode="map", output=str(tmp_path / "exp"), database="db",
        input=str(tmp_path / "meta.tsv"), force_remap=False, minnontrnasize=20,
        local=False, threads=1, skipcheck=False, bamdir=str(tmp_path / "bam"),
        quiet=True, dedup=False, keep_prededup=False, dedup_method="directional",
    )
    for key, value in overrides.items():
        setattr(args, key, value)
    with patch.object(toolsMap, "trnadatabase"), patch.object(toolsMap, "expdatabase"):
        mapper = toolsMap.MapSamples(args)
    # expdatabase is mocked, so give the two paths the dedup phase actually writes to real
    # values rather than Mock objects.
    mapper.expinfo.resultsdir = str(tmp_path / "results")
    mapper.expinfo.dedupinfo = str(tmp_path / "results" / "exp-dedupinfo.txt")
    return mapper


def test_mapper_reads_the_dedup_settings_off_args(tmp_path):
    mapper = _make_mapper(tmp_path, dedup=True, keep_prededup=True, dedup_method="unique")

    assert mapper.dedup is True
    assert mapper.keep_prededup is True
    assert mapper.dedup_method == "unique"


def test_mapper_defaults_to_no_dedup_when_args_predate_the_flag(tmp_path):
    """`analyze addsplit` and the test suite build args namespaces by hand."""
    args = SimpleNamespace(
        mode="map", output=str(tmp_path / "exp"), database="db",
        input=str(tmp_path / "meta.tsv"), force_remap=False, minnontrnasize=20,
        local=False, threads=1, skipcheck=False, bamdir=str(tmp_path / "bam"), quiet=True,
    )
    with patch.object(toolsMap, "trnadatabase"), patch.object(toolsMap, "expdatabase"):
        mapper = toolsMap.MapSamples(args)

    assert mapper.dedup is False


def _sample_names(tmp_path, names):
    """Writes a minimal metadata file and the bams `map` would have produced."""
    bamdir = tmp_path / "bam"
    bamdir.mkdir(parents=True, exist_ok=True)
    rows = ["fastq\tsample\tgroup"]
    for name in names:
        (bamdir / f"{name}.bam").write_bytes(b"")
        rows.append(f"{name}.fastq.gz\t{name}\tGroupA")
    (tmp_path / "meta.tsv").write_text("\n".join(rows) + "\n")


def test_dedup_runs_once_per_sample_after_mapping(tmp_path):
    _sample_names(tmp_path, ["SampleA", "SampleB"])
    mapper = _make_mapper(tmp_path, dedup=True, dedup_method="unique", keep_prededup=True)

    with patch.object(mapper, "mapsamples"), \
         patch.object(toolsMap.toolsDedup, "dedup_sample") as mock_dedup, \
         patch.object(toolsMap.toolsDedup, "detect_umi_separator", return_value="_"):
        mapper.main()

    assert mock_dedup.call_count == 2
    for call in mock_dedup.call_args_list:
        assert call.kwargs["method"] == "unique"
        assert call.kwargs["keep_prededup"] is True
    deduped = sorted(os.path.basename(call.args[0]) for call in mock_dedup.call_args_list)
    assert deduped == ["SampleA.bam", "SampleB.bam"]


def test_nothing_is_deduplicated_without_the_flag(tmp_path):
    _sample_names(tmp_path, ["SampleA"])
    mapper = _make_mapper(tmp_path, dedup=False)

    with patch.object(mapper, "mapsamples"), \
         patch.object(toolsMap.toolsDedup, "dedup_sample") as mock_dedup:
        mapper.main()

    mock_dedup.assert_not_called()


def test_dedup_run_records_its_provenance(tmp_path):
    """What was deduplicated, how, and against which read-name shape.

    `map` writes no runinfo of its own (that is an `analyze build` artifact), so deduplication
    -- which silently changes every downstream count -- would otherwise leave no record that it
    happened at all. Read counts deliberately live in umi_tools' own per-sample log rather than
    here: recovering them would mean a second full pass over every bam.
    """
    _sample_names(tmp_path, ["SampleA", "SampleB"])
    mapper = _make_mapper(tmp_path, dedup=True, dedup_method="unique")

    with patch.object(mapper, "mapsamples"), \
         patch.object(toolsMap.toolsDedup, "dedup_sample"), \
         patch.object(toolsMap.toolsDedup, "detect_umi_separator", return_value="_"), \
         patch.object(toolsMap.toolsDedup, "umi_tools_version", return_value="1.1.6"):
        mapper.main()

    written = open(mapper.expinfo.dedupinfo).read()
    assert "unique" in written
    assert "1.1.6" in written
    assert "SampleA" in written and "SampleB" in written


def test_dedup_concurrency_follows_threads_and_memory(tmp_path):
    """umi_tools is single-threaded, so samples are the only unit of parallelism -- but memory,
    not cores, is what actually limits it: roughly 7GB per concurrent job on a human-scale bam.

    Concurrency is therefore the smallest of half the thread count, what system memory affords,
    and the number of samples. Serially this phase took 2.4h on a nine-sample human dataset
    against a ~32min slowest sample, so the win is most of the wall clock.
    """
    _sample_names(tmp_path, ["S1", "S2", "S3", "S4", "S5", "S6"])
    mapper = _make_mapper(tmp_path, dedup=True, threads=8)

    seen = {}

    def fake_pool(processes=None, **kwargs):
        seen["processes"] = processes
        raise AssertionError("pool constructed")

    with patch.object(mapper, "mapsamples"), \
         patch.object(toolsMap.toolsDedup, "detect_umi_separator", return_value="_"), \
         patch.object(toolsMap, "_memory_job_limit", return_value=99), \
         patch.object(toolsMap, "Pool", side_effect=fake_pool):
        with pytest.raises(AssertionError):
            mapper.main()

    assert seen["processes"] == 4, "8 threads should allow 4 concurrent dedup jobs"


def test_dedup_concurrency_is_capped_by_available_memory(tmp_path):
    """A machine with many cores but little RAM must not launch a job per core.

    Without this bound, 40 threads would start 20 jobs at ~7GB each -- 140GB -- and be killed
    partway through, after some samples had already been deduplicated.
    """
    _sample_names(tmp_path, ["S1", "S2", "S3", "S4", "S5", "S6"])
    mapper = _make_mapper(tmp_path, dedup=True, threads=40)

    seen = {}

    def fake_pool(processes=None, **kwargs):
        seen["processes"] = processes
        raise AssertionError("pool constructed")

    with patch.object(mapper, "mapsamples"), \
         patch.object(toolsMap.toolsDedup, "detect_umi_separator", return_value="_"), \
         patch.object(toolsMap, "_memory_job_limit", return_value=2), \
         patch.object(toolsMap, "Pool", side_effect=fake_pool):
        with pytest.raises(AssertionError):
            mapper.main()

    assert seen["processes"] == 2, "memory, not cores, should bind here"


def test_dedup_concurrency_never_exceeds_the_sample_count(tmp_path):
    _sample_names(tmp_path, ["S1", "S2"])
    mapper = _make_mapper(tmp_path, dedup=True, threads=64)

    seen = {}

    def fake_pool(processes=None, **kwargs):
        seen["processes"] = processes
        raise AssertionError("pool constructed")

    with patch.object(mapper, "mapsamples"), \
         patch.object(toolsMap.toolsDedup, "detect_umi_separator", return_value="_"), \
         patch.object(toolsMap, "_memory_job_limit", return_value=99), \
         patch.object(toolsMap, "Pool", side_effect=fake_pool):
        with pytest.raises(AssertionError):
            mapper.main()

    assert seen["processes"] == 2


def test_memory_job_limit_scales_with_the_largest_bam(tmp_path):
    """The estimate is derived from the bams, since real usage depends on read count and UMI
    length rather than being a fixed figure.

    Sized off the largest rather than the mean: the pool hands out jobs in no guaranteed order,
    so the biggest ones can be resident together and an average would under-estimate exactly
    when it matters.
    """
    small = tmp_path / "small.bam"
    large = tmp_path / "large.bam"
    small.write_bytes(b"x" * 1024)
    large.write_bytes(b"x" * 4096)

    only_small = toolsMap._memory_job_limit([str(small)], multiple=1)
    with_large = toolsMap._memory_job_limit([str(small), str(large)], multiple=1)

    # Four times the bytes per job means at most a quarter as many concurrent jobs.
    assert only_small == with_large * 4


def test_memory_job_limit_is_absent_rather_than_wrong_when_undeterminable(tmp_path):
    assert toolsMap._memory_job_limit([str(tmp_path / "missing.bam")]) is None
    assert toolsMap._memory_job_limit([]) is None


def test_dedup_phase_sizes_its_pool_without_patching(tmp_path):
    """Exercises the real sizing path end to end.

    The other concurrency tests patch _memory_job_limit, which meant they still passed when the
    pool was sized from a variable that had not been assigned yet. This one does not patch it, so
    the wiring itself is covered.

    threads=2 forces the serial branch so dedup_sample can be observed in this process -- a real
    pool would run it in children where the patch does not apply. The sizing block runs before
    that branch either way, which is the part being covered.
    """
    _sample_names(tmp_path, ["S1", "S2"])
    mapper = _make_mapper(tmp_path, dedup=True, threads=2)

    with patch.object(mapper, "mapsamples"), \
         patch.object(toolsMap.toolsDedup, "detect_umi_separator", return_value="_"), \
         patch.object(toolsMap.toolsDedup, "dedup_sample") as mock_dedup:
        mapper.main()

    assert mock_dedup.call_count == 2
