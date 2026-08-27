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
