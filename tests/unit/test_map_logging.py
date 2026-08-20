"""Regression tests for toolsMap.py's progress reporting (roadmap.md Phase 2: "tqdm" stage 2) --
`MapSamples.mapsamples()` previously drove its multiprocessing pool with a bare
`pool.imap_unordered()` loop and no progress feedback at all for long-running per-sample mapping.
It now reuses the shared `toolsTG.progress_iterator()` helper (same pattern as toolsTrim.py's
`process()`), verified here via pytest's caplog fixture against the milestone-logging path
(non-tty, the default under pytest's captured streams) rather than mocking bowtie2/samtools
directly."""
import logging
from types import SimpleNamespace
from unittest.mock import patch

from trnagraph.modules.toolsMap import MapSamples, MapInfo, TRNAMapInfo


def _make_manifest(tmp_path, samplenames):
    fastqs = []
    for name in samplenames:
        fq = tmp_path / f"{name}.fastq"
        fq.write_text("@read\nACGT\n+\nIIII\n")
        fastqs.append(fq)
    manifest = tmp_path / "manifest.txt"
    lines = ["fastq\tsample\tgroup"]
    lines += [f"{fq}\t{name}\t{name}" for fq, name in zip(fastqs, samplenames)]
    manifest.write_text("\n".join(lines) + "\n")
    return str(manifest)


def _make_args(tmp_path, samplenames, quiet=False):
    return SimpleNamespace(
        mode='map', output=str(tmp_path / "exp"), database=str(tmp_path / "db"),
        input=_make_manifest(tmp_path, samplenames), lazy=False, minnontrnasize=20,
        local=False, threads=1, skipcheck=True, bamdir=str(tmp_path / "bam"), quiet=quiet,
    )


def _fake_map_info(samplename, *_args, **_kwargs):
    trnamapinfo = TRNAMapInfo(0, 0, 0, 2, 1, 0)
    return MapInfo(1, 0, 1, 2, "bowtie output", samplename, bowtiecommand="bowtie2 ...", trnamapinfo=trnamapinfo)


def _fake_map_sample(samplename, fastqfile, bamfile, expname, logfile=None):
    return _fake_map_info(samplename)


def test_mapsamples_logs_progress_milestones(tmp_path, caplog):
    samplenames = [f"sample{i}" for i in range(3)]
    mapper = MapSamples(_make_args(tmp_path, samplenames))

    with patch("trnagraph.modules.toolsMap.MapReads.map_sample", side_effect=_fake_map_sample), \
         caplog.at_level(logging.INFO, logger="trnagraph.modules.toolsMap"):
        mapper.main()

    messages = "\n".join(r.message for r in caplog.records)
    assert "Mapping samples" in messages
    assert "100%" in messages


def test_mapsamples_quiet_still_emits_milestones_for_file_persistence(tmp_path, caplog):
    samplenames = [f"sample{i}" for i in range(3)]
    mapper = MapSamples(_make_args(tmp_path, samplenames, quiet=True))

    with patch("trnagraph.modules.toolsMap.MapReads.map_sample", side_effect=_fake_map_sample), \
         caplog.at_level(logging.INFO, logger="trnagraph.modules.toolsMap"):
        mapper.main()

    messages = "\n".join(r.message for r in caplog.records)
    assert "Mapping samples" in messages
    assert "100%" in messages


def test_mapsamples_still_populates_mapinfo_regardless_of_pool_completion_order(tmp_path):
    samplenames = [f"sample{i}" for i in range(3)]
    mapper = MapSamples(_make_args(tmp_path, samplenames))

    with patch("trnagraph.modules.toolsMap.MapReads.map_sample", side_effect=_fake_map_sample):
        mapper.main()

    with open(mapper.expinfo.mapinfo) as f:
        header = f.readline().split()
    assert set(header) == set(samplenames)
