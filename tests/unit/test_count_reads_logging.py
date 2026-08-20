"""Regression tests for toolsCountReads.py's progress reporting (roadmap.md Phase 2: "tqdm"
stage 3) -- `countreads_main()` and `main()` previously drove their multiprocessing pools with a
blocking `Pool.map()` and no progress feedback at all for long-running per-sample read
classification/counting. Both now reuse the shared `toolsTG.progress_iterator()` helper (same
pattern as toolsTrim.py/toolsMap.py), verified here via pytest's caplog fixture against the
milestone-logging path (non-tty, the default under pytest's captured streams). The actual BAM
parsing (`getbamcounts`/`counttypereads`) is mocked out -- these tests exercise the progress
wiring and sample-to-result bookkeeping, not the counting algorithm itself, which is unchanged."""
import logging
from unittest.mock import patch

from trnagraph.modules import toolsCountReads
from trnagraph.modules.toolsCountReads import featurecount, counttypes


def _make_fixtures(tmp_path, samplenames):
    trnatable = tmp_path / "trnatable.txt"
    trnatable.write_text("tRNA1\tlocus1\tAla\tAGC\n")

    maturetrnas = tmp_path / "maturetRNAs.bed"
    maturetrnas.write_text("chr1\t0\t100\ttRNA1\t0\t+\n")

    trnaloci = tmp_path / "trnaloci.bed"
    trnaloci.write_text("chr1\t0\t100\tlocus1\t0\t+\n")

    samplefile = tmp_path / "samples.txt"
    lines = ["fastq\tsample\tgroup"]
    lines += [f"{name}.fastq\t{name}\t{name}" for name in samplenames]
    samplefile.write_text("\n".join(lines) + "\n")

    return {
        "trnatable": str(trnatable),
        "maturetrnas": [str(maturetrnas)],
        "trnaloci": [str(trnaloci)],
        "samplefile": str(samplefile),
    }


def _fake_getbamcounts(bamfile, samplename, *args, **kwargs):
    return featurecount(samplename, bamfile)


def _fake_counttypereads(bamfile, samplename, *args, **kwargs):
    return counttypes(samplename, bamfile)


def test_countreads_main_logs_progress_milestones_with_multiple_cores(tmp_path, caplog):
    samplenames = [f"sample{i}" for i in range(3)]
    fixtures = _make_fixtures(tmp_path, samplenames)
    countfile = tmp_path / "counts.txt"

    with patch("trnagraph.modules.toolsCountReads.getbamcounts", side_effect=_fake_getbamcounts), \
         caplog.at_level(logging.INFO, logger="trnagraph.modules.toolsCountReads"):
        toolsCountReads.countreads_main(
            bamdir=str(tmp_path), cores=2, countfile=str(countfile),
            genetypefile=None, bedfile=[], quiet=False, **fixtures,
        )

    messages = "\n".join(r.message for r in caplog.records)
    assert "Counting reads" in messages
    assert "100%" in messages


def test_countreads_main_quiet_still_emits_milestones_for_file_persistence(tmp_path, caplog):
    samplenames = [f"sample{i}" for i in range(3)]
    fixtures = _make_fixtures(tmp_path, samplenames)
    countfile = tmp_path / "counts.txt"

    with patch("trnagraph.modules.toolsCountReads.getbamcounts", side_effect=_fake_getbamcounts), \
         caplog.at_level(logging.INFO, logger="trnagraph.modules.toolsCountReads"):
        toolsCountReads.countreads_main(
            bamdir=str(tmp_path), cores=2, countfile=str(countfile),
            genetypefile=None, bedfile=[], quiet=True, **fixtures,
        )

    messages = "\n".join(r.message for r in caplog.records)
    assert "Counting reads" in messages
    assert "100%" in messages


def test_countreads_main_assigns_results_to_correct_sample_under_pool_mode(tmp_path):
    """Regression test for switching Pool.map() (order-preserving) to imap() wrapped in
    progress_iterator: the written count file's per-sample columns must still line up with the
    correct sample's counts afterwards. Each fake sample gets a distinct trnacount (10, 20, 30,
    40) so a column swap/misalignment would be caught by the assertion below; results cross the
    Pool's process boundary via their (pickled) return value, so this checks that value -- not a
    side channel, which wouldn't be visible back in the parent process."""
    samplenames = [f"sample{i}" for i in range(4)]
    fixtures = _make_fixtures(tmp_path, samplenames)
    countfile = tmp_path / "counts.txt"

    def _distinct_getbamcounts(bamfile, samplename, *args, **kwargs):
        result = featurecount(samplename, bamfile)
        index = samplenames.index(samplename)
        result.trnacounts["tRNA1"] = (index + 1) * 10
        return result

    with patch("trnagraph.modules.toolsCountReads.getbamcounts", side_effect=_distinct_getbamcounts):
        toolsCountReads.countreads_main(
            bamdir=str(tmp_path), cores=2, countfile=str(countfile), nofrag=True,
            genetypefile=None, bedfile=[], quiet=True, **fixtures,
        )

    lines = countfile.read_text().strip().split("\n")
    header = lines[0].split("\t")
    row = next(l for l in lines[1:] if l.startswith("tRNA1\t"))
    values = dict(zip(header, row.split("\t")[1:]))

    for index, name in enumerate(samplenames):
        assert values[name] == str((index + 1) * 10)


def test_counttypereads_main_logs_progress_milestones_with_multiple_cores(tmp_path, caplog):
    samplenames = [f"sample{i}" for i in range(3)]
    fixtures = _make_fixtures(tmp_path, samplenames)
    countfile = tmp_path / "typecounts.txt"
    realcountfile = tmp_path / "typerealcounts.txt"

    with patch("trnagraph.modules.toolsCountReads.counttypereads", side_effect=_fake_counttypereads), \
         caplog.at_level(logging.INFO, logger="trnagraph.modules.toolsCountReads"):
        toolsCountReads.main(
            bamdir=str(tmp_path), cores=2, countfile=str(countfile), realcountfile=str(realcountfile),
            bedfile=[], quiet=False, **fixtures,
        )

    messages = "\n".join(r.message for r in caplog.records)
    assert "Counting reads" in messages
    assert "100%" in messages


def test_counttypereads_main_quiet_still_emits_milestones(tmp_path, caplog):
    samplenames = [f"sample{i}" for i in range(3)]
    fixtures = _make_fixtures(tmp_path, samplenames)
    countfile = tmp_path / "typecounts.txt"
    realcountfile = tmp_path / "typerealcounts.txt"

    with patch("trnagraph.modules.toolsCountReads.counttypereads", side_effect=_fake_counttypereads), \
         caplog.at_level(logging.INFO, logger="trnagraph.modules.toolsCountReads"):
        toolsCountReads.main(
            bamdir=str(tmp_path), cores=2, countfile=str(countfile), realcountfile=str(realcountfile),
            bedfile=[], quiet=True, **fixtures,
        )

    messages = "\n".join(r.message for r in caplog.records)
    assert "Counting reads" in messages
    assert "100%" in messages
