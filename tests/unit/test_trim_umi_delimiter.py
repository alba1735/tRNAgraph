"""Tests that fastp appends UMIs with the same delimiter `umi_tools` uses.

tRNAgraph has two UMI extraction paths and they had drifted apart:

  * `-u N`          -- fastp's own `--umi`, whose `--umi_delim` default is `:`
  * `-u N --umi3`   -- `umi_tools extract`, which uses `_`

`umi_tools dedup` (`preprocess map --dedup`) defaults to `--umi-separator=_`, and a colon is
additionally ambiguous because Illumina read names are already colon-delimited. Pinning fastp to
`_` makes both tRNAgraph paths produce one read-name shape, so a dataset trimmed by tRNAgraph is
always deduplicable without the user having to know which UMI path they took.

This asserts on the fastp argv because that is the boundary with the external binary; the
end-to-end proof is that a trimmed FASTQ maps and deduplicates, which the tutorial's yeast run
exercises for real.
"""
import pathlib
from types import SimpleNamespace

import pytest

from trnagraph.modules.toolsTrim import FastpTrimmer


def _write_fastq(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("@read\nACGT\n+\nIIII\n")
    return str(path)


def _make_trimmer(tmp_path, umilength, umi3):
    r1 = _write_fastq(tmp_path / "raw" / "SampleA_R1.fastq")
    manifest = tmp_path / "manifest.txt"
    manifest.write_text(f"SampleA\t{r1}\n")
    args = SimpleNamespace(
        input=str(manifest), adapter1=None, adapter2=None, length=15,
        umilength=umilength, umi3=umi3, threads=1, style=None,
        log=None, quiet=False, verbose=False,
    )
    return FastpTrimmer(args)


def _command_for(trimmer):
    prefix, files = next(iter(trimmer.samples.items()))
    cmd, _primary = trimmer._construct_command(prefix, files)
    return cmd


@pytest.mark.parametrize("umilength,umi3,expect_delim", [
    (7, False, True),   # fastp does the extraction -- delimiter must be pinned
    (0, False, False),  # no UMI requested -- no UMI flags at all
    (7, True, False),   # umi_tools does the extraction; fastp gets no UMI flags
])
def test_umi_delimiter_is_set_only_when_fastp_extracts_the_umi(
    tmp_path, umilength, umi3, expect_delim
):
    trimmer = _make_trimmer(tmp_path, umilength=umilength, umi3=umi3)

    cmd = _command_for(trimmer)

    if expect_delim:
        assert "--umi_delim" in cmd
        assert cmd[cmd.index("--umi_delim") + 1] == "_"
    else:
        assert "--umi_delim" not in cmd
