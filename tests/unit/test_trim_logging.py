"""Regression tests for toolsTrim.py's logging (roadmap.md Phase 2: "Trimmer log/quiet" +
"Logging") -- fastp's per-sample stderr was previously discarded entirely on success, only ever
surfaced via print() on failure. FastpTrimmer now logs via a stable logging.getLogger(__name__)
logger (no handlers of its own -- see test_cli_logging.py's configure_logging(), the shared
place --log/--quiet handlers get attached), verified here via pytest's caplog fixture, which
captures emitted log records independent of handler configuration."""
import logging
from types import SimpleNamespace
from unittest.mock import patch

from trnagraph.modules.toolsTrim import FastpTrimmer


def _make_manifest(tmp_path):
    r1 = tmp_path / "sample_R1.fastq"
    r1.write_text("@read\nACGT\n+\nIIII\n")
    manifest = tmp_path / "manifest.txt"
    manifest.write_text(f"sample1\t{r1}\n")
    return str(manifest)


def _make_args(tmp_path):
    return SimpleNamespace(
        input=_make_manifest(tmp_path), adapter1=None, adapter2=None, length=15,
        umilength=0, umi3=False, threads=1, colormap=None,
        log=None, quiet=False, verbose=False,
    )


class _FakeCompletedProcess:
    def __init__(self, returncode, stderr):
        self.returncode = returncode
        self.stderr = stderr


def test_run_process_returns_fastp_stderr_on_success_instead_of_discarding_it(tmp_path):
    trimmer = FastpTrimmer(_make_args(tmp_path))

    with patch("subprocess.run", return_value=_FakeCompletedProcess(0, "fastp summary stats here")):
        result = trimmer._run_process("sample1", {"r1": "r1.fastq", "r2": None})

    assert result == ("sample1", True, "fastp summary stats here")


def test_process_logs_fastp_stderr_at_info_level_on_success(tmp_path, caplog):
    trimmer = FastpTrimmer(_make_args(tmp_path))

    with patch("subprocess.run", return_value=_FakeCompletedProcess(0, "fastp summary stats here")), \
         patch.object(FastpTrimmer, "_generate_summary", lambda self: None), \
         caplog.at_level(logging.INFO, logger="trnagraph.modules.toolsTrim"):
        trimmer.process()

    info_messages = "\n".join(r.message for r in caplog.records if r.levelno == logging.INFO)
    assert "fastp summary stats here" in info_messages


def test_process_logs_failure_at_error_level_with_stderr(tmp_path, caplog):
    trimmer = FastpTrimmer(_make_args(tmp_path))

    with patch("subprocess.run", return_value=_FakeCompletedProcess(1, "fastp exploded")), \
         patch.object(FastpTrimmer, "_generate_summary", lambda self: None), \
         caplog.at_level(logging.INFO, logger="trnagraph.modules.toolsTrim"):
        trimmer.process()

    error_records = [r for r in caplog.records if r.levelno == logging.ERROR]
    assert any("fastp exploded" in r.message for r in error_records)


def test_construct_command_omits_umi_flags_from_fastp_for_umi3():
    """Regression test for roadmap.md's "Trimmer UMI Handling" item -- fastp has no
    tail-anchored --umi_loc value, so --umi3 must not pass ANY --umi/--umi_loc flags to fastp at
    all (previously this branch was a bare `pass`, a silent no-op that looked like it worked but
    extracted nothing). The actual 3' extraction happens via a separate umi_tools post-processing
    step -- see the tests below."""
    args = SimpleNamespace(umilength=6, umi3=True, adapter1=None, adapter2=None, length=15, threads=1)
    trimmer = FastpTrimmer.__new__(FastpTrimmer)
    trimmer.args = args
    trimmer.fastp_threads = 1

    cmd, _ = trimmer._construct_command("prefix", {"r1": "r1.fastq", "r2": None})

    assert "--umi" not in cmd
    assert "--umi_loc" not in cmd
    assert "--umi_len" not in cmd


def test_construct_command_keeps_umi_flags_for_default_5prime_case():
    args = SimpleNamespace(umilength=6, umi3=False, adapter1=None, adapter2=None, length=15, threads=1)
    trimmer = FastpTrimmer.__new__(FastpTrimmer)
    trimmer.args = args
    trimmer.fastp_threads = 1

    cmd, _ = trimmer._construct_command("prefix", {"r1": "r1.fastq", "r2": None})

    assert "--umi" in cmd
    assert cmd[cmd.index("--umi_loc") + 1] == "read1"


def test_run_process_invokes_umi_tools_extract_3prime_after_fastp_for_umi3(tmp_path):
    """The real bug this item fixes: for --umi3, a umi_tools extract --3prime post-processing
    step must actually run on fastp's own trimmed/merged output (mirroring tRAX's original
    approach), replacing it in place -- not a silent no-op."""
    args = _make_args(tmp_path)
    args.umilength = 6
    args.umi3 = True
    trimmer = FastpTrimmer(args)

    fastp_result = _FakeCompletedProcess(0, "fastp stats")
    umi_result = _FakeCompletedProcess(0, "umi_tools stats")

    with patch("subprocess.run", side_effect=[fastp_result, umi_result]) as mock_run, \
         patch("os.replace") as mock_replace:
        result = trimmer._run_process("sample1", {"r1": "r1.fastq", "r2": None})

    assert result == ("sample1", True, "fastp statsumi_tools stats")
    assert mock_run.call_count == 2
    umi_cmd = mock_run.call_args_list[1].args[0]
    assert umi_cmd[:3] == ["umi_tools", "extract", "--3prime"]
    assert "NNNNNN" in umi_cmd
    mock_replace.assert_called_once()


def test_run_process_does_not_invoke_umi_tools_for_default_5prime_case(tmp_path):
    args = _make_args(tmp_path)
    args.umilength = 6
    args.umi3 = False
    trimmer = FastpTrimmer(args)

    with patch("subprocess.run", return_value=_FakeCompletedProcess(0, "fastp stats")) as mock_run:
        result = trimmer._run_process("sample1", {"r1": "r1.fastq", "r2": None})

    assert result == ("sample1", True, "fastp stats")
    assert mock_run.call_count == 1


def test_run_process_reports_failure_when_umi_tools_extraction_fails(tmp_path):
    args = _make_args(tmp_path)
    args.umilength = 6
    args.umi3 = True
    trimmer = FastpTrimmer(args)

    fastp_result = _FakeCompletedProcess(0, "fastp stats")
    umi_result = _FakeCompletedProcess(1, "umi_tools exploded")

    with patch("subprocess.run", side_effect=[fastp_result, umi_result]):
        result = trimmer._run_process("sample1", {"r1": "r1.fastq", "r2": None})

    assert result == ("sample1", False, "fastp statsumi_tools exploded")


def test_repeated_construction_reuses_the_same_module_logger(tmp_path):
    """
    Regression test for the fixed leak: FastpTrimmer must use a stable, module-scoped logger
    name (logging.getLogger(__name__)), not one keyed on id(self) -- the logging module never
    garbage-collects loggers, so an instance-keyed name would register a new logger (and leak
    any handler it held) on every construction.
    """
    trimmer1 = FastpTrimmer(_make_args(tmp_path))
    trimmer2 = FastpTrimmer(_make_args(tmp_path))

    assert trimmer1.logger is trimmer2.logger
    assert trimmer1.logger.name == "trnagraph.modules.toolsTrim"
