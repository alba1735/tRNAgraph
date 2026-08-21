"""Regression tests for adataBuild.py's build-progress phase list assembly (roadmap.md Phase 2:
"tqdm" -- post-Stage-3 design follow-up). `_analysis_pipeline_phase_names()`/
`_full_build_phase_names()` compute the exact, fixed phase sequence a `toolsTG.PhaseTracker` needs
upfront -- everything here is knowable from the CLI args before any work starts (vst strategy,
readlengthsplit), so the shared tracker's percentage is accurate across the whole `analyze build`
command instead of just one class's slice of it. DESeq2 size factors are always computed now
(the `--nosizefactors` flag was removed -- confirmed broken, see docs/roadmap.md), so the phase
list no longer varies on that axis."""
from types import SimpleNamespace

from trnagraph.modules.adataBuild import _analysis_pipeline_phase_names, _full_build_phase_names


def test_analysis_pipeline_phases_include_deseq2_steps():
    assert _analysis_pipeline_phase_names() == [
        "Counting Reads", "Analyzing counts", "Counting Read Types",
        "Analyzing unique counts", "Generating Read Coverage plots",
    ]


def _args(**overrides):
    defaults = dict(vst="vst", readlengthsplit=None)
    defaults.update(overrides)
    return SimpleNamespace(**defaults)


def test_full_build_phases_simple_case_no_split():
    phases = _full_build_phase_names(_args())
    assert phases == [
        "Counting Reads", "Analyzing counts", "Counting Read Types", "Analyzing unique counts",
        "Generating Read Coverage plots", "Building AnnData object", "Computing VST", "Writing h5ad",
    ]


def test_full_build_phases_skip_vst_when_vst_strategy_is_none():
    phases = _full_build_phase_names(_args(vst="none"))
    assert "Computing VST" not in phases
    assert phases[-1] == "Writing h5ad"


def test_full_build_phases_repeat_analysis_block_twice_when_readlengthsplit_set():
    phases = _full_build_phase_names(_args(readlengthsplit=60))
    analysis_block = ["Counting Reads", "Analyzing counts", "Counting Read Types", "Analyzing unique counts", "Generating Read Coverage plots"]
    assert phases == (
        analysis_block + ["Building AnnData object", "Computing VST"]
        + analysis_block + analysis_block
        + ["Writing h5ad"]
    )
