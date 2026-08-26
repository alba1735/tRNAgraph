"""Characterization tests for `makedb`'s per-organism covariance model selection.

These pin behavior that already worked, ahead of a refactor: the model paths used to be built
from `tRNADatabaseBuilder.script_dir` (its own package-relative `__file__` walk, duplicating what
`toolsTG.assets_dir()` now does) and were chosen by an inline if/elif chain buried in the middle
of `main()`, so nothing could check the mapping without running a full database build.

The expected filenames below are transcribed from that pre-refactor chain, which is the ground
truth here: the point of these assertions is that moving the lookup onto the shared helper does
not silently change which model an organism mode gets. That matters more than usual because
docs/roadmap.md records `arch` and `mito` as modes no run has ever exercised end to end -- a
wrong model there would not be caught by the bacterial test suite.
"""
from pathlib import Path

import pytest

from trnagraph.modules import toolsTDatabase, toolsTG


EXPECTED = {
    # orgmode: (mature model, tRNA model, prokaryotic mode)
    "euk": ("trnamature-euk.cm", "TRNAinf-euk.cm", False),
    "arch": ("trnamature-arch.cm", "TRNAinf-arch.cm", True),
    "mito": ("TRNAMatureMitoinf.cm", "TRNAinf.cm", False),
    "bact": ("trnamature-bact.cm", "TRNAinf-bact.cm", True),
}


@pytest.mark.parametrize("orgmode", sorted(EXPECTED))
def test_covariance_models_match_the_documented_per_domain_mapping(orgmode):
    mature, trna, prok_mode = toolsTDatabase._covariance_models(orgmode)
    expected_mature, expected_trna, expected_prok = EXPECTED[orgmode]
    assert Path(mature).name == expected_mature
    assert Path(trna).name == expected_trna
    assert prok_mode is expected_prok


@pytest.mark.parametrize("orgmode", sorted(EXPECTED))
def test_covariance_models_resolve_to_packaged_files(orgmode):
    """Both models must be real files inside the installed package's `assets/cm/`, so `makedb`
    works from a wheel install and not only from a source checkout."""
    cm_dir = Path(toolsTG.assets_dir()) / "cm"
    for path in toolsTDatabase._covariance_models(orgmode)[:2]:
        assert Path(path).parent == cm_dir
        assert Path(path).is_file(), f"missing packaged covariance model: {path}"


def test_covariance_models_rejects_an_unknown_orgmode():
    """The inline chain left `mature_model`/`trna_model` unbound for an unrecognised mode, which
    would surface as an UnboundLocalError far from the cause. The constructor already validates
    orgmode, so this is defence in depth rather than a reachable path today."""
    with pytest.raises(ValueError, match="orgmode"):
        toolsTDatabase._covariance_models("plant")
