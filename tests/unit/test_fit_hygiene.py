"""A dispersion trend is fitted across FEATURES, so a narrow feature axis cannot support one.

The threshold here is about matrix ROWS, not replicates: three replicates per group is standard
and was never the problem. What broke was a caller subsetting the feature axis hard enough to
hand PyDESeq2 a 9-sample x 1-feature matrix, where the dispersion prior's median-absolute-
deviation is taken over an empty slice and the Cox-Reid log-determinant degenerates. That run
emitted 78 `dispersion trend curve fitting did not converge` UserWarnings and 156 NumPy
RuntimeWarnings, then discarded every number it had computed.

Measured across 9 samples x N features, `fit_type='mean'` is clean from 2 features upward while
`parametric` warned at every size tried. PyDESeq2 already falls back to mean on non-convergence
and says so; choosing it up front is the same statistics without the announcement.
"""
import logging
import warnings
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd

import trnagraph.modules.toolsTG as toolsTG
from trnagraph.modules.toolsTG import adataLog2FC

READTYPE = "nreads_total_unique_norm"
RAW_READTYPE = "nreads_total_unique_raw"


def _adata(n_features, n_per_group=3, seed=0, runinfo=None):
    rng = np.random.default_rng(seed)
    groups = ["A", "B"]
    samples = [f"{g}_rep{i}" for g in groups for i in range(n_per_group)]
    rows = []
    for f in range(n_features):
        for sample in samples:
            raw = rng.negative_binomial(n=10, p=10 / (10 + 400))
            rows.append({"trna": f"trna{f}", "sample": sample,
                         "group": sample.rsplit("_rep", 1)[0],
                         RAW_READTYPE: raw, READTYPE: raw})
    obs = pd.DataFrame(rows)
    obs.index = [f"{r.trna}_{r['sample']}" for _, r in obs.iterrows()]
    adata = ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs)
    if runinfo is not None:
        adata.uns["trnagraphruninfo"] = runinfo
    return adata


def _fit_type_used(adata, **kwargs):
    log2fc = adataLog2FC(adata, compare="group", readtype=READTYPE, readcount_cutoff=0, **kwargs)
    with patch("trnagraph.modules.toolsTG.DeseqDataSet", wraps=toolsTG.DeseqDataSet) as dds:
        log2fc.log2fc_df()
    assert dds.call_args is not None, "expected a fit"
    return dds.call_args.kwargs.get("fit_type")


def test_a_narrow_feature_axis_forces_a_mean_trend():
    assert _fit_type_used(_adata(5)) == "mean"


def test_a_wide_feature_axis_keeps_the_parametric_trend():
    """The whole-object case -- hundreds of tRNAs -- must be untouched by this, or the change
    would silently alter every volcano and heatmap in the project."""
    assert _fit_type_used(_adata(toolsTG.MIN_FEATURES_FOR_PARAMETRIC_TREND)) == "parametric"


def test_a_single_feature_is_skipped_rather_than_fitted(caplog):
    """One feature cannot define a dispersion prior at all, so there is nothing to downgrade
    to -- the fit is skipped with a reason naming the subset."""
    adata = _adata(1)
    log2fc = adataLog2FC(adata, compare="group", readtype=READTYPE, readcount_cutoff=0)

    with patch("trnagraph.modules.toolsTG.DeseqDataSet") as dds, \
         caplog.at_level(logging.WARNING, logger="trnagraph.modules.toolsTG"):
        df, pairs = log2fc.log2fc_df()

    assert dds.call_args is None, "a one-feature matrix must never reach PyDESeq2"
    assert pairs == [("A", "B")]
    assert list(df.columns) == ["log2_A-B", "pval_A-B"], "callers index these unconditionally"
    assert "at least 2" in caplog.text


def test_the_skip_is_silent_about_a_genuinely_empty_result(caplog):
    """Zero surviving features is the ordinary cutoff case, not a too-narrow subset, and warning
    about it would fire on every run with a strict --cutoff."""
    adata = _adata(3)
    log2fc = adataLog2FC(adata, compare="group", readtype=READTYPE, readcount_cutoff=10_000)

    with caplog.at_level(logging.WARNING, logger="trnagraph.modules.toolsTG"):
        log2fc.log2fc_df()

    assert "at least" not in caplog.text


def test_a_narrow_fit_emits_no_convergence_warnings():
    """The point of the change, stated as behaviour rather than as a setting."""
    adata = _adata(6)
    log2fc = adataLog2FC(adata, compare="group", readtype=READTYPE, readcount_cutoff=0)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        log2fc.log2fc_df()

    messages = [str(w.message) for w in caught]
    assert not [m for m in messages if "dispersion trend curve fitting did not converge" in m]


# --- inheriting --dispfittype ---------------------------------------------------------

def test_the_fit_inherits_the_dispfittype_the_object_was_built_with():
    """`--dispfittype` reached the build-time DESeq2 and VST but never these fits, so one object
    got a parametric trend at build and a parametric-with-silent-fallback trend at graph time,
    with no way to say otherwise."""
    adata = _adata(toolsTG.MIN_FEATURES_FOR_PARAMETRIC_TREND,
                   runinfo={"flags": {"dispfittype": "mean"}})

    assert _fit_type_used(adata) == "mean"


def test_an_explicit_dispfittype_beats_the_recorded_one():
    adata = _adata(toolsTG.MIN_FEATURES_FOR_PARAMETRIC_TREND,
                   runinfo={"flags": {"dispfittype": "mean"}})

    assert _fit_type_used(adata, dispfittype="parametric") == "parametric"


def test_an_object_with_no_provenance_reads_as_the_pydeseq2_default():
    """An object built before the value was recorded, or a slim carrier holding only a frame,
    must keep behaving exactly as it did rather than failing to resolve."""
    assert toolsTG.recorded_dispfittype(_adata(2)) == "parametric"
    assert _fit_type_used(_adata(toolsTG.MIN_FEATURES_FOR_PARAMETRIC_TREND)) == "parametric"


def test_the_narrow_axis_override_beats_an_explicit_parametric_request():
    """A parametric trend through too few points is not a preference, it is a fit that cannot
    be made -- so the feature count wins over the flag rather than the other way round."""
    assert _fit_type_used(_adata(4), dispfittype="parametric") == "mean"


# --- warnings reach the log file ------------------------------------------------------

def test_library_warnings_are_written_to_the_run_log(tmp_path, monkeypatch):
    """Third-party warnings go through the `warnings` module to stderr, which the log file never
    saw. Someone reading an archived log therefore got a materially different story from someone
    watching the terminal -- which is how a several-hundred-line warning storm went unnoticed
    until a human happened to be looking at one."""
    from trnagraph import cli

    monkeypatch.chdir(tmp_path)
    with cli.handle_output(quiet=True, tool="unittest"):
        warnings.warn("a library said something", UserWarning)

    logs = list((tmp_path / ".log").glob("*_unittest.log"))
    assert len(logs) == 1
    assert "a library said something" in logs[0].read_text()


def test_capturing_is_released_afterwards(tmp_path, monkeypatch):
    """The handlers are shared objects and get closed on teardown, so leaving 'py.warnings'
    holding them would make the NEXT invocation write through a closed stream."""
    from trnagraph import cli

    monkeypatch.chdir(tmp_path)
    with cli.handle_output(quiet=True, tool="unittest"):
        pass

    assert logging.getLogger("py.warnings").handlers == []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        warnings.warn("after teardown", UserWarning)
    assert [str(w.message) for w in caught] == ["after teardown"]
