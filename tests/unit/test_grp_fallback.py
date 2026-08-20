"""Regression tests for roadmap.md Phase 2 "Parameter Fallback": graphing parameter columns
that don't exist in adata.obs should fall back to 'sample' with a warning, not raise ValueError."""
import anndata as ad
import numpy as np
import pandas as pd
from types import SimpleNamespace

from trnagraph.modules import toolsTG
from trnagraph.modules.adataGraph import anndataGrapher


def _make_adata(n_obs=4):
    obs = pd.DataFrame(
        {
            "sample": [f"s{i}" for i in range(n_obs)],
            "group": ["A" if i % 2 == 0 else "B" for i in range(n_obs)],
        },
        index=[f"obs{i}" for i in range(n_obs)],
    )
    return ad.AnnData(X=np.zeros((n_obs, 1)), obs=obs)


def _make_grapher(args, adata):
    grapher = anndataGrapher.__new__(anndataGrapher)
    grapher.args = args
    grapher.adata = adata
    return grapher


def test_resolve_grp_column_passes_through_existing_column():
    adata = _make_adata()
    assert toolsTG.resolve_grp_column(adata, "group", "comparegrp1") == "group"


def test_resolve_grp_column_falls_back_to_sample_and_warns(caplog):
    import logging
    adata = _make_adata()
    with caplog.at_level(logging.WARNING, logger="trnagraph.modules.toolsTG"):
        result = toolsTG.resolve_grp_column(adata, "nonexistent", "comparegrp1")
    assert result == "sample"
    warning_text = "\n".join(r.message for r in caplog.records)
    assert "nonexistent" in warning_text
    assert "comparegrp1" in warning_text
    assert "sample" in warning_text


def test_resolve_grp_column_respects_custom_default():
    adata = _make_adata()
    result = toolsTG.resolve_grp_column(adata, "nonexistent", "coveragegrp", default="group")
    assert result == "group"


def test_anndata_grapher_resolves_all_grp_args_up_front():
    """
    self.args.covgrp/comparegrp1/comparegrp2/heatgrp/volgrp are all read directly, multiple
    times, by anndataGrapher (the log2FC precompute in __init__, cmap_dict, and each plot
    module call in plot()) -- resolving only inside the individual plots*.py modules leaves
    the __init__-time precompute step (which runs first and reads the raw, unvalidated
    self.args values) still crashing on a typo'd --heatgrp/--volgrp. This must be resolved
    once, centrally, so every call site sees the same corrected value.
    """
    adata = _make_adata()
    args = SimpleNamespace(
        covgrp="typo_cov", comparegrp1="typo_cmp1", comparegrp2="group",
        heatgrp="typo_heat", volgrp="group",
    )
    grapher = _make_grapher(args, adata)
    grapher._resolve_grp_args()

    assert args.covgrp == "sample"
    assert args.comparegrp1 == "sample"
    assert args.comparegrp2 == "group"
    assert args.heatgrp == "sample"
    assert args.volgrp == "group"
