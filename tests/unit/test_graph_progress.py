"""Regression tests for adataGraph.py's progress reporting (roadmap.md Phase 2: "tqdm" graphing
objective). `anndataGrapher.main()` previously drove its outer cross-graph-type multiprocessing
pool with a blocking `pool.starmap()` and no progress feedback at all. It now builds a shared
`toolsTG.PhaseTracker` (one phase per graph type, weighted by a best-effort upfront plot-count
estimate -- coverage can produce hundreds/thousands of plots, most other types only a handful),
switches the outer pool to an ORDERED `pool.imap()` (not imap_unordered -- phase-tracker
completion must be registered in the same order the phases were declared) wrapped in a
`with self.phase_tracker.phase(): pass` per completed type, and threads the tracker into
coverage's own visualizer so its per-plot completions tick the outer bar in lockstep instead of
only completing atomically."""
import logging
from types import SimpleNamespace
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd

from trnagraph.modules.adataGraph import anndataGrapher


def _make_adata(n_obs=6):
    obs = pd.DataFrame(
        {
            "sample": [f"s{i}" for i in range(n_obs)],
            "group": ["A", "B", "C"] * (n_obs // 3),
            "trna": [f"tRNA{i}" for i in range(n_obs)],
        },
        index=[f"obs{i}" for i in range(n_obs)],
    )
    return ad.AnnData(X=np.zeros((n_obs, 1)), obs=obs)


def _make_grapher(args, adata):
    grapher = anndataGrapher.__new__(anndataGrapher)
    grapher.logger = logging.getLogger("trnagraph.modules.adataGraph")
    grapher.args = args
    grapher.adata = adata
    grapher.quiet = getattr(args, "quiet", False)
    return grapher


def test_compute_graph_weight_coverage_uses_covobs_nunique():
    adata = _make_adata(n_obs=6)
    args = SimpleNamespace(covobs="trna", pcareadtypes=[], radarmethod=[], diffrts=[], volgrp="group",
                           heatgrp="group", clustergrp="group", comparegrp2="group", logogrp="group")
    grapher = _make_grapher(args, adata)

    assert grapher._compute_graph_weight("coverage") == 6


def test_compute_graph_weight_pca_uses_readtype_count_plus_extras():
    adata = _make_adata()
    args = SimpleNamespace(covobs="trna", pcareadtypes=["a", "b", "c"], radarmethod=[], diffrts=[],
                           volgrp="group", heatgrp="group", clustergrp="group", comparegrp2="group", logogrp="group")
    grapher = _make_grapher(args, adata)

    assert grapher._compute_graph_weight("pca") == 5  # len(pcareadtypes) + 2


def test_compute_graph_weight_never_raises_and_falls_back(caplog):
    """Weight computation is a rough estimator, never allowed to block the actual graphing run --
    an unexpected/missing column must fall back to a modest constant, not raise."""
    adata = _make_adata()
    args = SimpleNamespace(covobs="does_not_exist", pcareadtypes=[], radarmethod=[], diffrts=[],
                           volgrp="group", heatgrp="group", clustergrp="group", comparegrp2="group", logogrp="group")
    grapher = _make_grapher(args, adata)

    weight = grapher._compute_graph_weight("coverage")
    assert isinstance(weight, int) and weight >= 1


def _fake_plot(self, gt, threaded=None):
    return f"{gt} done" if threaded else None


_fake_plot.__name__ = "plot"


def test_main_builds_phase_tracker_covering_every_graphtype_single_threaded(caplog):
    adata = _make_adata()
    args = SimpleNamespace(
        graphtypes=["pca", "count"], threads=1, verbose=False, covobs="trna", covgrp=None,
        pcareadtypes=["a"], radarmethod=[], diffrts=[], volgrp="group", heatgrp="group",
        clustergrp="group", comparegrp2="group", logogrp="group", quiet=False,
    )
    grapher = _make_grapher(args, adata)

    with patch.object(anndataGrapher, "plot", _fake_plot), \
         caplog.at_level(logging.INFO, logger="trnagraph.modules.adataGraph"):
        grapher.main()

    assert grapher.phase_tracker.phases == ["pca", "count"]
    phase_messages = [r.message for r in caplog.records if r.message.startswith("Graphing phase")]
    assert len(phase_messages) == 2
    assert "100%" in phase_messages[-1]


def test_main_passes_shared_phase_tracker_into_coverage_visualizer():
    adata = _make_adata()
    args = SimpleNamespace(
        graphtypes=["coverage"], threads=1, verbose=False, covobs="trna", covgrp=None, covtype="coverage",
        covgap=[False], covmethod="mean", combinedpdfonly=True, colormap=None, output="/tmp/graphtest",
        pcareadtypes=[], radarmethod=[], diffrts=[], volgrp="group", heatgrp="group",
        clustergrp="group", comparegrp2="group", logogrp="group", quiet=False,
    )
    grapher = _make_grapher(args, adata)
    grapher.cmap_dict = {}
    grapher.variant_spec = SimpleNamespace(raw="norm:full")

    captured = {}

    class _FakeVisualizer:
        def __init__(self, *a, **kwargs):
            captured["phase_tracker"] = kwargs.get("phase_tracker")
            captured["quiet"] = kwargs.get("quiet")

        def generate_split(self):
            pass

        def generate_combine(self):
            pass

    with patch("trnagraph.modules.adataGraph.plotsCoverage.visualizer", _FakeVisualizer), \
         patch("trnagraph.modules.adataGraph.toolsTG.builder", lambda *a, **k: ""):
        grapher.main()

    assert captured["phase_tracker"] is grapher.phase_tracker
    assert captured["quiet"] is False
