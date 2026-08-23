"""Regression test for roadmap.md's "Add Threads flag back to Cluster/Graph" item. `analyze
cluster` had no --threads flag at all, and even where thread counts are meaningful (HDBSCAN's
core_dist_n_jobs is unaffected by a fixed --randomstate seed; UMAP's own n_jobs is only
overridden to 1 by UMAP itself when a seed IS set), nothing in adataCluster.py ever passed a
thread count into either library call."""
import logging
from unittest.mock import MagicMock, patch

import anndata as ad
import numpy as np
import pandas as pd

from trnagraph.modules.adataCluster import anndataCluster


def _make_cluster(threads=6):
    cluster = anndataCluster.__new__(anndataCluster)
    cluster.logger = logging.getLogger("trnagraph.modules.adataCluster")
    cluster.randomstate = 42
    cluster.mindist = 0.1
    cluster.variance_threshold = 0.0
    cluster.stats_metrics_umap = "euclidean"
    cluster.stats_metrics_hdbscan = "euclidean"
    cluster.threads = threads
    return cluster


def _make_adata(n_obs=8, n_var=5):
    rng = np.random.default_rng(0)
    return ad.AnnData(X=rng.normal(size=(n_obs, n_var)), obs=pd.DataFrame(index=[f"o{i}" for i in range(n_obs)]))


def test_adata_cluster_passes_threads_to_umap_n_jobs():
    cluster = _make_cluster(threads=6)
    adata = _make_adata()

    fake_umap_instance = MagicMock()
    fake_umap_instance.fit_transform.return_value = np.zeros((adata.n_obs, 2))
    fake_umap_instance.graph_ = MagicMock()
    fake_hdbscan_instance = MagicMock()
    fake_hdbscan_instance.fit_predict.return_value = np.zeros(adata.n_obs, dtype=int)

    with patch("trnagraph.modules.adataCluster.umap.UMAP", return_value=fake_umap_instance) as mock_umap, \
         patch("trnagraph.modules.adataCluster.hdbscan.HDBSCAN", return_value=fake_hdbscan_instance):
        cluster.adataCluster(adata, neighbors_plot=5, neighbors_cluster=5, min_samples=2, min_cluster_size=2, n_components=2)

    for call in mock_umap.call_args_list:
        assert call.kwargs.get("n_jobs") == 6


def test_adata_cluster_passes_threads_to_hdbscan_core_dist_n_jobs():
    cluster = _make_cluster(threads=6)
    adata = _make_adata()

    fake_umap_instance = MagicMock()
    fake_umap_instance.fit_transform.return_value = np.zeros((adata.n_obs, 2))
    fake_hdbscan_instance = MagicMock()
    fake_hdbscan_instance.fit_predict.return_value = np.zeros(adata.n_obs, dtype=int)

    with patch("trnagraph.modules.adataCluster.umap.UMAP", return_value=fake_umap_instance), \
         patch("trnagraph.modules.adataCluster.hdbscan.HDBSCAN", return_value=fake_hdbscan_instance) as mock_hdbscan:
        cluster.adataCluster(adata, neighbors_plot=5, neighbors_cluster=5, min_samples=2, min_cluster_size=2, n_components=2)

    assert mock_hdbscan.call_args.kwargs.get("core_dist_n_jobs") == 6
