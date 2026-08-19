"""Regression tests for anndataCluster.adataCombine's AnnData slot placement (roadmap.md Phase 1: "obs/obsm/obsp/uns mapping audit")."""
import anndata as ad
import numpy as np
import pandas as pd

from trnagraph.modules.adataCluster import anndataCluster
from trnagraph.modules.toolsSchemas import VariantTag


def _make_adata(n_obs=6):
    obs = pd.DataFrame(
        {
            "sample": [f"s{i}" for i in range(n_obs)],
            "group": ["A" if i % 2 == 0 else "B" for i in range(n_obs)],
            "trna": [f"trna{i % 3}" for i in range(n_obs)],
        },
        index=[f"obs{i}" for i in range(n_obs)],
    )
    return ad.AnnData(X=np.zeros((n_obs, 1)), obs=obs)


def _make_cluster(tag="full"):
    cluster = anndataCluster.__new__(anndataCluster)
    cluster.variant_spec = VariantTag(raw=f"norm:{tag}", norm="norm", tag=tag)
    return cluster


def test_sample_results_go_to_obsm_not_uns_full_variant():
    """
    'sample'-group cluster/UMAP results are obs-aligned (a subset of adata.obs_names, since
    low-coverage/'Und' samples get filtered before clustering) and must be stored in obsm,
    reindexed onto the full obs axis -- not uns, which wouldn't get resliced if the AnnData
    object is subset later, and previously left the data effectively write-only since nothing
    downstream read it back out of uns.
    """
    adata = _make_adata(n_obs=6)
    cluster = _make_cluster("full")

    # Simulate a subset (2 of 6 samples filtered out by adataPreprocess's coverage/Und filters)
    kept = adata.obs.index[:4]
    df = pd.DataFrame(
        {
            "standard_umap1": np.arange(4, dtype=float),
            "standard_umap2": np.arange(4, dtype=float) * 2,
            "cluster_umap1": np.arange(4, dtype=float) * 3,
            "cluster_hdbscan": [0, 0, 1, -1],
        },
        index=kept,
    )

    result = cluster.adataCombine(adata, df, "sample")

    assert "sample_cluster_umap" in result.obsm
    assert "sample_cluster_umap" not in result.uns
    obsm_df = result.obsm["sample_cluster_umap"]
    assert obsm_df.shape[0] == adata.n_obs  # reindexed onto full obs axis
    assert obsm_df.loc[kept, "cluster_hdbscan"].tolist() == [0, 0, 1, -1]
    assert obsm_df.loc[adata.obs.index[4:], "cluster_hdbscan"].isna().all()  # filtered-out rows are NaN, not dropped

    # Scalar per-obs summary columns still land in .obs as before.
    assert "sample_cluster" in result.obs.columns
    assert "sample_umap1" in result.obs.columns


def test_sample_results_go_to_namespaced_obsm_for_split_variant():
    adata = _make_adata(n_obs=4)
    cluster = _make_cluster("u60")

    df = pd.DataFrame(
        {
            "standard_umap1": np.arange(4, dtype=float),
            "standard_umap2": np.arange(4, dtype=float),
            "cluster_umap1": np.arange(4, dtype=float),
            "cluster_hdbscan": [0, 0, 1, 1],
        },
        index=adata.obs.index,
    )

    result = cluster.adataCombine(adata, df, "sample")

    assert "sample_cluster_umap_u60" in result.obsm
    assert "sample_cluster_umap" not in result.obsm
    split_uns = result.uns.get("size_splits", {}).get("u60", {})
    assert "sample_cluster_umap" not in split_uns


def test_group_results_stay_in_uns_since_they_are_not_obs_aligned():
    """
    'group'-group results collapse onto a trna x group axis with different cardinality than
    adata's obs axis, so they genuinely can't be represented in obsm and belong in uns.
    """
    adata = _make_adata(n_obs=6)
    cluster = _make_cluster("full")

    # Index is 'trna_group' composite keys -- NOT a subset of adata.obs_names.
    df = pd.DataFrame(
        {
            "standard_umap1": [0.0, 1.0],
            "standard_umap2": [0.0, 1.0],
            "cluster_umap1": [0.0, 1.0],
            "cluster_hdbscan": [0, 1],
        },
        index=["trna0_A", "trna1_B"],
    )

    result = cluster.adataCombine(adata, df, "group")

    assert "group_cluster_umap" in result.uns
    assert "group_cluster_umap" not in result.obsm
    pd.testing.assert_frame_equal(result.uns["group_cluster_umap"], df)


def test_stale_uns_entry_is_cleared_on_reclustering():
    """A pre-existing uns['sample_cluster_umap'] from before this moved to obsm should be
    dropped when re-clustering, instead of lingering as an orphaned duplicate."""
    adata = _make_adata(n_obs=3)
    adata.uns["sample_cluster_umap"] = pd.DataFrame({"stale": [1, 2, 3]})
    cluster = _make_cluster("full")

    df = pd.DataFrame(
        {
            "standard_umap1": [0.0, 1.0, 2.0],
            "standard_umap2": [0.0, 1.0, 2.0],
            "cluster_umap1": [0.0, 1.0, 2.0],
            "cluster_hdbscan": [0, 0, 1],
        },
        index=adata.obs.index,
    )

    result = cluster.adataCombine(adata, df, "sample")

    assert "sample_cluster_umap" not in result.uns
    assert "sample_cluster_umap" in result.obsm
