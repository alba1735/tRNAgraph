"""Regression tests for toolsTG.adataLog2FC.log2fc_df (roadmap.md Phase 1: "Replace manual log2FC with native DESeq2 output")."""
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd

import trnagraph.modules.toolsTG as toolsTG
from trnagraph.modules.toolsTG import adataLog2FC

READTYPE = "nreads_total_unique_norm"
RAW_READTYPE = "nreads_total_unique_raw"


def _make_adata(trna_group_means, n_per_group=4, seed=0):
    """
    trna_group_means: {trna_name: {group_name: mean_raw_count}}. Builds a synthetic obs
    dataframe (one row per trna x sample) shaped like the real pipeline's adata.obs, with a
    NORMALIZED readtype column equal to the raw column (sizefactor=1, irrelevant here since
    PyDESeq2 refits its own size factors from raw counts regardless).
    """
    rng = np.random.default_rng(seed)
    groups = sorted({g for means in trna_group_means.values() for g in means})
    samples = [f"{g}_rep{i}" for g in groups for i in range(n_per_group)]
    sample_group = {s: s.rsplit("_rep", 1)[0] for s in samples}

    rows = []
    for trna, group_means in trna_group_means.items():
        for sample in samples:
            group = sample_group[sample]
            mean = group_means[group]
            raw = rng.negative_binomial(n=10, p=10 / (10 + mean))
            rows.append({"trna": trna, "sample": sample, "group": group, RAW_READTYPE: raw, READTYPE: raw})

    obs = pd.DataFrame(rows)
    obs.index = [f"{r.trna}_{r['sample']}" for _, r in obs.iterrows()]
    adata = ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs)
    return adata


def test_log2fc_df_keeps_pair_columns_when_no_trna_passes_cutoff():
    """
    pairs (and therefore the log2_<pair>/pval_<pair> column set) come from self.compare's
    column labels, independent of which/how many trna rows survive the readcount_cutoff
    filter. Callers (e.g. plotsVolcano.py) index those columns unconditionally, so an empty
    result must still have them -- a column-less empty df caused a KeyError downstream.
    """
    adata = _make_adata({"trnaLow": {"A": 5, "B": 6}}, n_per_group=3)
    log2fc = adataLog2FC(adata, compare="group", readtype=READTYPE, readcount_cutoff=80)

    df, pairs = log2fc.log2fc_df()

    assert pairs == [("A", "B")]
    assert list(df.columns) == ["log2_A-B", "pval_A-B"]
    assert len(df) == 0


def test_log2fc_df_detects_real_signal_and_reports_valid_pvalues():
    """
    Sanity check for the PyDESeq2-based replacement: a trna with a strong, consistent
    between-group difference should get a log2FC of the right sign/magnitude, a flat trna
    should get a log2FC near zero, and every padj must be a valid probability (or NaN, which
    PyDESeq2 uses for independent-filtering-excluded features -- never out of [0, 1]).
    """
    # DESeq2-style size-factor estimation assumes most features are NOT differentially
    # expressed, so it needs several stable ("flat") features alongside the one differential
    # feature under test -- with too few features and one wildly differential, the size
    # factors themselves get skewed by it, which would bias even the "flat" gene's apparent
    # fold change. This mirrors having a real multi-hundred-feature tRNA panel.
    trna_group_means = {"trnaHigh": {"A": 800, "B": 50}}  # ~16x, A > B
    trna_group_means.update({f"trnaFlat{i}": {"A": 400, "B": 400} for i in range(8)})
    adata = _make_adata(trna_group_means, n_per_group=5)
    log2fc = adataLog2FC(adata, compare="group", readtype=READTYPE, readcount_cutoff=80)

    df, pairs = log2fc.log2fc_df()

    assert pairs == [("A", "B")]
    assert set(df.index) == set(trna_group_means)

    # contrast=[condition, "B", "A"] -> log2FoldChange = log2(B/A); B << A so this is strongly negative.
    assert df.loc["trnaHigh", "log2_A-B"] < -1.5
    for i in range(8):
        assert abs(df.loc[f"trnaFlat{i}", "log2_A-B"]) < 1.0

    pvals = df["pval_A-B"].dropna()
    assert not pvals.empty
    assert ((pvals >= 0) & (pvals <= 1)).all()


def test_log2fc_df_pins_deseq2_to_a_single_cpu():
    """log2fc_df() runs inside adataGraph.py's own multiprocessing.Pool worker processes
    (plotsHeatmap.py and plotsVolcano.py both call it at graph time). DeseqDataSet defaults to
    n_cpus=None, which PyDESeq2 resolves to "all available CPUs" via joblib's loky backend --
    spawning a second real process pool from inside an already-forked worker, a known
    deadlock-prone pattern (confirmed live on this project: real trnagraph graph processes were
    found hung 24+ hours later with orphaned loky resource-tracker children attached). joblib's
    n_jobs=1 is the one setting that avoids creating any nested pool at all (it runs sequentially
    in-process, not just with a smaller pool), so log2fc_df() must always pin n_cpus=1."""
    # Two features, not one: a single-feature matrix is now skipped before the fit (a
    # dispersion prior cannot be estimated from it), so a one-tRNA fixture would never
    # construct a DeseqDataSet at all and this would pass without testing anything.
    adata = _make_adata({"trnaA": {"A": 400, "B": 400},
                         "trnaB": {"A": 300, "B": 300}}, n_per_group=3)
    log2fc = adataLog2FC(adata, compare="group", readtype=READTYPE, readcount_cutoff=80)

    with patch("trnagraph.modules.toolsTG.DeseqDataSet", wraps=toolsTG.DeseqDataSet) as mock_dds:
        log2fc.log2fc_df()

    assert mock_dds.call_args is not None
    assert mock_dds.call_args.kwargs.get("n_cpus") == 1


def _make_amino_adata(amino_group_means, n_per_group=3, seed=1):
    """As _make_adata, but each amino carries two tRNA isodecoders, so an amino-level fit has
    to aggregate them rather than pick one. This is the shape `-g compare` works at: it plots
    fold changes per amino acid and per isoacceptor, never per tRNA."""
    rng = np.random.default_rng(seed)
    groups = sorted({g for means in amino_group_means.values() for g in means})
    samples = [f"{g}_rep{i}" for g in groups for i in range(n_per_group)]
    sample_group = {s: s.rsplit("_rep", 1)[0] for s in samples}

    rows = []
    for amino, group_means in amino_group_means.items():
        for copy in (1, 2):
            for sample in samples:
                raw = rng.negative_binomial(n=10, p=10 / (10 + group_means[sample_group[sample]]))
                rows.append({
                    "trna": f"tRNA-{amino}-AGC-{copy}", "amino": amino, "iso": f"{amino}-AGC",
                    "sample": sample, "group": sample_group[sample],
                    RAW_READTYPE: raw, READTYPE: raw,
                })
    obs = pd.DataFrame(rows)
    obs.index = [f"{r.trna}_{r['sample']}" for _, r in obs.iterrows()]
    return ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs)


def test_log2fc_df_can_fit_on_a_feature_axis_other_than_trna():
    """`-g compare` works at amino/iso level while both pivot sites here were hardwired to
    index='trna', which is why compare had to keep its own hand-rolled t-test engine."""
    adata = _make_amino_adata({"Ala": {"A": 800, "B": 50}, "Arg": {"A": 400, "B": 400},
                               "Asn": {"A": 300, "B": 300}, "Asp": {"A": 250, "B": 250}})
    log2fc = adataLog2FC(adata, compare="group", readtype=READTYPE, readcount_cutoff=0)

    df, pairs = log2fc.log2fc_df(index_col="amino")

    assert pairs == [("A", "B")]
    assert sorted(df.index) == ["Ala", "Arg", "Asn", "Asp"], "expected one row per amino acid"
    # Ala is 800 -> 50 across the contrast, i.e. strongly down from A to B.
    assert df.loc["Ala", "log2_A-B"] < -1


def test_log2fc_df_still_defaults_to_the_trna_axis():
    adata = _make_amino_adata({"Ala": {"A": 800, "B": 50}, "Arg": {"A": 400, "B": 400},
                               "Asn": {"A": 300, "B": 300}, "Asp": {"A": 250, "B": 250}})
    log2fc = adataLog2FC(adata, compare="group", readtype=READTYPE, readcount_cutoff=0)

    df, _ = log2fc.log2fc_df()

    assert all(name.startswith("tRNA-") for name in df.index)
