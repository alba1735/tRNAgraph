"""Tests for `plotsMismatch`, the native port of tRAX's `boxplotmismatches.R` and
`plotreadmismatches.R`.

The rate deliberately diverges from tRAX on one point. tRAX computes
`mismatchedbases / (coverage + 10)` from the *size-factor-normalized* coverage table
(`boxplotmismatches.R:57`, reading `-coverage.txt`), while gating on the raw base-composition
columns. The size factor cancels in `m/c` but not in `m/(c + 10)`, so the pseudocount acts
like `10 x sizefactor` in raw terms -- a range of roughly 4 to 41 across the hg38 dataset's
size factors. Since the aggregation is a maximum across samples, that systematically favors
the lowest-size-factor sample and makes the y-axis partly a function of library size.
tRNAgraph computes the rate from raw counts instead, where the pseudocount means the same
thing in every sample.
"""
import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")

from trnagraph.modules import plotsMismatch


COVERAGE_TYPES = ["coverage", "mismatchedbases", "deletedbases",
                  "adenines", "thymines", "cytosines", "guanines"]


def _build_adata(*, positions=("1", "2"), trnas=("tRNA-Ala-AGC-1", "tRNA-Gly-GCC-1"),
                 samples=("s1", "s2"), values=None, raw_values=None, gaps=None):
    """Assemble a minimal AnnData shaped like a real build: obs = tRNA x sample."""
    obs_rows = [(t, s) for t in trnas for s in samples]
    obs = pd.DataFrame({
        "trna": [t for t, _ in obs_rows],
        "sample": [s for _, s in obs_rows],
        "group": ["g1" for _ in obs_rows],
        "amino": [t.split("-")[1] for t, _ in obs_rows],
        "iso": [t.split("-")[2] for t, _ in obs_rows],
    }, index=[f"{s}_{t}" for t, s in obs_rows])

    var = pd.DataFrame({
        "positions": [p for c in COVERAGE_TYPES for p in positions],
        "coverage": [c for c in COVERAGE_TYPES for _ in positions],
        "gap": [(gaps or {}).get(p, False) for _ in COVERAGE_TYPES for p in positions],
    }, index=[f"{c}_{p}" for c in COVERAGE_TYPES for p in positions])

    shape = (len(obs), len(var))
    raw = np.zeros(shape) if raw_values is None else np.asarray(raw_values, dtype=float)
    norm = np.full(shape, -999.0) if values is None else np.asarray(values, dtype=float)

    adata = anndata.AnnData(X=norm, obs=obs, var=var)
    adata.layers["raw"] = raw
    return adata


def _set(adata, obs_label, coverage_type, position, value):
    col = f"{coverage_type}_{position}"
    adata.layers["raw"][adata.obs_names.get_loc(obs_label), adata.var_names.get_loc(col)] = value


def _passing_bases(adata, obs_label, position, total=100):
    """Put enough base-composition counts behind a position to clear tRAX's >50 gate."""
    for base in ("adenines", "thymines", "cytosines", "guanines"):
        _set(adata, obs_label, base, position, total / 4)


# --- rate computation ---------------------------------------------------------------

def test_rate_is_raw_mismatches_over_raw_coverage_plus_pseudocount():
    adata = _build_adata()
    _set(adata, "s1_tRNA-Ala-AGC-1", "coverage", "1", 90)
    _set(adata, "s1_tRNA-Ala-AGC-1", "mismatchedbases", "1", 30)
    _passing_bases(adata, "s1_tRNA-Ala-AGC-1", "1")

    df = plotsMismatch.position_mismatch_rates(adata, pseudocount=10)

    row = df[(df["trna"] == "tRNA-Ala-AGC-1") & (df["position"] == "1")]
    assert row["rate"].iloc[0] == pytest.approx(30 / 100)


def test_rate_ignores_the_normalized_layer():
    """X carries normalized values; the rate must not read them."""
    adata = _build_adata()
    _set(adata, "s1_tRNA-Ala-AGC-1", "coverage", "1", 90)
    _set(adata, "s1_tRNA-Ala-AGC-1", "mismatchedbases", "1", 30)
    _passing_bases(adata, "s1_tRNA-Ala-AGC-1", "1")

    df = plotsMismatch.position_mismatch_rates(adata, pseudocount=10)

    assert (df["rate"] >= 0).all()
    assert df["rate"].max() == pytest.approx(0.3)


def test_rate_is_the_maximum_across_samples():
    adata = _build_adata()
    for sample, (cov, mis) in {"s1": (990, 10), "s2": (90, 30)}.items():
        label = f"{sample}_tRNA-Ala-AGC-1"
        _set(adata, label, "coverage", "1", cov)
        _set(adata, label, "mismatchedbases", "1", mis)
        _passing_bases(adata, label, "1")

    df = plotsMismatch.position_mismatch_rates(adata, pseudocount=10)

    row = df[(df["trna"] == "tRNA-Ala-AGC-1") & (df["position"] == "1")]
    assert row["rate"].iloc[0] == pytest.approx(0.3)
    assert len(row) == 1, "samples collapse into one row per (feature, position)"


def test_pseudocount_is_honoured():
    adata = _build_adata()
    _set(adata, "s1_tRNA-Ala-AGC-1", "coverage", "1", 90)
    _set(adata, "s1_tRNA-Ala-AGC-1", "mismatchedbases", "1", 9)
    _passing_bases(adata, "s1_tRNA-Ala-AGC-1", "1")

    lo = plotsMismatch.position_mismatch_rates(adata, pseudocount=10)
    hi = plotsMismatch.position_mismatch_rates(adata, pseudocount=910)

    assert lo["rate"].max() == pytest.approx(9 / 100)
    assert hi["rate"].max() == pytest.approx(9 / 1000)


# --- filtering ----------------------------------------------------------------------

def test_positions_below_the_base_composition_gate_are_dropped():
    """boxplotmismatches.R:110 -- adenines+thymines+cytosines+guanines > 50."""
    adata = _build_adata()
    for position, total in (("1", 100), ("2", 40)):
        _set(adata, "s1_tRNA-Ala-AGC-1", "coverage", position, 90)
        _set(adata, "s1_tRNA-Ala-AGC-1", "mismatchedbases", position, 30)
        _passing_bases(adata, "s1_tRNA-Ala-AGC-1", position, total=total)

    df = plotsMismatch.position_mismatch_rates(adata, pseudocount=10)
    kept = df[df["trna"] == "tRNA-Ala-AGC-1"]["position"].tolist()

    assert "1" in kept
    assert "2" not in kept


def test_gap_positions_are_dropped():
    adata = _build_adata(gaps={"2": True})
    for position in ("1", "2"):
        _set(adata, "s1_tRNA-Ala-AGC-1", "coverage", position, 90)
        _set(adata, "s1_tRNA-Ala-AGC-1", "mismatchedbases", position, 30)
        _passing_bases(adata, "s1_tRNA-Ala-AGC-1", position)

    df = plotsMismatch.position_mismatch_rates(adata, pseudocount=10)

    assert "2" not in df["position"].tolist()


def test_amino_and_anticodon_travel_with_each_row():
    adata = _build_adata()
    _set(adata, "s1_tRNA-Gly-GCC-1", "coverage", "1", 90)
    _set(adata, "s1_tRNA-Gly-GCC-1", "mismatchedbases", "1", 30)
    _passing_bases(adata, "s1_tRNA-Gly-GCC-1", "1")

    df = plotsMismatch.position_mismatch_rates(adata, pseudocount=10)
    row = df[df["trna"] == "tRNA-Gly-GCC-1"].iloc[0]

    assert row["amino"] == "Gly"
    assert row["iso"] == "GCC"


def test_position_order_follows_the_var_axis():
    """var is already in Sprinzl order, so it is the ordering authority -- not sorted()."""
    adata = _build_adata(positions=("10", "9"))
    for position in ("10", "9"):
        _set(adata, "s1_tRNA-Ala-AGC-1", "coverage", position, 90)
        _set(adata, "s1_tRNA-Ala-AGC-1", "mismatchedbases", position, 30)
        _passing_bases(adata, "s1_tRNA-Ala-AGC-1", position)

    df = plotsMismatch.position_mismatch_rates(adata, pseudocount=10)

    assert list(df["position"].cat.categories) == ["10", "9"]


def test_deletions_read_the_deletedbases_coverage_type():
    adata = _build_adata()
    _set(adata, "s1_tRNA-Ala-AGC-1", "coverage", "1", 90)
    _set(adata, "s1_tRNA-Ala-AGC-1", "mismatchedbases", "1", 0)
    _set(adata, "s1_tRNA-Ala-AGC-1", "deletedbases", "1", 30)
    _passing_bases(adata, "s1_tRNA-Ala-AGC-1", "1")

    mismatches = plotsMismatch.position_mismatch_rates(adata, pseudocount=10)
    deletions = plotsMismatch.position_mismatch_rates(adata, pseudocount=10,
                                                      numerator="deletedbases")

    assert mismatches["rate"].max() == pytest.approx(0.0)
    assert deletions["rate"].max() == pytest.approx(0.3)


# --- read-level histogram -----------------------------------------------------------

def test_histogram_is_skipped_when_the_uns_key_is_absent(caplog):
    """Objects built before uns['mismatch_counts'] existed must not fail the graph run."""
    adata = _build_adata()

    result = plotsMismatch.mismatch_count_shares(adata)

    assert result is None


def test_histogram_returns_per_sample_shares():
    adata = _build_adata()
    adata.uns["mismatch_counts"] = pd.DataFrame({
        "count": [0, 1, 0, 1],
        "type": ["trna", "trna", "nontrna", "nontrna"],
        "s1": [75, 25, 40, 10],
        "s2": [30, 10, 5, 5],
    })

    df = plotsMismatch.mismatch_count_shares(adata)

    trna_s1 = df[(df["type"] == "trna") & (df["Sample"] == "s1")]
    assert trna_s1["share"].tolist() == pytest.approx([0.75, 0.25])
    nontrna_s2 = df[(df["type"] == "nontrna") & (df["Sample"] == "s2")]
    assert nontrna_s2["share"].tolist() == pytest.approx([0.5, 0.5])


def test_histogram_tolerates_a_type_with_no_reads():
    adata = _build_adata()
    adata.uns["mismatch_counts"] = pd.DataFrame({
        "count": [0, 1],
        "type": ["nontrna", "nontrna"],
        "s1": [0, 0],
    })

    df = plotsMismatch.mismatch_count_shares(adata)

    assert df["share"].tolist() == [0.0, 0.0]


# --- rendering ----------------------------------------------------------------------

def test_per_amino_plots_render_when_anticodon_is_categorical(tmp_path):
    """seaborn reads hue levels from a categorical's full category list, not the observed
    values, so a per-amino subset demanded a palette covering every anticodon in the
    dataset and raised `The palette dictionary is missing keys`."""
    import matplotlib
    matplotlib.use("Agg")

    adata = _build_adata(trnas=("tRNA-Ala-AGC-1", "tRNA-Gly-GCC-1"))
    adata.obs["amino"] = adata.obs["amino"].astype("category")
    adata.obs["iso"] = adata.obs["iso"].astype("category")
    for label in adata.obs_names:
        _set(adata, label, "coverage", "1", 90)
        _set(adata, label, "mismatchedbases", "1", 30)
        _set(adata, label, "deletedbases", "1", 5)
        _passing_bases(adata, label, "1")

    output = str(tmp_path) + "/"
    plotsMismatch.visualizer(adata, {}, output, 10, threaded=False).generate_plots()

    assert (tmp_path / "positionmismatches.pdf").exists()
    assert (tmp_path / "positiondeletions.pdf").exists()
    assert (tmp_path / "individual" / "Ala_positionmismatches.pdf").exists()
    assert (tmp_path / "individual" / "Gly_positionmismatches.pdf").exists()
    # The multi-page roll-up sits at the base next to the overview, not in individual/ --
    # plotsCoverage's convention, and it stops the roll-up reading as a copy of the overview.
    assert (tmp_path / "combined_positionmismatches_by_amino.pdf").exists()
    assert not (tmp_path / "individual" / "positionmismatches_combined.pdf").exists()


def test_histogram_renders_and_a_user_colormap_is_honoured(tmp_path):
    import matplotlib
    matplotlib.use("Agg")

    adata = _build_adata()
    for label in adata.obs_names:
        _set(adata, label, "coverage", "1", 90)
        _set(adata, label, "mismatchedbases", "1", 30)
        _passing_bases(adata, label, "1")
    adata.uns["mismatch_counts"] = pd.DataFrame({
        "count": [0, 1], "type": ["trna", "trna"], "s1": [75, 25], "s2": [60, 40],
    })

    colormap = {"amino": {"Ala": "#FFAE49", "Gly": "#44B7C2"}, "mismatch": {"0": "#024B7A"}}
    output = str(tmp_path) + "/"
    plotsMismatch.visualizer(adata, colormap, output, 10, threaded=False).generate_plots()

    assert (tmp_path / "mismatchcounts.pdf").exists()


def test_combined_roll_up_has_one_page_per_amino(tmp_path):
    """It looked like 'a copy with only Ala' -- pin the page count so that stays checkable."""
    import matplotlib
    matplotlib.use("Agg")

    trnas = ("tRNA-Ala-AGC-1", "tRNA-Gly-GCC-1", "tRNA-Val-CAC-1")
    adata = _build_adata(trnas=trnas)
    adata.obs["amino"] = adata.obs["amino"].astype("category")
    adata.obs["iso"] = adata.obs["iso"].astype("category")
    for label in adata.obs_names:
        _set(adata, label, "coverage", "1", 90)
        _set(adata, label, "mismatchedbases", "1", 30)
        _passing_bases(adata, label, "1")

    plotsMismatch.visualizer(adata, {}, str(tmp_path) + "/", 10, threaded=False).generate_plots()

    combined = (tmp_path / "combined_positionmismatches_by_amino.pdf").read_bytes()
    pages = combined.count(b"/Type /Page\n") or combined.count(b"/Type /Page")
    assert pages >= len(trnas), f"expected one page per amino acid, found {pages}"
    for amino in ("Ala", "Gly", "Val"):
        assert (tmp_path / "individual" / f"{amino}_positionmismatches.pdf").exists()
