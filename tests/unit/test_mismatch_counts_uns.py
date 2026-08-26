"""Tests for `adata.uns['mismatch_counts']`, the read-level mismatch histogram.

The per-position mismatch data has always been in the object (the `mismatchedbases`,
`deletedbases`, `adenines`/`thymines`/`cytosines`/`guanines` coverage types in `var`),
but the read-level histogram -- how many reads carried 0, 1, 2 ... mismatches -- existed
only as the flat `<exp>-mismatches.txt`. `graph -g mismatch` renders it, and plots in this
project are built from the AnnData object rather than from result files, so it has to
live in `uns`.

The on-disk file is size-factor-normalized (`toolsCountReads.printmismatchcounts` divides
by the sample's size factor). The object stores raw counts instead, recovered by undoing
that division.
"""
import pytest

from trnagraph.modules.adataBuild import _load_mismatch_counts


def _write(path, header, rows):
    path.write_text("\n".join(["\t".join(header)] + ["\t".join(map(str, r)) for r in rows]) + "\n")


def test_raw_counts_are_recovered_from_the_normalized_file(tmp_path):
    size_factors = {"s1": 2.0, "s2": 0.5}
    raw = {"s1": [400, 300], "s2": [50, 25]}

    path = tmp_path / "exp-mismatches.txt"
    _write(path, ["count", "type", "s1", "s2"], [
        [0, "trna", raw["s1"][0] / size_factors["s1"], raw["s2"][0] / size_factors["s2"]],
        [1, "trna", raw["s1"][1] / size_factors["s1"], raw["s2"][1] / size_factors["s2"]],
    ])

    df = _load_mismatch_counts(str(path), size_factors)

    assert list(df.columns) == ["count", "type", "s1", "s2"]
    assert df["s1"].tolist() == raw["s1"]
    assert df["s2"].tolist() == raw["s2"]


def test_recovered_counts_are_integers(tmp_path):
    """Real magnitudes: vibrChol1's largest bin against its own primary size factor."""
    size_factor = 0.9103082380957978
    raw = 60077

    path = tmp_path / "exp-mismatches.txt"
    _write(path, ["count", "type", "s1"], [[0, "trna", repr(raw / size_factor)]])

    df = _load_mismatch_counts(str(path), {"s1": size_factor})

    assert df["s1"].dtype.kind in "iu"
    assert df["s1"].tolist() == [raw]


def test_count_and_type_columns_are_preserved(tmp_path):
    path = tmp_path / "exp-mismatches.txt"
    _write(path, ["count", "type", "s1"], [
        [0, "trna", 10.0], [0, "nontrna", 4.0], [1, "trna", 6.0], [1, "nontrna", 2.0],
    ])

    df = _load_mismatch_counts(str(path), {"s1": 1.0})

    assert df["count"].tolist() == [0, 0, 1, 1]
    assert df["type"].tolist() == ["trna", "nontrna", "trna", "nontrna"]


def test_missing_size_factor_leaves_the_column_unscaled(tmp_path):
    """A sample absent from the size-factor set is a build inconsistency, not a crash."""
    path = tmp_path / "exp-mismatches.txt"
    _write(path, ["count", "type", "s1", "s2"], [[0, "trna", 10.0, 8.0]])

    df = _load_mismatch_counts(str(path), {"s1": 2.0})

    assert df["s1"].tolist() == [20]
    assert df["s2"].tolist() == [8]


def test_absent_file_returns_none(tmp_path):
    assert _load_mismatch_counts(str(tmp_path / "nope.txt"), {"s1": 1.0}) is None


def test_percentages_are_unchanged_by_the_rescale(tmp_path):
    """The histogram plot renders per-sample percentages, which the size factor cannot move."""
    path = tmp_path / "exp-mismatches.txt"
    _write(path, ["count", "type", "s1"], [[0, "trna", 30.0], [1, "trna", 10.0]])

    df = _load_mismatch_counts(str(path), {"s1": 3.0})

    share = df["s1"] / df["s1"].sum()
    assert share.tolist() == pytest.approx([0.75, 0.25])


# --- split variants ----------------------------------------------------------------

def test_split_variants_drop_the_nontrna_series():
    """Read-length splits exclude non-tRNA features, so the histogram must too."""
    import pandas as pd
    from trnagraph.modules.adataBuild import _trna_only_mismatch_counts

    df = pd.DataFrame({
        "count": [0, 0, 1, 1],
        "type": ["trna", "nontrna", "trna", "nontrna"],
        "s1": [10, 4, 6, 2],
    })

    out = _trna_only_mismatch_counts(df)

    assert out["type"].tolist() == ["trna", "trna"]
    assert out["s1"].tolist() == [10, 6]
    assert out.index.tolist() == [0, 1]


def test_split_variant_helper_passes_none_through():
    from trnagraph.modules.adataBuild import _trna_only_mismatch_counts

    assert _trna_only_mismatch_counts(None) is None
