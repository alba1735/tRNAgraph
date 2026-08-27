"""Regression tests for toolsSchemas.py's pydantic validation of tRNAgraph's four user-facing
input files (roadmap.md Phase 2: "Pydantic validation across input files")."""
import pytest
from pydantic import ValidationError

from trnagraph.modules.toolsSchemas import GraphFilterConfig, StyleFile, MetadataFile, PairsFile


def test_graph_filter_config_accepts_minimal_valid_dict():
    cfg = GraphFilterConfig.model_validate({"name": "analysis_subset"})
    assert cfg.name == "analysis_subset"
    assert cfg.obs is None
    assert cfg.var is None


def test_graph_filter_config_requires_name():
    with pytest.raises(ValidationError):
        GraphFilterConfig.model_validate({"obs": {"treatment": ["treatment_A"]}})


def test_graph_filter_config_rejects_path_traversal_in_name():
    with pytest.raises(ValidationError):
        GraphFilterConfig.model_validate({"name": "../../etc"})


def test_graph_filter_config_allows_var_r_without_var():
    """
    Previously a `var_r`-without-`var` config crashed downstream as `KeyError: 'var'` since the
    code did `filter_dict = self.args.config['var']` unconditionally whenever 'var_r' was
    present. The model itself must accept this combination -- 'var' legitimately stays None,
    and it's on the caller to build the reverse filter from adata rather than from `.var`.
    """
    cfg = GraphFilterConfig.model_validate({"name": "x", "var_r": {"coverage": ["unique"]}})
    assert cfg.var is None
    assert cfg.var_r == {"coverage": ["unique"]}


def test_graph_filter_config_rejects_unknown_keys():
    with pytest.raises(ValidationError):
        GraphFilterConfig.model_validate({"name": "x", "typo_key": {"a": [1]}})


def test_style_file_accepts_a_legacy_colormap_shaped_file():
    """ColormapFile was folded into StyleFile; the old file shape must still load."""
    style = StyleFile.model_validate({"group": {"VC_24h": "#FFAE49", "VC_log": "#44B7C2"}})
    assert style.colors_for("group")["VC_24h"] == "#FFAE49"


def test_style_file_rejects_non_string_color_value():
    with pytest.raises(ValidationError):
        StyleFile.model_validate({"colors": {"group": {"VC_24h": 12345}}})


def test_metadata_file_accepts_valid_rows():
    meta = MetadataFile(
        path="metadata.tsv",
        header=["sample", "group"],
        rows=[["s1", "A"], ["s2", "B"]],
    )
    assert len(meta.rows) == 2


def test_metadata_file_rejects_duplicate_sample_names():
    with pytest.raises(ValidationError):
        MetadataFile(
            path="metadata.tsv",
            header=["sample", "group"],
            rows=[["s1", "A"], ["s1", "B"]],
        )


def test_metadata_file_rejects_row_length_mismatch():
    with pytest.raises(ValidationError):
        MetadataFile(
            path="metadata.tsv",
            header=["sample", "group", "timepoint"],
            rows=[["s1", "A", "0h"], ["s2", "B"]],
        )


def test_metadata_file_rejects_empty_file():
    with pytest.raises(ValidationError):
        MetadataFile(path="metadata.tsv", header=["sample", "group"], rows=[])


def test_pairs_file_accepts_valid_pairs():
    pairs = PairsFile(path="pairs.txt", pairs=[("s1", "s2"), ("s3", "s4")])
    assert len(pairs.pairs) == 2


def test_pairs_file_rejects_empty_pairs():
    with pytest.raises(ValidationError):
        PairsFile(path="pairs.txt", pairs=[])
