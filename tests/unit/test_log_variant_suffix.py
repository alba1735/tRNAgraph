"""Log filenames must name the variant a run was for.

`graph`, `cluster` and `tools log2fc` all take `--variant`, and running one of them in a
loop over `norm:full`, `norm:u60`, `norm:o60` produced three logs whose names differed
only by timestamp -- so telling them apart meant opening each one and reading its contents.
The variant is part of what a run *was*, so it belongs in the name alongside the object.
"""
from trnagraph.cli import _log_suffix


def test_default_variant_is_omitted():
    """norm:full is the default, so naming it would add noise to every ordinary run."""
    assert _log_suffix("results/trnagraph.h5ad", "norm:full") == "trnagraph"


def test_split_variant_is_named():
    assert _log_suffix("results/trnagraph.h5ad", "norm:u60") == "trnagraph_norm_u60"


def test_normalization_variant_is_named():
    assert _log_suffix("results/trnagraph.h5ad", "vst:full") == "trnagraph_vst_full"


def test_every_variant_of_one_object_gets_a_distinct_name():
    """The actual complaint: a loop over variants produced indistinguishable log names."""
    variants = ["norm:full", "norm:u60", "norm:o60", "vst:u60", "allfeatures:full"]
    names = [_log_suffix("trnagraph.h5ad", v) for v in variants]

    assert len(set(names)) == len(variants)


def test_no_variant_falls_back_to_the_object_name():
    """Commands without a --variant flag keep the name they already had."""
    assert _log_suffix("results/trnagraph.h5ad") == "trnagraph"
    assert _log_suffix("results/trnagraph.h5ad", None) == "trnagraph"


def test_the_colon_is_not_carried_into_a_filename():
    suffix = _log_suffix("trnagraph.h5ad", "norm:u60")

    assert ":" not in suffix
    assert "/" not in suffix
