"""Regression test for toolsTG.parse_variant() actually constructing a real VariantTag on its
success path. Added after a real bug: an edit to toolsTG.py accidentally deleted the module's
`from .toolsSchemas import VariantTag` import, and nothing caught it -- every existing test
exercised parse_variant()'s validation-error branches (which never reach the VariantTag(...)
call), so a plain `trnagraph graph` invocation was the first thing to hit the resulting
NameError. This locks in the success path so a future editing slip fails fast in CI/unit tests
instead of only at the CLI."""
import anndata as ad
import numpy as np
import pandas as pd

from trnagraph.modules.toolsTG import parse_variant
from trnagraph.modules.toolsSchemas import VariantTag


def _make_adata():
    obs = pd.DataFrame({"sample": ["s0", "s1"]}, index=["obs0", "obs1"])
    return ad.AnnData(X=np.zeros((2, 1)), obs=obs)


def test_parse_variant_constructs_a_real_varianttag_for_the_default():
    result = parse_variant(_make_adata(), "norm:full")

    assert isinstance(result, VariantTag)
    assert result.norm == "norm"
    assert result.tag == "full"


def test_parse_variant_constructs_a_real_varianttag_for_an_existing_split():
    adata = _make_adata()
    adata.uns["size_splits"] = {"u60": {}}

    result = parse_variant(adata, "raw:u60")

    assert isinstance(result, VariantTag)
    assert result.norm == "raw"
    assert result.tag == "u60"
