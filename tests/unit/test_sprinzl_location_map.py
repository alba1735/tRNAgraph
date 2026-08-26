"""Regression tests for the Sprinzl position -> structural region annotation (`adata.var`).

`_adata_build_` built its position->region dictionaries inline, so nothing could check them
without a full AnnData build. They were wrong: the anticodon-stem entries were keyed by `int`
(`range(27,32)` / `range(39,44)`) while every other entry used `str(i)` and `var['positions']`
holds strings, so those ten positions matched nothing and came out NaN. Measured on the
vibrChol1 build, 150 of 1395 var rows carried a null `location` AND a null `half` -- the
anticodon stem, one of the four canonical arms, was unlabeled in every object tRNAgraph
had ever written. The anticodon LOOP was unaffected, because it used `str(i)`.
"""
import pytest

from trnagraph.modules import toolsGetCoverage
from trnagraph.modules.adataBuild import _sprinzl_location_maps


ANTICODON_STEM = {
    "27": "fiveprime_anticodonstem", "28": "fiveprime_anticodonstem",
    "29": "fiveprime_anticodonstem", "30": "fiveprime_anticodonstem",
    "31": "fiveprime_anticodonstem",
    "39": "threeprime_anticodonstem", "40": "threeprime_anticodonstem",
    "41": "threeprime_anticodonstem", "42": "threeprime_anticodonstem",
    "43": "threeprime_anticodonstem",
}


def test_every_key_is_a_string():
    """The defect in one assertion: `var['positions']` holds strings, so an int key can never
    match. Nothing may be keyed by int, regardless of which region it belongs to."""
    loc, half = _sprinzl_location_maps()
    assert all(isinstance(k, str) for k in loc), sorted(k for k in loc if not isinstance(k, str))
    assert all(isinstance(k, str) for k in half), sorted(k for k in half if not isinstance(k, str))


@pytest.mark.parametrize("position,expected", sorted(ANTICODON_STEM.items()))
def test_anticodon_stem_positions_are_labeled(position, expected):
    loc, _half = _sprinzl_location_maps()
    assert loc.get(position) == expected


@pytest.mark.parametrize("position", sorted(ANTICODON_STEM))
def test_anticodon_stem_counts_as_the_centre_half(position):
    """`half` splits a tRNA into 5'/centre/3' for fragment typing. The anticodon LOOP was
    already 'center'; the stem flanking it silently was not."""
    _loc, half = _sprinzl_location_maps()
    assert half.get(position) == "center"


@pytest.mark.parametrize("orgmode", sorted(toolsGetCoverage.POSITION_TABLES))
def test_every_sprinzl_position_of_every_organism_mode_is_mapped(orgmode):
    """No position that the coverage step can actually emit may be left without a region. This
    is the assertion that would have caught the original defect for any arm, not just this one."""
    loc, half = _sprinzl_location_maps()
    unmapped_loc = [str(p) for p in toolsGetCoverage.POSITION_TABLES[orgmode] if str(p) not in loc]
    unmapped_half = [str(p) for p in toolsGetCoverage.POSITION_TABLES[orgmode] if str(p) not in half]
    assert unmapped_loc == [], f"{orgmode}: positions with no location: {unmapped_loc}"
    assert unmapped_half == [], f"{orgmode}: positions with no half: {unmapped_half}"
