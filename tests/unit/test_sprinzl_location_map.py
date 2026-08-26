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
    loc, half, _code = _sprinzl_location_maps()
    assert all(isinstance(k, str) for k in loc), sorted(k for k in loc if not isinstance(k, str))
    assert all(isinstance(k, str) for k in half), sorted(k for k in half if not isinstance(k, str))


@pytest.mark.parametrize("position,expected", sorted(ANTICODON_STEM.items()))
def test_anticodon_stem_positions_are_labeled(position, expected):
    loc, _half, _code = _sprinzl_location_maps()
    assert loc.get(position) == expected


@pytest.mark.parametrize("position", sorted(ANTICODON_STEM))
def test_anticodon_stem_counts_as_the_centre_half(position):
    """`half` splits a tRNA into 5'/centre/3' for fragment typing. The anticodon LOOP was
    already 'center'; the stem flanking it silently was not."""
    _loc, half, _code = _sprinzl_location_maps()
    assert half.get(position) == "center"


@pytest.mark.parametrize("orgmode", sorted(toolsGetCoverage.POSITION_TABLES))
def test_every_sprinzl_position_of_every_organism_mode_is_mapped(orgmode):
    """No position that the coverage step can actually emit may be left without a region. This
    is the assertion that would have caught the original defect for any arm, not just this one."""
    loc, half, _code = _sprinzl_location_maps()
    unmapped_loc = [str(p) for p in toolsGetCoverage.POSITION_TABLES[orgmode] if str(p) not in loc]
    unmapped_half = [str(p) for p in toolsGetCoverage.POSITION_TABLES[orgmode] if str(p) not in half]
    assert unmapped_loc == [], f"{orgmode}: positions with no location: {unmapped_loc}"
    assert unmapped_half == [], f"{orgmode}: positions with no half: {unmapped_half}"


# ---------------------------------------------------------------------------
# Region vocabulary, ported from tRNAscan-SE's SprinzlPos.pm (Patricia Chan, UCSC).
# Its `ss_pos` table is the reference: the acceptor stem is 1-7 paired with 66-72
# ("A1:U72 base pair at the end of the acceptor stem", tRNAscan-SE 2.0, NAR 2021),
# 73-76 are their own region L7 "3p end", 44-48 are L5 "Variable Loop", and the
# e-series positions are P4 "Variable Stem" -- the reverse of what tRNAgraph used
# to call them.
# ---------------------------------------------------------------------------

EXPECTED_REGIONS = {
    "fiveprime_extra": (["-1"], None),
    "fiveprime_acceptorstem": ([str(i) for i in range(1, 8)], "5P1"),
    "a_to_d_internal": (["8", "9"], "L1"),
    "fiveprime_dstem": ([str(i) for i in range(10, 14)], "5P2"),
    "dloop": ([str(i) for i in range(14, 22)] + ["17a", "20a", "20b"], "L2"),
    "threeprime_dstem": ([str(i) for i in range(22, 26)], "3P2"),
    "d_to_anticodon_internal": (["26"], "L3"),
    "fiveprime_anticodonstem": ([str(i) for i in range(27, 32)], "5P3"),
    "anticodonloop": ([str(i) for i in range(32, 39)], "L4"),
    "threeprime_anticodonstem": ([str(i) for i in range(39, 44)], "3P3"),
    "variableloop": ([str(i) for i in range(44, 49)], "L5"),
    "variablestem": (["e%d" % i for i in range(1, 20)], "P4"),
    "fiveprime_tstem": ([str(i) for i in range(49, 54)], "5P5"),
    "tloop": ([str(i) for i in range(54, 61)], "L6"),
    "threeprime_tstem": ([str(i) for i in range(61, 66)], "3P5"),
    "threeprime_acceptorstem": ([str(i) for i in range(66, 73)], "3P1"),
    "threeprime_end": ([str(i) for i in range(73, 77)], "L7"),
}


@pytest.mark.parametrize("region", sorted(EXPECTED_REGIONS))
def test_region_membership_matches_sprinzlpos(region):
    loc, _half, _code = _sprinzl_location_maps()
    positions, _short = EXPECTED_REGIONS[region]
    assert sorted(p for p, r in loc.items() if r == region) == sorted(positions)


def test_acceptor_stem_is_symmetric():
    """The defect this vocabulary change exists to fix: the 5' side labeled 8 positions
    (-1 through 7) and the 3' side 11 (66-76), folding the discriminator base and the CCA
    tail into a seven-base-pair stem."""
    loc, _half, _code = _sprinzl_location_maps()
    five = [p for p, r in loc.items() if r == "fiveprime_acceptorstem"]
    three = [p for p, r in loc.items() if r == "threeprime_acceptorstem"]
    assert len(five) == len(three) == 7


def test_phantom_zero_position_is_gone():
    """Sprinzl numbering jumps -1 -> 1; '0' was generated by range(-1, 8) and matched nothing."""
    loc, half, code = _sprinzl_location_maps()
    assert "0" not in loc and "0" not in half and "0" not in code


@pytest.mark.parametrize("region", sorted(EXPECTED_REGIONS))
def test_region_short_codes_match_sprinzlpos(region):
    _loc, _half, code = _sprinzl_location_maps()
    positions, short = EXPECTED_REGIONS[region]
    assert {code.get(p) for p in positions} == {short}


def test_minus_one_is_the_only_position_without_a_canonical_code():
    """SprinzlPos.pm's position list starts at 1, so the 5' extra base (the G-1 of tRNA-His)
    has no region code in the canonical scheme. Inventing one would be worse than leaving it
    unset, so it is the single documented exception."""
    loc, _half, code = _sprinzl_location_maps()
    assert sorted(p for p in loc if code.get(p) is None) == ["-1"]


def test_fragment_heuristic_regions_preserve_pre_rename_membership():
    """obs['fragment'] is classified from mean coverage over the acceptor-stem regions. The
    3' side used to include 73-76, so the CCA tail -- where 3' fragments pile up -- counted
    toward that mean. Splitting 73-76 into their own region must not silently change fragment
    calls, so the heuristic names both regions explicitly rather than tracking the rename."""
    from trnagraph.modules.adataBuild import (
        FRAGMENT_FIVEPRIME_REGIONS, FRAGMENT_THREEPRIME_REGIONS, _sprinzl_location_maps,
    )
    loc, _half, _code = _sprinzl_location_maps()
    five = {p for p, r in loc.items() if r in FRAGMENT_FIVEPRIME_REGIONS}
    three = {p for p, r in loc.items() if r in FRAGMENT_THREEPRIME_REGIONS}
    assert five == {str(i) for i in range(-1, 8)} - {"0"}
    assert three == {str(i) for i in range(66, 77)}
