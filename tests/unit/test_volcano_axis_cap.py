"""The volcano x-axis was scaled to the single most extreme log2FC.

`ax.set_xlim(-1.1 * max(abs(x)), 1.1 * max(abs(x)))` let one outlier set the axis for every
gene. Measured on the hg38 object, 27 of 36 comparison/readtype combinations had the bulk of
the data (median |LFC|) occupying under 20% of the half-axis -- worst case 11.6%, where the
median was 1.12 against a max of 8.78.

The fix caps the axis at a percentile of |LFC| and pins out-of-range points to the boundary
as triangles, so nothing is dropped and the reader can see that something lies beyond. The
cap only engages when there is a real tail: a dataset whose max is close to its percentile
keeps the old max-based axis rather than gaining a pointless edge marker.

Shrinkage is not an alternative to this. Measured on the same data, apeGLM moved the median
|LFC| from 1.00 to 0.62 but the max from 8.28 to 9.02 -- the extremes are real (69 of the 82
features above |LFC| 5 are significant at padj<0.05), which is exactly what a heavy-tailed
prior is designed to preserve.
"""
import numpy as np
import pandas as pd
import pytest

from trnagraph.modules.plotsVolcano import (
    VOLCANO_CAP_ENGAGE_RATIO,
    VOLCANO_CAP_PERCENTILE,
    VOLCANO_MIN_HALF_WIDTH,
    resolve_x_limit,
)


def test_default_percentile_is_the_measured_choice():
    """p99 barely moved the hg38 axis; p95 more than doubled the bulk's share of it."""
    assert VOLCANO_CAP_PERCENTILE == 95


def test_a_long_tail_engages_the_cap():
    values = pd.Series(list(np.full(95, 1.0)) + list(np.linspace(5, 9, 5)))

    cap, capped = resolve_x_limit(values)

    assert capped
    assert cap < values.abs().max()


def test_a_tame_distribution_keeps_the_max_based_axis():
    """vibrChol1's worst comparison: max 3.13 against a p95 near 2 -- no tail to tame."""
    values = pd.Series(list(np.full(50, 1.0)) + [2.0, 2.5, 3.13])

    cap, capped = resolve_x_limit(values)

    assert not capped
    assert cap == pytest.approx(values.abs().max())


def test_the_cap_never_falls_below_the_minimum_half_width():
    """Preserves the old behaviour of squaring up small plots to +-3."""
    values = pd.Series(np.full(200, 0.2))

    cap, _ = resolve_x_limit(values)

    assert cap >= VOLCANO_MIN_HALF_WIDTH


def test_explicit_limit_overrides_everything():
    values = pd.Series(list(np.full(95, 1.0)) + list(np.linspace(5, 9, 5)))

    cap, capped = resolve_x_limit(values, explicit=2.5)

    assert cap == 2.5
    assert capped, "points beyond an explicit limit are still off-scale"


def test_explicit_limit_wider_than_the_data_does_not_claim_to_cap():
    values = pd.Series([0.5, 1.0, 1.5])

    cap, capped = resolve_x_limit(values, explicit=10)

    assert cap == 10
    assert not capped


def test_engage_ratio_is_the_documented_threshold():
    """A max just under ratio x percentile must not engage; just over must."""
    base = list(np.full(99, 1.0))
    p95 = pd.Series(base + [4.0]).abs().quantile(0.95)
    just_under = pd.Series(base + [max(VOLCANO_MIN_HALF_WIDTH, p95) * (VOLCANO_CAP_ENGAGE_RATIO - 0.2)])
    just_over = pd.Series(base + [max(VOLCANO_MIN_HALF_WIDTH, p95) * (VOLCANO_CAP_ENGAGE_RATIO + 0.5)])

    assert not resolve_x_limit(just_under)[1]
    assert resolve_x_limit(just_over)[1]


def test_negative_values_are_measured_by_magnitude():
    positive = pd.Series(list(np.full(95, 1.0)) + list(np.linspace(5, 9, 5)))
    negative = -positive

    assert resolve_x_limit(positive) == resolve_x_limit(negative)


def test_empty_input_is_survivable():
    cap, capped = resolve_x_limit(pd.Series(dtype=float))

    assert cap >= VOLCANO_MIN_HALF_WIDTH
    assert not capped


def test_all_nan_input_is_survivable():
    cap, capped = resolve_x_limit(pd.Series([np.nan, np.nan]))

    assert cap >= VOLCANO_MIN_HALF_WIDTH
    assert not capped


def test_zero_or_negative_explicit_limit_is_rejected():
    with pytest.raises(ValueError):
        resolve_x_limit(pd.Series([1.0, 2.0]), explicit=0)
    with pytest.raises(ValueError):
        resolve_x_limit(pd.Series([1.0, 2.0]), explicit=-3)


def test_hg38_worst_case_more_than_doubles_the_bulks_share_of_the_axis():
    """The regression this exists for, using the real measured shape."""
    rng = np.random.default_rng(0)
    bulk = np.abs(rng.normal(1.12, 0.6, 100))
    tail = np.linspace(4.6, 8.78, 5)
    values = pd.Series(np.concatenate([bulk, tail]))

    cap, capped = resolve_x_limit(values)
    assert capped

    before = values.median() / (1.1 * values.abs().max())
    after = values.median() / cap
    assert after > 2 * before
