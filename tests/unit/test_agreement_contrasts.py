"""Which contrasts an agreement figure is built from.

The universe is reference-anchored, not all-pairwise: with a reference and three other levels
that is three contrasts, not the six pairwise combinations. The reference figure's legend names
exactly this shape -- 0-14, 0-35 and 0-70 for a four-timepoint column.

One contrast is drawn on the x-axis; colour encodes how many of the whole universe a feature is
significant in. So the universe has to be derivable before anything is plotted, and it has to
stay stable regardless of which contrast is being drawn.
"""
import pytest

from trnagraph.modules import plotsAgreement


def test_every_contrast_is_anchored_on_the_reference():
    universe = plotsAgreement.contrast_universe(['Day 0', 'Day 14', 'Day 35', 'Day 70'], 'Day 0')

    assert universe == [('Day 0', 'Day 14'), ('Day 0', 'Day 35'), ('Day 0', 'Day 70')]


def test_the_reference_is_not_contrasted_with_itself():
    universe = plotsAgreement.contrast_universe(['a', 'b'], 'a')

    assert universe == [('a', 'b')]


def test_level_order_is_preserved():
    """An ordered category exists so figures read in experimental order, not alphabetical."""
    universe = plotsAgreement.contrast_universe(['Day 0', 'Day 70', 'Day 35'], 'Day 0')

    assert universe == [('Day 0', 'Day 70'), ('Day 0', 'Day 35')]


def test_a_reference_absent_from_the_levels_is_an_error():
    """Naming a level that is not there produces an empty universe and a blank figure, which
    looks like a result. The message names what was asked for and what exists."""
    with pytest.raises(plotsAgreement.InvalidAgreementError) as excinfo:
        plotsAgreement.contrast_universe(['Day 0', 'Day 70'], 'Day 99')

    assert 'Day 99' in str(excinfo.value)
    assert 'Day 0' in str(excinfo.value)


def test_a_single_level_cannot_be_contrasted():
    with pytest.raises(plotsAgreement.InvalidAgreementError):
        plotsAgreement.contrast_universe(['Day 0'], 'Day 0')
