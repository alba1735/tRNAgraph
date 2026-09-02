"""Declared category order becomes an ordered Categorical on adata.obs.

One source of truth, two consumers. Plot legends and axes read the category order; DESeq2
reads the FIRST category as the reference level every contrast is taken against. Storing it on
the object rather than passing it per call is what keeps those two from disagreeing -- a
transient flag would have to be re-supplied identically on every plot and DE invocation, and
would not reach the statistics at all.

The error paths matter more than the happy path here. A level present in the data but missing
from the declaration cannot be silently appended: the declaration would then no longer describe
the order, and the figure would be wrong in a way nothing announces. A column named in the
declaration but absent from obs is almost always a typo, and doing nothing about it leaves the
user staring at an unchanged figure.
"""
import pandas as pd
import pytest

from trnagraph.modules import toolsTG


def _obs():
    return pd.DataFrame({
        'timepoint': ['Day 35', 'Day 0', 'Day 70', 'Day 0'],
        'sample': ['a', 'b', 'c', 'd'],
    })


def test_declared_order_becomes_an_ordered_categorical():
    obs = toolsTG.apply_category_order(_obs(), {'timepoint': ['Day 0', 'Day 35', 'Day 70']})

    assert isinstance(obs['timepoint'].dtype, pd.CategoricalDtype)
    assert obs['timepoint'].cat.ordered
    assert list(obs['timepoint'].cat.categories) == ['Day 0', 'Day 35', 'Day 70']


def test_declared_order_survives_where_alphabetical_would_be_wrong():
    """D0/D7/D14/D35 sorts alphabetically to D0, D14, D35, D7 -- chronology and alphabet
    genuinely disagree, which is the whole reason this key exists. (D0/D14/D35/D70 would NOT
    demonstrate it: there the two happen to coincide.)"""
    obs = pd.DataFrame({'timepoint': ['D35', 'D0', 'D14', 'D7']})
    declared = ['D0', 'D7', 'D14', 'D35']

    obs = toolsTG.apply_category_order(obs, {'timepoint': declared})

    assert list(obs['timepoint'].cat.categories) == declared
    assert list(obs['timepoint'].cat.categories) != sorted(declared)
    assert obs['timepoint'].cat.categories[0] == 'D0', 'the DESeq2 reference level'


def test_a_level_in_the_data_but_not_declared_is_an_error_naming_both():
    """Appending it silently would make the declaration stop describing the order."""
    with pytest.raises(Exception) as excinfo:
        toolsTG.apply_category_order(_obs(), {'timepoint': ['Day 0', 'Day 35']})

    message = str(excinfo.value)
    assert 'Day 70' in message, 'names the undeclared level'
    assert 'Day 0' in message and 'Day 35' in message, 'names what WAS declared'


def test_a_column_that_does_not_exist_is_an_error_naming_it():
    """Almost always a typo; doing nothing leaves the user staring at an unchanged figure."""
    with pytest.raises(Exception) as excinfo:
        toolsTG.apply_category_order(_obs(), {'timepont': ['Day 0']})

    message = str(excinfo.value)
    assert 'timepont' in message
    assert 'timepoint' in message, 'lists the available columns so the typo is visible'


def test_a_declared_level_absent_from_the_data_is_kept():
    """A filtered subset legitimately lacks a timepoint; dropping it would make the order
    depend on the subset rather than on the declaration."""
    obs = toolsTG.apply_category_order(_obs(), {'timepoint': ['Day 0', 'Day 14', 'Day 35', 'Day 70']})

    assert 'Day 14' in list(obs['timepoint'].cat.categories)
