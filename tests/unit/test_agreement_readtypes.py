"""Which read types an agreement figure carries, and which level it is anchored on.

`total` is excluded BY CONSTRUCTION rather than by documentation: it is the sum of the other
four, so plotting it beside them draws the same reads twice and a tRNA appears to agree with
itself. The heatmap's total-vs-combine naming documents the same trap.
"""
import pandas as pd
import pytest

from trnagraph.modules import plotsAgreement
from trnagraph.modules.toolsSchemas import MultivariateConfig


def test_total_is_dropped_even_when_asked_for():
    kept = plotsAgreement.agreement_readtypes(
        ['nreads_fiveprime_unique_norm', 'nreads_total_unique_norm',
         'nreads_wholecounts_unique_norm'])

    assert kept == ['nreads_fiveprime_unique_norm', 'nreads_wholecounts_unique_norm']


def test_dropping_total_never_empties_the_figure():
    """Asking for total alone is a request that cannot be honoured, and saying so beats
    rendering an empty axes."""
    with pytest.raises(plotsAgreement.InvalidAgreementError) as excinfo:
        plotsAgreement.agreement_readtypes(['nreads_total_unique_norm'])

    assert 'total' in str(excinfo.value)


def test_the_bare_readtype_name_is_matched_not_the_whole_column():
    """Callers hold resolved obs columns, which carry unique/norm suffixes."""
    kept = plotsAgreement.agreement_readtypes(
        ['nreads_total_norm', 'nreads_threeprime_norm'])

    assert kept == ['nreads_threeprime_norm']


def test_the_reference_defaults_to_the_first_level():
    levels = pd.Categorical(['Day 70', 'Day 0', 'Day 35'],
                            categories=['Day 0', 'Day 35', 'Day 70'], ordered=True)

    assert plotsAgreement.resolve_reference(levels, MultivariateConfig()) == 'Day 0'


def test_an_ordered_category_decides_the_first_level_not_the_data_order():
    """`order` in --config exists so a timecourse reads in experimental order; the reference
    follows it, which is also the DESeq2 reference level."""
    levels = pd.Categorical(['Day 70', 'Day 0'], categories=['Day 70', 'Day 0'], ordered=True)

    assert plotsAgreement.resolve_reference(levels, MultivariateConfig()) == 'Day 70'


def test_an_explicit_reference_wins():
    levels = pd.Categorical(['Day 0', 'Day 70'], categories=['Day 0', 'Day 70'], ordered=True)

    assert plotsAgreement.resolve_reference(
        levels, MultivariateConfig(reference='Day 70')) == 'Day 70'


def test_a_reference_that_is_not_a_level_is_refused_by_name():
    levels = pd.Categorical(['Day 0', 'Day 70'], categories=['Day 0', 'Day 70'], ordered=True)

    with pytest.raises(plotsAgreement.InvalidAgreementError) as excinfo:
        plotsAgreement.resolve_reference(levels, MultivariateConfig(reference='Day 99'))

    assert 'Day 99' in str(excinfo.value)
