"""One shape per read type, shared by every plot that mixes them.

The combined heatmap stacks several read types under the SAME tRNA names, so a tRNA can
appear four times with nothing to say which row is which. A marker beside the label
disambiguates them.

The vocabulary is not new: it is the one already used in the project's published volcano
figure -- circle 5', square 3', triangle other, diamond whole, cross pre-tRNA -- and an
abandoned block in plotsHeatmap tried to draw exactly these before. That attempt hardcoded
obs column names and got them wrong (`nreads_whole_unique_norm`, which does not exist; the
real column is `nreads_wholecounts_*`), so whole-count rows would have gone unmarked. The
lookup therefore normalizes a full obs column name down to its bare read type rather than
matching literal strings.
"""
import pytest

from trnagraph.modules import plotsPalette


@pytest.mark.parametrize('readtype,glyph,label', [
    ('fiveprime', '●', "5' counts"),
    ('threeprime', '■', "3' counts"),
    ('other', '▲', 'Other counts'),
    ('wholecounts', '◆', 'Whole counts'),
])
def test_each_read_type_has_its_published_glyph(readtype, glyph, label):
    marker = plotsPalette.readtype_marker(readtype)

    assert marker['glyph'] == glyph
    assert marker['label'] == label


@pytest.mark.parametrize('column', [
    'nreads_wholecounts_norm',
    'nreads_wholecounts_unique_norm',
    'wholecounts',
])
def test_a_full_obs_column_resolves_to_the_same_marker_as_its_bare_read_type(column):
    """The abandoned attempt listed each basis spelling separately and missed one entirely.
    Normalizing instead means a new basis cannot desynchronise the map."""
    assert plotsPalette.readtype_marker(column)['glyph'] == '◆'


@pytest.mark.parametrize('readtype', ['wholeprecounts', 'partialprecounts', 'trailercounts'])
def test_the_pre_trna_family_shares_one_glyph(readtype):
    """The published figure labels a single cross 'Pre-tRNA counts'. Three near-identical
    shapes nobody can distinguish at 8pt would be worse than one honest family marker."""
    marker = plotsPalette.readtype_marker(readtype)

    assert marker['glyph'] == '✕'
    assert marker['label'] == 'Pre-tRNA counts'


def test_an_unmapped_read_type_gets_a_neutral_fallback_rather_than_an_error():
    """`antisense` is legitimately requestable through --diffrts, so it must render as
    something rather than aborting a plot."""
    marker = plotsPalette.readtype_marker('antisense')

    assert marker['glyph'] == '○'
    assert marker['filled'] is False
    assert 'antisense' in marker['label'].lower()


def test_every_marker_carries_a_matplotlib_marker_code_for_legends():
    """The legend is built from matplotlib markers while the labels carry unicode glyphs, so
    the two have to be kept beside each other or they drift into different shapes."""
    for readtype in ('fiveprime', 'threeprime', 'other', 'wholecounts', 'trailercounts',
                     'antisense'):
        assert plotsPalette.readtype_marker(readtype)['marker']
