"""The stacked type plot gains a variant collapsing every non-tRNA type into one gray band.

`type_counts` carries 14 categories on the demo data, 9 of which are non-tRNA. For a figure
whose point is the tRNA composition, those 9 bands are noise competing for the reader's
attention with the categories that matter.

The tRNA/non-tRNA boundary is `adataBuild.TRNA_TYPE_LABELS`, reused rather than redefined --
it is the same definition the split-variant filters already use, and it deliberately keeps
`Mt_tRNA` on the tRNA side. The collapsed plot is emitted *in addition* to the full one.
"""
import pandas as pd

from trnagraph.modules.adataBuild import TRNA_TYPE_LABELS
from trnagraph.modules.plotsCount import NONTRNA_COLLAPSED_COLOR, NONTRNA_COLLAPSED_LABEL, collapse_nontrna_types


def _counts():
    return pd.DataFrame(
        {'s1': [100, 50, 10, 5, 3, 2, 1, 1], 's2': [200, 60, 20, 6, 4, 2, 2, 1]},
        index=['tRNA', 'pretRNA', 'Mt_tRNA', 'rRNA', 'miRNA', 'snRNA', 'snoRNA', 'other'],
    )


def test_non_trna_rows_collapse_into_one():
    collapsed = collapse_nontrna_types(_counts())

    assert list(collapsed.index) == ['tRNA', 'pretRNA', 'Mt_tRNA', NONTRNA_COLLAPSED_LABEL]


def test_the_collapsed_row_sums_the_rows_it_replaced():
    original = _counts()

    collapsed = collapse_nontrna_types(original)

    expected = original.loc[['rRNA', 'miRNA', 'snRNA', 'snoRNA', 'other']].sum()
    assert collapsed.loc[NONTRNA_COLLAPSED_LABEL].tolist() == expected.tolist()


def test_totals_are_preserved():
    """Collapsing is a regrouping, not a filter -- no reads may go missing."""
    original = _counts()

    collapsed = collapse_nontrna_types(original)

    assert collapsed.sum().tolist() == original.sum().tolist()


def test_mitochondrial_trna_stays_on_the_trna_side():
    collapsed = collapse_nontrna_types(_counts())

    assert 'Mt_tRNA' in collapsed.index
    assert collapsed.loc['Mt_tRNA'].tolist() == [10, 20]


def test_trna_rows_keep_their_original_order():
    """The collapsed band is appended, so the tRNA bands stack as they did before."""
    original = _counts()

    collapsed = collapse_nontrna_types(original)

    trna_rows = [i for i in collapsed.index if i != NONTRNA_COLLAPSED_LABEL]
    assert trna_rows == [i for i in original.index if i in TRNA_TYPE_LABELS]


def test_a_frame_with_no_non_trna_rows_is_unchanged():
    """vibrChol1 built without a GTF has nothing to collapse."""
    original = _counts().loc[['tRNA', 'pretRNA', 'Mt_tRNA']]

    collapsed = collapse_nontrna_types(original)

    pd.testing.assert_frame_equal(collapsed, original)


def test_a_frame_with_only_non_trna_rows_collapses_to_one_row():
    original = _counts().loc[['rRNA', 'miRNA']]

    collapsed = collapse_nontrna_types(original)

    assert list(collapsed.index) == [NONTRNA_COLLAPSED_LABEL]
    assert collapsed.loc[NONTRNA_COLLAPSED_LABEL].tolist() == original.sum().tolist()


def test_an_empty_frame_survives():
    original = pd.DataFrame(columns=['s1'])

    assert collapse_nontrna_types(original).empty


def test_the_collapsed_band_is_gray():
    assert NONTRNA_COLLAPSED_COLOR.lower() in ('#b0b0b0', '#808080', 'gray', 'grey')


def test_the_label_says_non_trna_not_small_rna():
    """Project terminology: these are non-tRNA features, never 'small RNA'."""
    assert 'small' not in NONTRNA_COLLAPSED_LABEL.lower()
    assert 'non' in NONTRNA_COLLAPSED_LABEL.lower()
