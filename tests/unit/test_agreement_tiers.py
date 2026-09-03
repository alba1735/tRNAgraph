"""How many contrasts a feature agrees across, and which way it moved.

Colour on the agreement volcano encodes a count: significant in 3 of 3, 2 of 3, 1 of 3. The
count only means something if the contrasts agree in DIRECTION -- a tRNA up at one timepoint and
down at another has not been consistently anything, and colouring it as a 2-of-3 responder would
say it had.

Direction comes from the contrast on the axis, which is what puts the point on its side of zero.
"""
import numpy as np
import pandas as pd

from trnagraph.modules import plotsAgreement


REF = 'Day 0'


def _frame(features, per_contrast):
    """A log2FC frame shaped like toolsTG's: log2_{a}-{b} and pval_{a}-{b} columns."""
    data = {}
    for (a, b), (log2fc, pval) in per_contrast.items():
        data[f'log2_{a}-{b}'] = pd.Series(log2fc, index=features)
        data[f'pval_{a}-{b}'] = pd.Series(pval, index=features)
    return pd.DataFrame(data, index=features)


def _table(frame, contrasts, drawn, **kw):
    return plotsAgreement.agreement_table({'fiveprime': frame}, contrasts, drawn,
                                          log2fc=1.5, padj=0.001, **kw)


def test_a_feature_significant_in_every_contrast_tops_the_tiers():
    contrasts = [(REF, 'Day 35'), (REF, 'Day 70')]
    frame = _frame(['tRNA-A'], {
        (REF, 'Day 35'): ([3.0], [1e-6]),
        (REF, 'Day 70'): ([4.0], [1e-8]),
    })

    row = _table(frame, contrasts, (REF, 'Day 70')).iloc[0]

    assert row['n_agree'] == 2
    assert row['n_contrasts'] == 2


def test_a_contrast_that_moves_the_other_way_does_not_count_as_agreement():
    """The point of the figure: consistency, not merely repeated significance."""
    contrasts = [(REF, 'Day 35'), (REF, 'Day 70')]
    frame = _frame(['tRNA-A'], {
        (REF, 'Day 35'): ([-3.0], [1e-6]),   # significant, opposite direction
        (REF, 'Day 70'): ([4.0], [1e-8]),
    })

    row = _table(frame, contrasts, (REF, 'Day 70')).iloc[0]

    assert row['n_agree'] == 1


def test_a_contrast_below_either_threshold_does_not_count():
    contrasts = [(REF, 'Day 35'), (REF, 'Day 70')]
    frame = _frame(['tRNA-A'], {
        (REF, 'Day 35'): ([4.0], [0.01]),    # big effect, not significant
        (REF, 'Day 70'): ([4.0], [1e-8]),
    })

    row = _table(frame, contrasts, (REF, 'Day 70')).iloc[0]

    assert row['n_agree'] == 1


def test_tiers_at_three_contrasts_reproduce_the_reference_figure():
    """Three anchored contrasts give three tiers: 3 of 3, 2 of 3, 1 of 3."""
    contrasts = [(REF, 'Day 14'), (REF, 'Day 35'), (REF, 'Day 70')]
    sig, quiet = (4.0, 1e-8), (0.1, 0.9)
    frame = _frame(['all', 'two', 'one'], {
        (REF, 'Day 14'): ([sig[0], quiet[0], quiet[0]], [sig[1], quiet[1], quiet[1]]),
        (REF, 'Day 35'): ([sig[0], sig[0], quiet[0]], [sig[1], sig[1], quiet[1]]),
        (REF, 'Day 70'): ([sig[0], sig[0], sig[0]], [sig[1], sig[1], sig[1]]),
    })

    table = _table(frame, contrasts, (REF, 'Day 70')).set_index('feature')

    assert list(table.loc[['all', 'two', 'one'], 'n_agree']) == [3, 2, 1]
    assert set(table['n_contrasts']) == {3}


def test_tiers_at_five_contrasts_degrade_to_a_ramp_rather_than_mislabelled_buckets():
    contrasts = [(REF, f'L{i}') for i in range(1, 6)]
    features = [f'f{n}' for n in range(1, 6)]
    per_contrast = {}
    for i, contrast in enumerate(contrasts, start=1):
        # Feature fN is significant in the first N contrasts.
        log2fc = [4.0 if n >= i else 0.1 for n in range(1, 6)]
        pval = [1e-8 if n >= i else 0.9 for n in range(1, 6)]
        per_contrast[contrast] = (log2fc, pval)
    frame = _frame(features, per_contrast)

    table = _table(frame, contrasts, contrasts[0]).set_index('feature')

    assert list(table.loc[features, 'n_agree']) == [1, 2, 3, 4, 5]
    assert set(table['n_contrasts']) == {5}


def test_direction_comes_from_the_contrast_on_the_axis():
    contrasts = [(REF, 'Day 35'), (REF, 'Day 70')]
    frame = _frame(['up', 'down'], {
        (REF, 'Day 35'): ([4.0, -4.0], [1e-8, 1e-8]),
        (REF, 'Day 70'): ([4.0, -4.0], [1e-8, 1e-8]),
    })

    table = _table(frame, contrasts, (REF, 'Day 70')).set_index('feature')

    assert table.loc['up', 'direction'] == 'Day 70'
    assert table.loc['down', 'direction'] == REF


def test_a_feature_not_significant_in_the_drawn_contrast_has_no_direction():
    """It is drawn, in the neutral colour -- it is a real measurement, just not a call."""
    contrasts = [(REF, 'Day 35'), (REF, 'Day 70')]
    frame = _frame(['quiet'], {
        (REF, 'Day 35'): ([4.0], [1e-8]),
        (REF, 'Day 70'): ([0.1], [0.9]),
    })

    row = _table(frame, contrasts, (REF, 'Day 70')).iloc[0]

    assert row['direction'] is None
    assert row['n_agree'] == 0


def test_features_deseq2_could_not_call_are_dropped():
    """NaN means independent filtering excluded it, not that it was insignificant."""
    contrasts = [(REF, 'Day 70')]
    frame = _frame(['real', 'unfitted'], {
        (REF, 'Day 70'): ([4.0, np.nan], [1e-8, np.nan]),
    })

    table = _table(frame, contrasts, (REF, 'Day 70'))

    assert list(table['feature']) == ['real']


def test_every_read_type_given_appears_with_its_own_rows():
    contrasts = [(REF, 'Day 70')]
    frame = _frame(['tRNA-A'], {(REF, 'Day 70'): ([4.0], [1e-8])})
    table = plotsAgreement.agreement_table(
        {'fiveprime': frame, 'threeprime': frame}, contrasts, (REF, 'Day 70'),
        log2fc=1.5, padj=0.001)

    assert sorted(table['readtype']) == ['fiveprime', 'threeprime']


def test_the_axis_call_uses_the_loose_threshold_not_the_strict_one():
    """Two thresholds, two jobs. 0.05 decides whether a point is called on the contrast being
    drawn -- so nothing below that line is ever coloured -- while 0.001 decides whether a
    contrast counts toward agreement. A point between the lines is a real call on this axis
    that was not strong anywhere, and it is drawn as such."""
    contrasts = [(REF, 'Day 35'), (REF, 'Day 70')]
    frame = _frame(['between'], {
        (REF, 'Day 35'): ([4.0], [0.02]),
        (REF, 'Day 70'): ([4.0], [0.02]),
    })

    row = _table(frame, contrasts, (REF, 'Day 70')).iloc[0]

    assert row['direction'] == 'Day 70', 'clears 0.05 on the drawn axis, so it is called'
    assert row['n_agree'] == 0, 'clears 0.001 nowhere, so nothing agrees'


def test_a_point_below_the_loose_threshold_is_not_called_at_all():
    contrasts = [(REF, 'Day 70')]
    frame = _frame(['faint'], {(REF, 'Day 70'): ([4.0], [0.2])})

    row = _table(frame, contrasts, (REF, 'Day 70')).iloc[0]

    assert row['direction'] is None


def test_agreement_still_needs_the_strict_threshold():
    contrasts = [(REF, 'Day 35'), (REF, 'Day 70')]
    frame = _frame(['mixed'], {
        (REF, 'Day 35'): ([4.0], [0.02]),    # called, but not strong
        (REF, 'Day 70'): ([4.0], [1e-9]),    # strong
    })

    row = _table(frame, contrasts, (REF, 'Day 70')).iloc[0]

    assert row['n_agree'] == 1
