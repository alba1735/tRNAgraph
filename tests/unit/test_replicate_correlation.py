"""Replicate-correlation QC: do samples agree with their own group more than with others?

One number decides whether an experiment's grouping is supported by its data: the gap between
mean within-group and mean between-group correlation. A near-zero gap means the groups are not
separable at all; a sample whose correlation with its own group sits well below its peers' is a
candidate for exclusion.

It also gives deduplication an objective pass/fail. Deduplication removes a technical artifact,
so it should make replicates *more* alike and groups *more* distinct — on the human OTTR dataset
the gap widened from +0.0663 to +0.0744. A run where the gap narrows is a warning that
deduplication removed signal rather than noise.
"""
import numpy as np
import pandas as pd
import pytest

from trnagraph.modules.toolsTG import replicate_correlation


def _counts(values, samples, groups, features):
    """Long-form obs frame of per-feature counts, as adata.obs carries them."""
    rows = []
    for si, s in enumerate(samples):
        for fi, f in enumerate(features):
            rows.append({'sample': s, 'group': groups[si], 'trna': f,
                         'nreads_total_unique_norm': values[si][fi]})
    return pd.DataFrame(rows)


def test_separates_groups_that_genuinely_differ():
    features = [f'tRNA-{i}' for i in range(30)]
    rng = np.random.default_rng(0)
    a = rng.lognormal(3, 1, 30)
    b = rng.lognormal(3, 1, 30)
    # Two groups of three, replicates differing only by small noise.
    vals = [a * rng.uniform(0.95, 1.05, 30) for _ in range(3)] + \
           [b * rng.uniform(0.95, 1.05, 30) for _ in range(3)]
    obs = _counts(vals, ['A1','A2','A3','B1','B2','B3'], ['A','A','A','B','B','B'], features)

    res = replicate_correlation(obs, group_col='group')

    assert res['summary']['within_mean'] > res['summary']['between_mean']
    assert res['summary']['separation'] == pytest.approx(
        res['summary']['within_mean'] - res['summary']['between_mean'])
    assert set(res['per_sample']['sample']) == {'A1','A2','A3','B1','B2','B3'}


def test_reports_no_separation_when_groups_are_arbitrary():
    features = [f'tRNA-{i}' for i in range(30)]
    rng = np.random.default_rng(1)
    base = rng.lognormal(3, 1, 30)
    vals = [base * rng.uniform(0.9, 1.1, 30) for _ in range(6)]
    obs = _counts(vals, ['A1','A2','A3','B1','B2','B3'], ['A','A','A','B','B','B'], features)

    res = replicate_correlation(obs, group_col='group')

    # Labels carry no information here, so within and between should be comparable.
    assert abs(res['summary']['separation']) < 0.05


def test_single_sample_groups_yield_no_within_group_pairs():
    """Every sample its own group -- the default trim_metadata.tsv mistake."""
    features = [f'tRNA-{i}' for i in range(20)]
    rng = np.random.default_rng(2)
    vals = [rng.lognormal(3, 1, 20) for _ in range(3)]
    obs = _counts(vals, ['S1','S2','S3'], ['S1','S2','S3'], features)

    res = replicate_correlation(obs, group_col='group')

    assert res['summary']['n_within_pairs'] == 0
    assert np.isnan(res['summary']['within_mean'])
