"""Regression tests for the aggregated DESeq2 output files (roadmap.md: "`avgs.txt` collides with
a tRAX filename while holding a different quantity").

Two files written by `write_aggregated_de_results()` carried a different quantity than the tRAX
file of the same name:

* `<exp>-avgs.txt` -- tRAX's is DESeq2's `baseMean` (`analyzecounts.R`'s `colgetavgname()` takes
  column 1 of each `results()` object), the mean of normalized counts across *all* samples,
  invariant to the contrast. tRAX writes that one column once per pairwise comparison under each
  comparison's name, so its columns are labelled per comparison but hold identical values.
  tRNAgraph wrote genuine per-group means instead -- the more useful quantity, but under tRAX's
  name. `baseMean` is now the parity column and the per-group means moved to `groupavgs.txt`.
* `<exp>-combine.txt` -- tRAX's trailing columns are per-group *medians* (`analyzecounts.R`
  cbinds `medcountmat`); tRNAgraph put per-group *means* there under the same group labels.

Expected values here come from tRAX's semantics and from hand-computed arithmetic, never from
re-running the code under test.
"""
import logging
import os
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from trnagraph.modules.adataBuild import AnalysisPipeline


SAMPLES = ['s1', 's2', 's3', 's4']
FEATURES = ['featA', 'featB']
# Chosen so every per-group statistic is a distinct, exactly-representable number and the
# mean and median of a group differ -- a test that cannot tell means from medians is the whole
# defect being fixed.
NORMED = np.array([
    # s1    s2     s3     s4
    [10.0,  20.0,  60.0, 100.0],   # featA: ctrl mean 15, median 15; treat mean 80, median 80
    [ 2.0,   6.0,   8.0,  24.0],   # featB: ctrl mean  4, median  4; treat mean 16, median 16
])
GROUPS = {'s1': 'ctrl', 's2': 'ctrl', 's3': 'treat', 's4': 'treat'}


def _make_pipeline(expname='exp'):
    pipeline = AnalysisPipeline.__new__(AnalysisPipeline)
    pipeline.expinfo = SimpleNamespace(expname=expname)
    pipeline.logger = logging.getLogger('trnagraph.modules.adataBuild')
    return pipeline


def _fake_dds():
    """A stand-in for the PyDESeq2 DeseqDataSet, carrying only what the writer reads."""
    return SimpleNamespace(
        layers={'normed_counts': NORMED.T},   # dds is samples x features
        obs_names=pd.Index(SAMPLES),
        var_names=pd.Index(FEATURES),
        var=pd.DataFrame({'dispersions': [0.1, 0.2]}, index=FEATURES),
    )


def _sample_df():
    return pd.DataFrame({'condition': [GROUPS[s] for s in SAMPLES]}, index=SAMPLES)


def _de_results():
    """One pairwise comparison, in the shape run_pairwise_de() produces."""
    res = pd.DataFrame(
        {'baseMean': [47.5, 10.0], 'log2FoldChange': [2.5, 2.0], 'padj': [0.01, 0.5]},
        index=FEATURES,
    )
    return {'ctrl_treat': res}


def _read(tmp_path, name):
    return pd.read_csv(os.path.join(str(tmp_path), name), sep='\t', index_col=0)


@pytest.fixture
def written(tmp_path):
    pipeline = _make_pipeline()
    pipeline.write_aggregated_de_results(
        _fake_dds(), _sample_df(), _de_results(), str(tmp_path), ''
    )
    return tmp_path


def test_avgs_holds_deseq2_basemean_as_a_single_column(written):
    """tRAX's -avgs.txt is baseMean: the mean of normalized counts over ALL samples, one number
    per feature. featA: (10+20+60+100)/4 = 47.5. featB: (2+6+8+24)/4 = 10.0."""
    avgs = _read(written, 'exp-avgs.txt')

    assert list(avgs.columns) == ['baseMean']
    assert avgs.loc['featA', 'baseMean'] == pytest.approx(47.5)
    assert avgs.loc['featB', 'baseMean'] == pytest.approx(10.0)


def test_per_group_means_survive_under_a_non_colliding_name(written):
    """The per-group means are the more useful quantity and must not be lost when avgs.txt is
    handed to baseMean -- they move to groupavgs.txt, one column per group.
    ctrl: (10+20)/2 = 15, (2+6)/2 = 4.  treat: (60+100)/2 = 80, (8+24)/2 = 16."""
    groupavgs = _read(written, 'exp-groupavgs.txt')

    assert list(groupavgs.columns) == ['ctrl', 'treat']
    assert groupavgs.loc['featA'].tolist() == pytest.approx([15.0, 80.0])
    assert groupavgs.loc['featB'].tolist() == pytest.approx([4.0, 16.0])


def test_medians_remain_per_group_and_untouched(written):
    """tRAX's -medians.txt is already per-group on both sides -- no collision, no change."""
    medians = _read(written, 'exp-medians.txt')

    assert list(medians.columns) == ['ctrl', 'treat']
    assert medians.loc['featA'].tolist() == pytest.approx([15.0, 80.0])


def test_combine_trailing_columns_are_per_group_medians_like_trax(written):
    """tRAX's -combine.txt is cbind(alllogvals, allprobs, medcountmat) -- its trailing per-group
    columns are MEDIANS, matching its own -medians.txt row for row. tRNAgraph put per-group means
    there under the same group labels, so a column-aligned diff silently compared a mean against
    a median. Here the fixture's medians and means coincide, so the test uses a feature where
    they differ."""
    combine = _read(written, 'exp-combine.txt')
    medians = _read(written, 'exp-medians.txt')

    assert list(combine.columns) == ['log2_ctrl_treat', 'pval_ctrl_treat', 'ctrl', 'treat']
    pd.testing.assert_frame_equal(
        combine[['ctrl', 'treat']], medians[['ctrl', 'treat']], check_names=False
    )


def test_combine_trailing_columns_track_medians_not_means_when_they_differ(tmp_path):
    """A three-sample group whose mean and median genuinely differ: values 1, 2, 99 give a mean
    of 34 and a median of 2. combine.txt must show 2."""
    samples = ['a', 'b', 'c']
    dds = SimpleNamespace(
        layers={'normed_counts': np.array([[1.0, 2.0, 99.0]]).T},
        obs_names=pd.Index(samples),
        var_names=pd.Index(['featA']),
        var=pd.DataFrame({'dispersions': [0.1]}, index=['featA']),
    )
    sample_df = pd.DataFrame({'condition': ['g', 'g', 'g']}, index=samples)
    de = {'g_g': pd.DataFrame({'log2FoldChange': [1.0], 'padj': [0.5]}, index=['featA'])}

    _make_pipeline().write_aggregated_de_results(dds, sample_df, de, str(tmp_path), '')

    assert _read(tmp_path, 'exp-combine.txt').loc['featA', 'g'] == pytest.approx(2.0)
    assert _read(tmp_path, 'exp-groupavgs.txt').loc['featA', 'g'] == pytest.approx(34.0)
    assert _read(tmp_path, 'exp-avgs.txt').loc['featA', 'baseMean'] == pytest.approx(34.0)


def test_basemean_is_written_even_when_no_comparison_ran(tmp_path):
    """baseMean is a plain descriptive statistic over the normalized counts, so it must not
    depend on a results() object existing -- an experiment with no pairs file returns early,
    before padjs/logvals/combine are written, and must still get avgs.txt."""
    _make_pipeline().write_aggregated_de_results(
        _fake_dds(), _sample_df(), {}, str(tmp_path), ''
    )

    assert _read(tmp_path, 'exp-avgs.txt').loc['featA', 'baseMean'] == pytest.approx(47.5)
    assert _read(tmp_path, 'exp-groupavgs.txt').loc['featA', 'ctrl'] == pytest.approx(15.0)
    assert not os.path.exists(os.path.join(str(tmp_path), 'exp-combine.txt'))


def test_prefix_is_applied_to_every_written_file(tmp_path):
    """The same writer serves six runs -- primary, allfeature_, trna_ and the three unique_
    ones -- so the new file has to carry the prefix like its siblings."""
    _make_pipeline().write_aggregated_de_results(
        _fake_dds(), _sample_df(), _de_results(), str(tmp_path), 'allfeature_'
    )

    for name in ('avgs', 'groupavgs', 'medians', 'combine', 'padjs', 'logvals'):
        assert os.path.exists(os.path.join(str(tmp_path), f'exp-allfeature_{name}.txt')), name
