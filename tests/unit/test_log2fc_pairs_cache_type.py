"""Regression tests for the type of `pairs` on the two paths out of adataLog2FC.main().

A fresh computation puts a `list` of tuples at uns['log2FC'][config][compare]['pairs'], but an
h5ad round-trip turns that into a numpy array of shape (n, 2). Every later run reads the cached
value, so consumers saw a different type on a warm object than on a cold one -- and
plotsVolcano's `if pairs:` raised `ValueError: The truth value of an array with more than one
element is ambiguous` on every cached run, which meant the combined volcano overview PDF was only
ever written the FIRST time a given config_name/cutoff was graphed.

Tested at two seams: adataLog2FC.main()'s returned dict (the cause) and plotsVolcano.visualizer()
(the symptom). A test that only exercised a cold, in-process object would have passed against the
broken code, which is how this shipped.
"""
import glob
import os

import anndata as ad
import numpy as np
import pandas as pd

from trnagraph.modules.plotsVolcano import visualizer
from trnagraph.modules.toolsTG import adataLog2FC

READTYPES = ['nreads_total_unique_norm']
CUTOFF = 5


def _make_adata():
    # One differential tRNA among several flat ones, as in test_log2fc.py: a lone differential
    # feature destabilizes PyDESeq2's dispersion-trend fit, which assumes most are not.
    rng = np.random.default_rng(0)
    groups = ['A', 'B']
    samples = [f'{g}_rep{i}' for g in groups for i in range(3)]
    sample_group = {s: s.rsplit('_rep', 1)[0] for s in samples}
    trna_means = {'trnaHigh': {'A': 800, 'B': 50}}
    trna_means.update({f'trnaFlat{i}': {'A': 400, 'B': 400} for i in range(6)})

    rows = []
    for trna, group_means in trna_means.items():
        for sample in samples:
            group = sample_group[sample]
            mean = group_means[group]
            raw_tu = rng.negative_binomial(n=10, p=10 / (10 + mean))
            raw_t = rng.negative_binomial(n=10, p=10 / (10 + mean))
            rows.append({
                'trna': trna, 'sample': sample, 'group': group,
                'nreads_total_unique_raw': raw_tu, 'nreads_total_unique_norm': raw_tu,
                'nreads_total_raw': raw_t, 'nreads_total_norm': raw_t,
            })
    obs = pd.DataFrame(rows)
    obs.index = [f"{r.trna}_{r['sample']}" for _, r in obs.iterrows()]
    return ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs)


def _roundtrip(adata, tmp_path, name='obj.h5ad'):
    """Write and re-read, which is what turns the cached `pairs` list into an ndarray."""
    path = tmp_path / name
    adata.write_h5ad(path)
    return ad.read_h5ad(path)


def _overview_files(output_dir):
    """Only the top-level overview PDFs; the per-pair plots live in individual/."""
    return sorted(os.path.basename(p) for p in glob.glob(os.path.join(output_dir, '*.pdf')))


def test_cached_pairs_have_the_same_type_as_freshly_computed_pairs(tmp_path):
    adata = _make_adata()
    _, cold_dict = adataLog2FC(adata, 'group', READTYPES[0], readcount_cutoff=CUTOFF).main()
    cold_pairs = cold_dict['default']['group']['pairs']

    _, warm_dict = adataLog2FC(
        _roundtrip(adata, tmp_path), 'group', READTYPES[0], readcount_cutoff=CUTOFF
    ).main()
    warm_pairs = warm_dict['default']['group']['pairs']

    assert warm_pairs == cold_pairs, (
        f'cached pairs {warm_pairs!r} ({type(warm_pairs).__name__}) differ from freshly '
        f'computed {cold_pairs!r} ({type(cold_pairs).__name__})'
    )
    assert isinstance(warm_pairs, list)
    assert all(isinstance(pair, tuple) for pair in warm_pairs)
    assert all(isinstance(level, str) for pair in warm_pairs for level in pair)


def test_combined_overview_is_written_on_a_second_run_against_a_cached_object(tmp_path):
    adata = _make_adata()
    first = str(tmp_path / 'first') + os.sep
    visualizer(adata, 'group', READTYPES, CUTOFF, first, threaded=False, is_full_variant=True)
    assert _overview_files(first) == [f'group_combined_{CUTOFF}_volcano.pdf']

    second = str(tmp_path / 'second') + os.sep
    visualizer(_roundtrip(adata, tmp_path), 'group', READTYPES, CUTOFF, second,
               threaded=False, is_full_variant=True)
    assert _overview_files(second) == [f'group_combined_{CUTOFF}_volcano.pdf'], (
        'the combined overview was not written from the cached log2FC entry'
    )
