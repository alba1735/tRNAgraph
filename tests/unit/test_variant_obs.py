"""A read-length variant's obs, without materialising the variant.

toolsTG.build_variant_view() resolves a variant by taking a full adata.copy() and swapping .X
to the tag's layer. That is correct for `graph`, which resolves exactly one variant and then
plots from it -- and unusable for anything that needs two at once, because on the human object
it copies roughly 460 MB per call.

The cross-variant analyses need only the OBS side: identity columns, and that tag's per-obs
numeric columns, which already live in adata.obsm['size_split_<tag>'] under the same unsuffixed
names adata.obs uses for the full variant. So this reads them in place and leaves .X, the
layers and the rest of the object alone.

'full' is the reserved pseudo-tag meaning "the unsuffixed/default location", exactly as it is
for --variant, and is never a real key in uns['size_splits'].
"""
import anndata as ad
import numpy as np
import pandas as pd
import pytest

from trnagraph.modules import toolsTG


FULL_VALUES = [10.0, 20.0, 30.0, 40.0]
SPLIT_VALUES = [1.0, 2.0, 3.0, 4.0]


def _adata():
    obs = pd.DataFrame({
        'trna': ['tRNA-Ala-AGC-1', 'tRNA-Ala-AGC-1', 'tRNA-Gly-GCC-1', 'tRNA-Gly-GCC-1'],
        'sample': ['s1', 's2', 's1', 's2'],
        'group': ['a', 'b', 'a', 'b'],
        'nreads_total_norm': FULL_VALUES,
        'nreads_total_raw': FULL_VALUES,
    })
    obs.index = [f'o{i}' for i in range(4)]
    adata = ad.AnnData(X=np.zeros((4, 3), dtype='float32'), obs=obs)
    adata.obsm['size_split_u60'] = pd.DataFrame(
        {'nreads_total_norm': SPLIT_VALUES, 'nreads_total_raw': SPLIT_VALUES},
        index=obs.index)
    adata.uns['size_splits'] = {'u60': {}}
    return adata


def test_a_split_tag_takes_its_numeric_columns_from_obsm():
    obs = toolsTG.variant_obs(_adata(), 'u60')

    assert list(obs['nreads_total_norm']) == SPLIT_VALUES
    assert list(obs['nreads_total_raw']) == SPLIT_VALUES


def test_identity_columns_are_shared_across_variants():
    """trna/sample/group are not per-variant: a length cutoff changes which reads count, not
    which tRNAs or samples exist."""
    obs = toolsTG.variant_obs(_adata(), 'u60')

    assert list(obs['trna']) == ['tRNA-Ala-AGC-1', 'tRNA-Ala-AGC-1', 'tRNA-Gly-GCC-1', 'tRNA-Gly-GCC-1']
    assert list(obs['group']) == ['a', 'b', 'a', 'b']


def test_full_is_the_unsuffixed_default():
    obs = toolsTG.variant_obs(_adata(), 'full')

    assert list(obs['nreads_total_norm']) == FULL_VALUES


def test_an_unknown_tag_is_refused_naming_what_is_available():
    with pytest.raises(Exception) as excinfo:
        toolsTG.variant_obs(_adata(), 'u50')

    assert 'u50' in str(excinfo.value)
    assert 'u60' in str(excinfo.value), 'names the tags that DO exist'
