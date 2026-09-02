"""uns['multivariate'] stores membership and provenance, and survives the .h5ad.

Two things this store must get right, both of which have bitten this codebase before.

FIRST, uns round-trips lossily. Python sets do not survive at all and lists come back as numpy
arrays of numpy strings, so membership has to be normalised on the way in and on the way out.
adataLog2FC already carries _as_pair_list() for exactly this reason: a cached `pairs` written as
tuples returns as an ndarray, and a consumer doing `if pairs:` raised on the array form.

SECOND, a cache key must include everything that changes the value. The `shrink` key was added
to uns['log2FC'] on that principle -- an entry written before it existed reads as None and
correctly recomputes, rather than serving fold changes from the wrong estimator. Membership
depends on more parameters than a fold change does (the significance thresholds, the read-count
floor, the grouping and reference level), so a stored set whose provenance no longer matches
what is being asked for must MISS rather than be served.

The store deliberately holds no copy of the fold-change frames themselves: those stay
authoritative in uns['log2FC'], so there is exactly one place to invalidate.
"""
import os

import anndata as ad
import numpy as np
import pandas as pd

from trnagraph.modules import toolsTG


PROVENANCE = {'grouping': 'group', 'reference': 'ctrl', 'readtype': 'nreads_total_norm',
              'cutoff': 80.0, 'shrink': 'apeGLM', 'log2fc': 1.5, 'padj': 0.001}
SETS = {'u60': {'tRNA-Gly-GCC-1', 'tRNA-Ala-AGC-1'}, 'o60': {'tRNA-Ala-AGC-1'}}


def _adata():
    obs = pd.DataFrame({'trna': ['a', 'b'], 'sample': ['s1', 's2']}, index=['o0', 'o1'])
    return ad.AnnData(X=np.zeros((2, 1), dtype='float32'), obs=obs)


def test_sets_are_stored_as_sorted_lists():
    """Python sets do not survive h5ad, and an arbitrary order makes diffs unreadable."""
    adata = _adata()

    toolsTG.write_multivariate(adata, 'default', 'venn', 'full', SETS, PROVENANCE)

    stored = adata.uns['multivariate']['default']['venn']['full']['sets']
    assert stored['u60'] == ['tRNA-Ala-AGC-1', 'tRNA-Gly-GCC-1']
    assert isinstance(stored['u60'], list)


def test_membership_survives_a_write_read_cycle(tmp_path):
    adata = _adata()
    toolsTG.write_multivariate(adata, 'default', 'venn', 'full', SETS, PROVENANCE)
    path = os.path.join(tmp_path, 'mv.h5ad')
    adata.write_h5ad(path)

    reloaded = ad.read_h5ad(path)
    sets = toolsTG.read_multivariate(reloaded, 'default', 'venn', 'full', PROVENANCE)

    assert sets == {'u60': ['tRNA-Ala-AGC-1', 'tRNA-Gly-GCC-1'], 'o60': ['tRNA-Ala-AGC-1']}


def test_a_changed_parameter_misses_rather_than_serving_a_stale_set():
    """The whole point of storing provenance: a set computed at padj <= 0.001 is not an answer
    to a question asked at padj <= 0.05."""
    adata = _adata()
    toolsTG.write_multivariate(adata, 'default', 'venn', 'full', SETS, PROVENANCE)

    asked = dict(PROVENANCE, padj=0.05)
    assert toolsTG.read_multivariate(adata, 'default', 'venn', 'full', asked) is None


def test_an_absent_entry_misses():
    assert toolsTG.read_multivariate(_adata(), 'default', 'venn', 'full', PROVENANCE) is None


def test_tags_are_kept_apart():
    """'full' is a REAL key here, unlike in uns['size_splits'] where it is a reserved pseudo-tag
    meaning the unsuffixed location -- membership for the full variant has to live somewhere."""
    adata = _adata()

    toolsTG.write_multivariate(adata, 'default', 'venn', 'full', {'x': {'a'}}, PROVENANCE)
    toolsTG.write_multivariate(adata, 'default', 'venn', 'u60', {'x': {'b'}}, PROVENANCE)

    assert toolsTG.read_multivariate(adata, 'default', 'venn', 'full', PROVENANCE) == {'x': ['a']}
    assert toolsTG.read_multivariate(adata, 'default', 'venn', 'u60', PROVENANCE) == {'x': ['b']}
