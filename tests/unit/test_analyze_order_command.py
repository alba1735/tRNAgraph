"""`analyze order` applies a declared category order to an object already built.

The order is written into the object rather than resolved per call, so an object built before
the declaration existed -- or before the maintainer knew what the order should be -- would
otherwise need a full rebuild to gain it. The human dataset is exactly that case: rebuilding it
is not practical, and its BAM/FASTQ inputs are treated as read-only.

An in-place command on the `addsplit` pattern rather than a side effect of `graph`: rewriting
obs dtypes is a far larger mutation than appending a cache entry, and as a side effect it would
fire differently depending on which config a given run happened to pass.
"""
import os

import anndata as ad
import numpy as np
import pandas as pd
from types import SimpleNamespace

from trnagraph.modules import adataBuild


DECLARED = ['D0', 'D7', 'D14', 'D35']


def _write_object(path):
    obs = pd.DataFrame({'timepoint': ['D35', 'D0', 'D14', 'D7'],
                        'sample': ['a', 'b', 'c', 'd']})
    ad.AnnData(X=np.zeros((4, 1), dtype='float32'), obs=obs).write_h5ad(path)


def test_order_is_applied_in_place(tmp_path):
    path = os.path.join(tmp_path, 'obj.h5ad')
    _write_object(path)

    adataBuild.apply_order(SimpleNamespace(
        anndata=path, order={'timepoint': DECLARED}, output=None))

    reloaded = ad.read_h5ad(path)
    assert reloaded.obs['timepoint'].cat.ordered
    assert list(reloaded.obs['timepoint'].cat.categories) == DECLARED


def test_output_writes_a_copy_and_leaves_the_input_alone(tmp_path):
    source = os.path.join(tmp_path, 'src.h5ad')
    destination = os.path.join(tmp_path, 'dst.h5ad')
    _write_object(source)

    adataBuild.apply_order(SimpleNamespace(
        anndata=source, order={'timepoint': DECLARED}, output=destination))

    assert list(ad.read_h5ad(destination).obs['timepoint'].cat.categories) == DECLARED
    assert not isinstance(ad.read_h5ad(source).obs['timepoint'].dtype, pd.CategoricalDtype), \
        'the input is untouched when an explicit output was given'
