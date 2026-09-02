"""Declared order has to survive the .h5ad, or it is not a source of truth.

The order is applied once at build time and read back by every later plotting, clustering and
DE call, so the whole design rests on `ordered=True` and the category sequence surviving a
write/read cycle. If AnnData or its HDF5 encoding dropped either, the failure would be silent
and would surface only as a legend quietly reverting to alphabetical -- long after the run that
declared it.

A characterization test rather than a driven one: it pins behaviour this design depends on but
does not itself implement.
"""
import os

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from trnagraph.modules import toolsTG


DECLARED = ['D0', 'D7', 'D14', 'D35']


def _object():
    obs = pd.DataFrame({'timepoint': ['D35', 'D0', 'D14', 'D7'],
                        'sample': ['a', 'b', 'c', 'd']})
    toolsTG.apply_category_order(obs, {'timepoint': DECLARED})
    return ad.AnnData(X=np.zeros((4, 2), dtype='float32'), obs=obs)


def test_ordered_flag_and_category_sequence_survive_a_write_read_cycle(tmp_path):
    path = os.path.join(tmp_path, 'ordered.h5ad')
    _object().write_h5ad(path)

    reloaded = ad.read_h5ad(path)

    assert reloaded.obs['timepoint'].cat.ordered, 'ordered=True is what makes it an ORDER'
    assert list(reloaded.obs['timepoint'].cat.categories) == DECLARED
