"""Every recorded flag snapshot has to survive being written to the .h5ad.

`cli_specified` is bookkeeping for the --config merge -- a frozenset of the option names the
user actually typed -- and h5ad cannot write a frozenset. It rides on the args namespace, so any
code recording `vars(args)` into uns picks it up and the build dies at the moment of saving,
after all the real work is done.

This was found and fixed once, for uns['trnagraphruninfo']['flags'], and guarded by a test that
inspects the SOURCE of AnnDataBuilder.__init__ for the string "if k != 'cli_specified'". That
guard could not see the two other places that snapshot the same namespace --
`build_flags` for each read-length split variant -- so the identical bug survived there and
fired on the first `--readlengthsplit` build:

    IORegistryError: No method registered for writing <class 'frozenset'> into
    <class 'h5py._hl.group.Group'>
    [NOTE] Error raised while writing key 'cli_specified' ... to /uns/size_splits/u60/build_flags

So the guard here is behavioural -- actually write the object -- and applies to any
unwritable type, not just the one that happened to be found first.
"""
import os
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd

from trnagraph.modules import toolsTG


def _args():
    """An args namespace shaped like a real `analyze build`, bookkeeping included."""
    return SimpleNamespace(bamdir=None, readlengthsplit=60, threads=8,
                           cli_specified=frozenset({'bamdir', 'readlengthsplit'}))


def test_the_merge_bookkeeping_is_dropped():
    flags = toolsTG.sanitize_flags(_args())

    assert 'cli_specified' not in flags


def test_none_becomes_the_string_sentinel():
    """h5ad does not round-trip None reliably, so the existing convention is a sentinel."""
    assert toolsTG.sanitize_flags(_args())['bamdir'] == 'None'


def test_real_settings_survive():
    flags = toolsTG.sanitize_flags(_args())

    assert flags['readlengthsplit'] == 60
    assert flags['threads'] == 8


def test_a_split_variants_build_flags_can_actually_be_written(tmp_path):
    """The failure this exists to prevent, reproduced at the point it actually bites."""
    adata = ad.AnnData(X=np.zeros((1, 1), dtype='float32'),
                       obs=pd.DataFrame({'trna': ['t1']}, index=['o0']))
    adata.uns['size_splits'] = {'u60': {'build_flags': toolsTG.sanitize_flags(_args())}}

    adata.write_h5ad(os.path.join(tmp_path, 'split.h5ad'))


def test_every_flag_snapshot_in_adatabuild_goes_through_the_helper():
    """The generalisation of the source check this bug escaped: it is not enough for ONE
    function to filter, because the namespace is snapshotted in three places."""
    import inspect as _inspect
    import re

    from trnagraph.modules import adataBuild

    source = _inspect.getsource(adataBuild)
    snapshots = re.findall(r'^\s*(\w*flags\w*) = \{k: .*vars\(', source, re.MULTILINE)

    assert not snapshots, (
        f'{snapshots} build a flag snapshot from vars() inline instead of calling '
        f'toolsTG.sanitize_flags(), which is how the frozenset reached the h5ad')
