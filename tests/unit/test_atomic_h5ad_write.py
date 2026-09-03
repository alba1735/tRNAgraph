"""Writing an .h5ad must never destroy the existing one.

anndata's own `.write()` truncates the target and streams into it, so an interruption -- a
crash, a full disk, two processes colliding on the file lock -- leaves a file that is neither
the old object nor the new one. It is not merely stale: anndata cannot open it at all, and the
build that produced it may be hours of compute.

Every write goes to a temporary file beside the target and is then renamed over it. os.replace
is atomic within a filesystem, so a reader sees either the whole old object or the whole new
one, and a failure leaves the old one exactly as it was.
"""
import os
import stat

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from trnagraph.modules import toolsTG


def _adata(marker='original'):
    return ad.AnnData(X=np.zeros((2, 2)),
                      obs=pd.DataFrame({'tag': [marker, marker]}, index=['a', 'b']))


def test_it_writes_a_readable_object(tmp_path):
    path = tmp_path / 'x.h5ad'
    toolsTG.write_h5ad(_adata(), path)

    assert list(ad.read_h5ad(path).obs['tag']) == ['original', 'original']


def test_a_failed_write_leaves_the_existing_object_untouched(tmp_path):
    """The whole point. Previously the target was already truncated by the time this failed."""
    path = tmp_path / 'x.h5ad'
    toolsTG.write_h5ad(_adata('original'), path)
    before = path.read_bytes()

    class Exploding:
        def write(self, target):
            open(target, 'wb').write(b'partial')
            raise RuntimeError('disk full')

    with pytest.raises(RuntimeError):
        toolsTG.write_h5ad(Exploding(), path)

    assert path.read_bytes() == before
    assert list(ad.read_h5ad(path).obs['tag']) == ['original', 'original']


def test_a_failed_write_leaves_no_temporary_file_behind(tmp_path):
    path = tmp_path / 'x.h5ad'
    toolsTG.write_h5ad(_adata(), path)

    class Exploding:
        def write(self, target):
            raise RuntimeError('nope')

    with pytest.raises(RuntimeError):
        toolsTG.write_h5ad(Exploding(), path)

    assert [p.name for p in tmp_path.iterdir()] == ['x.h5ad']


def test_a_successful_write_leaves_no_temporary_file_behind(tmp_path):
    path = tmp_path / 'x.h5ad'
    toolsTG.write_h5ad(_adata(), path)

    assert [p.name for p in tmp_path.iterdir()] == ['x.h5ad']


def test_the_temporary_file_sits_beside_the_target(tmp_path):
    """os.replace is only atomic within one filesystem, so the temporary cannot go to /tmp."""
    seen = {}
    path = tmp_path / 'x.h5ad'

    class Watching:
        def write(self, target):
            seen['dir'] = os.path.dirname(os.path.abspath(target))
            _adata().write(target)

    toolsTG.write_h5ad(Watching(), path)

    assert seen['dir'] == str(tmp_path)


def test_an_overwrite_keeps_the_original_file_permissions(tmp_path):
    """mkstemp creates 0600. Renaming that over a group-readable object on a shared server
    would silently make it unreadable to everyone else."""
    path = tmp_path / 'x.h5ad'
    toolsTG.write_h5ad(_adata(), path)
    os.chmod(path, 0o644)

    toolsTG.write_h5ad(_adata('second'), path)

    assert stat.S_IMODE(os.stat(path).st_mode) == 0o644


def test_a_new_file_is_not_created_private(tmp_path):
    path = tmp_path / 'fresh.h5ad'
    toolsTG.write_h5ad(_adata(), path)

    assert stat.S_IMODE(os.stat(path).st_mode) & 0o044, 'should be readable, not mkstemp 0600'


def test_no_module_persists_an_anndata_object_with_a_bare_write():
    """A new call site that bypasses the helper reintroduces the corruption, silently and only
    under interruption. Scanned rather than trusted, since the failure never shows up in a
    passing run."""
    import ast
    import pathlib

    root = pathlib.Path(toolsTG.__file__).parent
    offenders = []
    for path in sorted(root.glob('*.py')) + [root.parent / 'cli.py']:
        if path.name == 'toolsTG.py':
            continue  # where the helper itself lives
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            if not (isinstance(func, ast.Attribute) and func.attr == 'write'):
                continue
            # `<something>adata<something>.write(path)` -- the AnnData persistence shape.
            target = ast.unparse(func.value)
            if 'adata' in target.lower():
                offenders.append(f'{path.name}:{node.lineno} {target}.write(...)')

    assert not offenders, (
        'these persist an AnnData object directly; route them through toolsTG.write_h5ad so an '
        f'interrupted write cannot corrupt the object: {offenders}')
