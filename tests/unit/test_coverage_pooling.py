"""Coverage's pools must not pickle the visualizer into every task.

`Pool.map(self.some_method, items)` and `Pool.map(partial(self.some_method, ...), items)` both
pickle `self` for each task, and this `self` carries the whole AnnData. Measured on the hg38
object that is **549 MB per task**, against 436 tasks for the per-tRNA split and 28 per
combined PDF -- so the run spends its time pushing gigabytes through a pipe to do about 1.3
seconds of drawing per page. `generate_combine` did not finish within 10 minutes; after the
fix it takes 50 seconds.

Under `fork` a worker already shares the parent's memory, so the visualizer is published to a
module global before the Pool is created and nothing large is pickled. The regression is easy
to reintroduce by writing the obvious thing, and expensive to notice (it looks like a hang on
large data and is invisible on small), so the shape is checked statically.
"""
import ast
import multiprocessing
import pathlib

import pytest

from trnagraph.modules import plotsCoverage

SOURCE = pathlib.Path('src/trnagraph/modules/plotsCoverage.py')
POOL_METHODS = {'map', 'imap', 'imap_unordered', 'starmap', 'apply_async', 'map_async'}


def _pool_call_first_args():
    """(lineno, first positional argument node) for every pool submission in the module."""
    found = []
    for node in ast.walk(ast.parse(SOURCE.read_text())):
        if not isinstance(node, ast.Call) or not node.args:
            continue
        func = node.func
        if isinstance(func, ast.Attribute) and func.attr in POOL_METHODS:
            found.append((node.lineno, node.args[0]))
    return found


def _carries_self(node):
    """True if this callable expression closes over `self`."""
    if isinstance(node, ast.Attribute):
        return any(isinstance(n, ast.Name) and n.id == 'self' for n in ast.walk(node))
    if isinstance(node, ast.Call):  # e.g. partial(self.method, ...)
        return any(isinstance(n, ast.Name) and n.id == 'self' for n in ast.walk(node))
    return False


def test_no_pool_task_closes_over_self():
    offenders = [line for line, arg in _pool_call_first_args() if _carries_self(arg)]

    assert not offenders, (
        f'plotsCoverage lines {offenders} submit a callable that carries `self` to a pool, so '
        f'the whole AnnData is pickled for every task. Publish the visualizer with '
        f'_share_with_workers() and submit a module-level worker function instead.'
    )


def test_the_check_recognises_both_broken_shapes():
    """Guard the guard: both forms that shipped must be rejected by the detector."""
    for src in ('p.map(self.generate_combine_page, items)',
                "p.map(partial(self.generate_combine_page, coverage_fill='ci'), items)",
                'p.imap_unordered(self.generate_split_single, ulist)'):
        call = ast.parse(src).body[0].value
        assert _carries_self(call.args[0]), f'detector missed {src}'


def test_a_module_level_worker_is_accepted():
    call = ast.parse('p.map(_combine_page_worker, items)').body[0].value

    assert not _carries_self(call.args[0])


def test_the_worker_functions_exist_and_read_the_shared_global():
    for name in ('_combine_page_worker', '_split_single_worker'):
        assert callable(getattr(plotsCoverage, name))
    source = SOURCE.read_text()
    assert source.count('_WORKER_VISUALIZER') >= 3


def test_sharing_reports_whether_workers_can_inherit():
    """
    Only `fork` gives a child the parent's globals. Under spawn/forkserver the worker would
    find None, so the caller has to fall back to rendering serially rather than fail obscurely
    inside a pool.
    """
    shared = plotsCoverage._share_with_workers(object())

    assert shared == (multiprocessing.get_start_method(allow_none=False) == 'fork')
    assert plotsCoverage._WORKER_VISUALIZER is not None


# The function that actually dispatches to the pool for each path: generate_combine delegates
# its page rendering to _render_pages, while generate_split dispatches inline.
@pytest.mark.parametrize('method', ['_render_pages', 'generate_split'])
def test_both_pooled_paths_have_a_serial_fallback(method):
    """The fallback is what keeps a non-fork platform correct rather than merely slower."""
    text = SOURCE.read_text()
    source = ast.get_source_segment(text, next(
        n for n in ast.walk(ast.parse(text))
        if isinstance(n, ast.FunctionDef) and n.name == method))

    assert 'not shared' in source and 'self.threads <= 1' in source, (
        f'{method} dispatches to a Pool with no serial path for platforms where workers '
        f'cannot inherit the shared visualizer.'
    )
