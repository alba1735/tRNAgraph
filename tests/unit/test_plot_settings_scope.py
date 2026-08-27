"""Every plot helper that reads `settings` must actually receive it.

Threading a new parameter through 12 plot modules produced a bug the unit tests could not
see: `plotsCorrelation._plot_corr_matrix` used `settings` while only its *caller* had been
given the parameter, so the whole graph run died with `NameError: name 'settings' is not
defined` -- but only for `-g all`, because that helper is reached through a multiprocessing
pool. A static check is cheaper and more reliable here than trying to exercise every plot
path, and it catches the same mistake for any parameter threaded through later.
"""
import ast
import pathlib

import pytest

MODULES = sorted(pathlib.Path('src/trnagraph/modules').glob('plots*.py'))
THREADED_PARAMS = ('settings', 'colormap', 'output')


def _out_of_scope_uses(path, name):
    tree = ast.parse(path.read_text())
    offenders = []
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        params = {a.arg for a in node.args.args} | {a.arg for a in node.args.kwonlyargs}
        if node.args.vararg:
            params.add(node.args.vararg.arg)
        if node.args.kwarg:
            params.add(node.args.kwarg.arg)
        if name in params:
            continue
        reads = any(isinstance(n, ast.Name) and n.id == name and isinstance(n.ctx, ast.Load)
                    for n in ast.walk(node))
        assigns = any(isinstance(n, ast.Name) and n.id == name and isinstance(n.ctx, ast.Store)
                      for n in ast.walk(node))
        # A method reaching self.<name> is fine; only a bare local read is the bug.
        if reads and not assigns:
            offenders.append(node.name)
    return offenders


@pytest.mark.parametrize('path', MODULES, ids=lambda p: p.name)
@pytest.mark.parametrize('name', THREADED_PARAMS)
def test_no_plot_helper_reads_a_parameter_it_was_not_given(path, name):
    offenders = _out_of_scope_uses(path, name)

    assert not offenders, (
        f"{path.name}: {offenders} read '{name}' without receiving it. Either add it to the "
        f"signature and thread it from the caller, or read it off self."
    )


def test_the_check_actually_detects_the_original_bug(tmp_path):
    """Guard the guard: a known-bad shape must fail, or this test proves nothing."""
    bad = tmp_path / 'plotsBad.py'
    bad.write_text('def helper(a):\n    return save(a, settings)\n')

    assert _out_of_scope_uses(bad, 'settings') == ['helper']
