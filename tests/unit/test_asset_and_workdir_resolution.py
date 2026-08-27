"""Regression tests for roadmap.md's "`pip install .` (non-editable) fails `tools test`" item.

`toolsTestSuite.demoPipeline.__init__` used to guess a repo root by walking four directories up
from `__file__` and then look for assets at `<repo_root>/src/trnagraph/assets`. That path only
exists in a src-layout source checkout: from a wheel install the same walk lands on
`<venv>/lib/pythonX.Y`, so the asset copy fails outright -- and, worse, the suite first
`os.chdir`s into `<venv>/lib/pythonX.Y/test_vibrChol1` and, under `--all`, recursively deletes
its contents, i.e. inside the Python installation itself.

The assets were never the problem: they are packaged correctly (`[tool.setuptools.package-data]`
ships `trnagraph/assets/**` into the wheel). Only the lookup was wrong. These tests pin the two
invariants that make the lookup install-layout independent:

  1. assets resolve *inside the installed package*, never via a computed repo root, and
  2. the default working directory derives from the invocation cwd, never from `__file__`.
"""
import os
from pathlib import Path

import trnagraph
from trnagraph.modules import toolsTG


def test_assets_dir_resolves_inside_the_installed_package():
    """The invariant that a wheel install breaks: the assets directory is a *package*
    subdirectory, so it must be located relative to the imported package -- not relative to a
    source-tree layout that only exists under `pip install -e .`."""
    resolved = Path(toolsTG.assets_dir())
    assert resolved == Path(trnagraph.__file__).parent / "assets"


def test_assets_dir_contains_the_files_its_consumers_ask_for():
    """Behavioral counterpart to the invariant above: the demo pipeline copies `*.txt`/`*.json`
    out of this directory and makedb reads its covariance models from `cm/`, so a path that
    merely exists is not enough."""
    resolved = Path(toolsTG.assets_dir())
    for relative in ("vibrChol1.manifest.txt", "style.json", "cm/TRNAinf-bact.cm"):
        assert (resolved / relative).is_file(), f"missing packaged asset: {relative}"


def test_assets_dir_returns_a_plain_filesystem_path():
    """`toolsTestSuite` interpolates this into a shell `cp {assets_dir}/*.txt` command, so it
    must be an ordinary path string, not an importlib.resources Traversable."""
    assert isinstance(toolsTG.assets_dir(), str)


def test_resolve_work_dir_honours_an_explicit_directory(tmp_path, monkeypatch):
    """`-d`/`--directory` wins outright, and a relative value is anchored to the invocation
    cwd rather than left relative (the suite chdir()s into the result)."""
    from trnagraph.modules import toolsTestSuite

    monkeypatch.chdir(tmp_path)
    assert Path(toolsTestSuite._resolve_work_dir(str(tmp_path / "elsewhere"))) == tmp_path / "elsewhere"
    assert Path(toolsTestSuite._resolve_work_dir("elsewhere")) == Path(os.getcwd()) / "elsewhere"


def test_resolve_work_dir_defaults_under_the_invocation_cwd(tmp_path, monkeypatch):
    """Omitting `-d` puts the workspace where the command was actually run, matching what
    docs/testSuite.md already advertises."""
    from trnagraph.modules import toolsTestSuite

    monkeypatch.chdir(tmp_path)
    assert Path(toolsTestSuite._resolve_work_dir(None)) == Path(os.getcwd()) / "test_vibrChol1"


def test_resolve_work_dir_never_lands_inside_the_installed_package(tmp_path, monkeypatch):
    """The regression itself. The old default was `<four dirs above toolsTestSuite.py>/
    test_vibrChol1`, which under a wheel install resolves inside the Python installation --
    a directory `--all` then recursively deletes. Nothing about the default may depend on
    where the package happens to be installed."""
    from trnagraph.modules import toolsTestSuite

    monkeypatch.chdir(tmp_path)
    resolved = Path(toolsTestSuite._resolve_work_dir(None))
    package_dir = Path(trnagraph.__file__).parent
    assert package_dir not in resolved.parents
    assert package_dir.parent not in resolved.parents
