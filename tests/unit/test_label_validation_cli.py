"""A typo'd label is a usage error, and must read like one.

toolsTG.UnknownLabelError carries the whole message a user needs -- which flags were wrong,
what the near matches are, and how to list the real values. Letting it propagate would bury
that under a traceback of tRNAgraph's own call stack, which tells the user nothing they can
act on. cli.py renders it the way toolsTestSuite.WorkspaceNotOwnedError is already rendered:
one ERROR line, exit 1, no traceback.
"""
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

from trnagraph import cli

runner = CliRunner()


def _write_adata(path, n_obs=4):
    obs = pd.DataFrame(
        {
            "sample": [f"s{i}" for i in range(n_obs)],
            "group": ["A" if i % 2 == 0 else "B" for i in range(n_obs)],
            "amino": ["Ala"] * n_obs,
            "trna": [f"tRNA-Ala-AGC-{i}" for i in range(n_obs)],
            "nreads_total_unique_norm": np.zeros(n_obs),
            "nreads_total_norm": np.zeros(n_obs),
        },
        index=[f"obs{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(
        {"coverage": ["uniquecoverage", "coverage"], "gap": [False, False]},
        index=["v0", "v1"],
    )
    ad.AnnData(X=np.zeros((n_obs, 2)), obs=obs, var=var).write_h5ad(path)
    return path


def test_graph_reports_a_typod_grouping_column_without_a_traceback(tmp_path):
    h5ad = _write_adata(tmp_path / "obj.h5ad")
    result = runner.invoke(cli.app, [
        "--skip-env-check", "--skip-update-check", "graph",
        "-i", str(h5ad), "-o", str(tmp_path / "figures"),
        "-g", "pca", "--covgrp", "goup",
    ])

    assert result.exit_code == 1
    output = result.output
    assert "goup" in output
    assert "group" in output, "the near match is the actionable half of the message"
    assert "tools info" in output
    assert "Traceback" not in output


def test_cluster_reports_a_typod_coverage_type_without_a_traceback(tmp_path):
    """--coveragetype names adata.var['coverage'] values. Unvalidated, a typo produced an
    empty feature matrix and UMAP failed somewhere far from the cause."""
    h5ad = _write_adata(tmp_path / "obj.h5ad")
    result = runner.invoke(cli.app, [
        "--skip-env-check", "--skip-update-check", "analyze", "cluster",
        "-i", str(h5ad), "-v", "uniquecoverag",
    ])

    assert result.exit_code == 1
    assert "uniquecoverage" in result.output
    assert "Traceback" not in result.output


def test_log2fc_reports_a_typod_group_without_a_traceback(tmp_path):
    """`tools log2fc -g` is the same kind of label as `graph --volgrp`, and reached PyDESeq2
    unvalidated."""
    h5ad = _write_adata(tmp_path / "obj.h5ad")
    result = runner.invoke(cli.app, [
        "--skip-env-check", "--skip-update-check", "tools", "log2fc",
        "-i", str(h5ad), "-g", "grp",
    ])

    assert result.exit_code == 1
    assert "group" in result.output
    assert "Traceback" not in result.output


def test_log2fc_readtypes_are_validated_against_their_own_basis_carrying_vocabulary(tmp_path):
    """`tools log2fc -r` takes readtypes that DO carry the basis ('total_unique'), because it
    builds 'nreads_{rt}_norm' directly -- unlike `graph --diffrts`, where the basis is set once
    by --allreads and a basis-carrying readtype is rejected. Validating one against the other's
    vocabulary would reject every correct value."""
    h5ad = _write_adata(tmp_path / "obj.h5ad")
    ok = runner.invoke(cli.app, [
        "--skip-env-check", "--skip-update-check", "tools", "log2fc",
        "-i", str(h5ad), "-g", "group", "-r", "total_unique",
    ])
    assert "unknown label" not in ok.output

    bad = runner.invoke(cli.app, [
        "--skip-env-check", "--skip-update-check", "tools", "log2fc",
        "-i", str(h5ad), "-g", "group", "-r", "total_uniqu",
    ])
    assert bad.exit_code == 1
    assert "total_unique" in bad.output
    assert "Traceback" not in bad.output


def test_the_error_is_reported_once_not_twice(tmp_path):
    """handle_output records the failure in the run's log file and the guard prints it to the
    terminal. Both writing to the console printed the whole paragraph twice, which reads like
    two separate failures."""
    h5ad = _write_adata(tmp_path / "obj.h5ad")
    result = runner.invoke(cli.app, [
        "--skip-env-check", "--skip-update-check", "graph",
        "-i", str(h5ad), "-o", str(tmp_path / "figures"),
        "-g", "pca", "--covgrp", "goup",
    ])

    assert result.output.count("unknown label") == 1


def test_compare_at_default_settings_reports_the_collision_without_a_traceback(tmp_path):
    """`-g compare` with no grouping flags is the single most likely way to reach this plot,
    and it used to die with `ValueError: Grouper for 'group' not 1-dimensional` from inside
    pandas -- a message naming neither flag nor what to do about it."""
    h5ad = _write_adata(tmp_path / "obj.h5ad")
    result = runner.invoke(cli.app, [
        "--skip-env-check", "--skip-update-check", "graph",
        "-i", str(h5ad), "-o", str(tmp_path / "figures"), "-g", "compare",
    ])

    assert result.exit_code == 1
    assert "--comparegrp1" in result.output and "--comparegrp2" in result.output
    assert "Grouper" not in result.output
    assert "Traceback" not in result.output


def test_an_ordinary_run_is_unaffected_by_the_compare_defaults(tmp_path):
    """The regression this guards: --comparegrp1/--comparegrp2 both default to 'group', so a
    collision check that ignored which graph types were requested would abort every run."""
    h5ad = _write_adata(tmp_path / "obj.h5ad")
    result = runner.invoke(cli.app, [
        "--skip-env-check", "--skip-update-check", "graph",
        "-i", str(h5ad), "-o", str(tmp_path / "figures"), "-g", "count",
    ])

    assert "comparegrp" not in result.output
