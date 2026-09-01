"""`trnagraph tools info` reports an AnnData object's own vocabulary.

Nothing else in the CLI answers "what can I type after --covgrp". Discovering the valid
values meant opening the .h5ad in Python, which is a large part of why a typo'd label
mattered so much -- the strict validator added alongside this command points users here by
name, so what it reports has to actually be the vocabulary the validator checks against.
"""
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd

import pytest

from trnagraph.modules import toolsInfo, toolsTG


def _make_adata(n_obs=6):
    obs = pd.DataFrame(
        {
            "sample": [f"s{i}" for i in range(n_obs)],
            "group": ["A", "B", "A", "B", "A", "B"],
            "amino": ["Ala"] * n_obs,
            "nreads_total_norm": np.zeros(n_obs),
        },
        index=[f"obs{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(
        {"coverage": ["uniquecoverage", "readstarts"], "gap": [False, False]},
        index=["v0", "v1"],
    )
    return ad.AnnData(X=np.zeros((n_obs, 2)), obs=obs, var=var)


def _inspector(**overrides):
    args = dict(anndata=None, column=None, json=False)
    args.update(overrides)
    return toolsInfo.AnnDataInspector(SimpleNamespace(**args))


def test_summary_reports_each_obs_column_with_its_unique_values():
    adata = _make_adata()
    summary = _inspector().summary(adata)

    columns = {c["name"]: c for c in summary["obs"]}
    assert set(columns) == {"sample", "group", "amino", "nreads_total_norm"}
    assert columns["group"]["n_unique"] == 2
    assert columns["group"]["values"] == ["A", "B"]
    assert columns["amino"]["n_unique"] == 1


def _make_wide_adata(n_obs=50):
    obs = pd.DataFrame(
        {
            "sample": [f"s{i}" for i in range(n_obs)],
            "trna": [f"tRNA-Ala-AGC-{i}" for i in range(n_obs)],
            "group": ["A" if i % 2 else "B" for i in range(n_obs)],
        },
        index=[f"obs{i}" for i in range(n_obs)],
    )
    return ad.AnnData(X=np.zeros((n_obs, 1)),
                      obs=obs,
                      var=pd.DataFrame({"coverage": ["coverage"]}, index=["v0"]))


def test_a_high_cardinality_column_is_truncated_but_still_counted():
    """'trna' is 436 values on a human build. Printing it in full would bury every column a
    user actually types, so the values are capped -- but the COUNT must stay honest, since
    that is what tells the reader the list is not the whole story."""
    summary = _inspector().summary(_make_wide_adata())
    trna = {c["name"]: c for c in summary["obs"]}["trna"]

    assert trna["n_unique"] == 50
    assert len(trna["values"]) == toolsInfo.VALUE_CAP
    assert trna["truncated"] is True


def test_column_prints_one_column_in_full_and_nothing_else():
    """The cap needs an escape hatch, and the escape hatch is per column: lifting it globally
    would reproduce exactly the unreadable output the cap exists to prevent."""
    summary = _inspector(column="trna").summary(_make_wide_adata())

    assert [c["name"] for c in summary["obs"]] == ["trna"]
    trna = summary["obs"][0]
    assert len(trna["values"]) == 50
    assert trna["truncated"] is False


def test_summary_covers_var_uns_layers_obsm_and_the_matrix_shape():
    """A grouping column is only half the vocabulary: --covtype comes from var, the log2FC
    cache and run provenance live in uns, and the normalization layers decide which
    --variant strings resolve."""
    adata = _make_adata()
    adata.layers["raw"] = np.zeros((6, 2))
    adata.obsm["sample_cluster_umap"] = np.zeros((6, 2))
    adata.uns["runinfo"] = {"version": "9.9.9"}

    summary = _inspector().summary(adata)

    assert summary["shape"] == (6, 2)
    assert {c["name"] for c in summary["var"]} == {"coverage", "gap"}
    assert {c["name"] for c in summary["var"] if c["name"] == "coverage"}
    assert summary["layers"] == ["raw"]
    assert summary["obsm"] == ["sample_cluster_umap"]
    assert "runinfo" in {entry["name"] for entry in summary["uns"]}


def test_reported_variants_are_only_those_the_object_can_actually_resolve():
    """--variant strings are the one vocabulary a user cannot guess: 'vst:full' resolves only
    if a vst layer was built, and a split tag only if that split was added. Listing the four
    normalizations unconditionally would advertise variants that abort on use."""
    adata = _make_adata()
    adata.layers["raw"] = np.zeros((6, 2))

    assert _inspector().summary(adata)["variants"] == ["norm:full", "raw:full"]

    adata.layers["vst"] = np.zeros((6, 2))
    adata.layers["norm_u60"] = np.zeros((6, 2))
    adata.uns["size_splits"] = {"u60": {}}
    variants = _inspector().summary(adata)["variants"]
    assert "vst:full" in variants
    assert "norm:u60" in variants
    assert "allfeatures:full" not in variants


def test_an_unknown_column_is_rejected_with_near_matches():
    """Silently reporting nothing for a mistyped --column would be the same silent-wrong-answer
    failure the strict validator was added to remove -- and this is the command that failure
    points users at, so it cannot exhibit it itself. --column accepts an obs OR a var column,
    so the suggestions have to be drawn from both."""
    with pytest.raises(toolsTG.UnknownLabelError) as raised:
        _inspector(column="grup").summary(_make_adata())
    assert "group" in str(raised.value)


def test_column_also_accepts_a_var_column():
    summary = _inspector(column="coverage").summary(_make_adata())
    assert [c["name"] for c in summary["var"]] == ["coverage"]
    assert summary["obs"] == []


def test_run_renders_a_human_report_that_names_columns_and_their_values(tmp_path):
    """The default view is read by a person deciding what to type next, so the column name and
    the values have to appear together on the page."""
    path = tmp_path / "obj.h5ad"
    _make_adata().write_h5ad(path)

    text = _inspector(anndata=str(path)).run()

    assert "group" in text
    assert "A" in text and "B" in text
    assert "uniquecoverage" in text, "var values are as much of the vocabulary as obs values"
    assert "norm:full" in text, "the --variant strings are the vocabulary a user cannot guess"


def test_json_output_carries_the_same_structure_as_the_summary(tmp_path):
    """--json exists so a wrapper can enumerate valid groupings. It is a view of the same
    summary the text renders, not a second implementation that could disagree with it."""
    import json as json_module

    path = tmp_path / "obj.h5ad"
    adata = _make_adata()
    adata.write_h5ad(path)

    payload = json_module.loads(_inspector(anndata=str(path), json=True).run())

    assert set(payload) == set(_inspector().summary(adata))
    assert {c["name"] for c in payload["obs"]} == {"sample", "group", "amino", "nreads_total_norm"}
    assert payload["shape"] == [6, 2]


def test_the_cli_command_prints_the_report(tmp_path):
    from typer.testing import CliRunner

    from trnagraph import cli

    path = tmp_path / "obj.h5ad"
    _make_adata().write_h5ad(path)

    result = CliRunner().invoke(cli.app, [
        "--skip-env-check", "--skip-update-check", "tools", "info", "-i", str(path),
    ])

    assert result.exit_code == 0
    assert "group" in result.output
    assert "norm:full" in result.output


def test_the_cli_command_rejects_an_unknown_column_without_a_traceback(tmp_path):
    from typer.testing import CliRunner

    from trnagraph import cli

    path = tmp_path / "obj.h5ad"
    _make_adata().write_h5ad(path)

    result = CliRunner().invoke(cli.app, [
        "--skip-env-check", "--skip-update-check", "tools", "info",
        "-i", str(path), "--column", "grup",
    ])

    assert result.exit_code == 1
    assert "group" in result.output
    assert "Traceback" not in result.output


def test_a_continuous_numeric_column_reports_a_range_instead_of_its_values():
    """Found on the real vibrChol1 object: obs carries 30 float count columns, and enumerating
    them printed hundreds of values like 0.8767311244266289 -- burying 'group' and 'amino',
    which are the columns anyone actually types. A count is not a label, so its range is the
    only summary of it worth printing."""
    n = 60
    obs = pd.DataFrame(
        {
            "group": ["A" if i % 2 else "B" for i in range(n)],
            "timepoint": [i % 3 for i in range(n)],
            "nreads_total_norm": np.arange(n, dtype=float) + 0.5,
        },
        index=[f"obs{i}" for i in range(n)],
    )
    adata = ad.AnnData(X=np.zeros((n, 1)), obs=obs,
                       var=pd.DataFrame({"coverage": ["coverage"]}, index=["v0"]))

    columns = {c["name"]: c for c in _inspector().summary(adata)["obs"]}

    counts = columns["nreads_total_norm"]
    assert counts["values"] == []
    assert counts["range"] == [0.5, 59.5]

    # A LOW-cardinality numeric column is a label, not a measurement: an ordered timepoint or
    # dose is exactly the sort of thing --volgrp is pointed at, and its values must stay
    # visible.
    assert columns["timepoint"]["values"] == ["0", "1", "2"]
    assert "range" not in columns["timepoint"]


def test_a_very_long_value_is_elided_unless_the_column_was_asked_for_by_name():
    """Also found on the real object: obs['refseq'] holds 75-90nt sequences, so twenty of them
    on one line ran to ~2000 characters and wrapped into an unreadable block. The cap on the
    NUMBER of values does nothing about the width of each one."""
    n = 4
    sequences = [f"GCCGACTTAGCTCAGTAGGTAGAGCAACTGACTTGTAATCAGTAGGCACCAGTTCGATTCCGGTAGTCGGCACC{i}"
                 for i in range(n)]
    obs = pd.DataFrame({"refseq": sequences, "group": ["A", "B", "A", "B"]},
                       index=[f"obs{i}" for i in range(n)])
    adata = ad.AnnData(X=np.zeros((n, 1)), obs=obs,
                       var=pd.DataFrame({"coverage": ["coverage"]}, index=["v0"]))

    elided = {c["name"]: c for c in _inspector().summary(adata)["obs"]}["refseq"]
    assert all(len(v) <= toolsInfo.VALUE_WIDTH for v in elided["values"])
    assert any(v.endswith("...") for v in elided["values"])

    full = _inspector(column="refseq").summary(adata)["obs"][0]
    assert sorted(full["values"]) == sorted(sequences)
