"""Unknown obs/var labels abort the run instead of being silently substituted.

`resolve_grp_column` used to warn and fall back to 'sample', so a typo'd `--covgrp` produced
a complete, plausible set of figures grouped by the wrong column -- a run that looks healthy
and whose output means nothing. It also only covered five of the twelve grouping parameters;
the rest reached pandas and failed with whatever pandas raised.

These tests pin the replacement: one shared validator, run up front, reporting every bad
label at once with its near matches, raising a dedicated error the CLI can render as a usage
message rather than a traceback.
"""
import logging
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from trnagraph.modules import adataGraph, toolsTG


def _make_adata(n_obs=4):
    obs = pd.DataFrame(
        {
            "sample": [f"s{i}" for i in range(n_obs)],
            "group": ["A" if i % 2 == 0 else "B" for i in range(n_obs)],
            "amino": ["Ala"] * n_obs,
            "trna": [f"tRNA-Ala-AGC-{i}" for i in range(n_obs)],
            "nreads_total_unique_norm": np.zeros(n_obs),
            "nreads_total_norm": np.zeros(n_obs),
            "nreads_fiveprime_norm": np.zeros(n_obs),
            "nreads_wholecounts_unique_norm": np.zeros(n_obs),
        },
        index=[f"obs{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(
        {"coverage": ["uniquecoverage", "readstarts"], "gap": [False, False]},
        index=["v0", "v1"],
    )
    return ad.AnnData(X=np.zeros((n_obs, 2)), obs=obs, var=var)


def test_an_unknown_obs_column_aborts_naming_the_parameter_and_a_near_match():
    """The whole point of the change: a typo must stop the run, and say what to type instead."""
    adata = _make_adata()
    with pytest.raises(toolsTG.UnknownLabelError) as raised:
        toolsTG.validate_labels(adata, [("covgrp", "goup", "obs")])
    message = str(raised.value)
    assert "--covgrp" in message
    assert "goup" in message
    assert "group" in message


def test_a_coverage_type_is_checked_against_var_not_obs():
    """--covtype names an adata.var['coverage'] value; checking it against obs columns
    would reject every valid coverage type and suggest nonsense for the invalid ones."""
    adata = _make_adata()
    toolsTG.validate_labels(adata, [("covtype", "readstarts", "coverage")])

    with pytest.raises(toolsTG.UnknownLabelError) as raised:
        toolsTG.validate_labels(adata, [("covtype", "readstart", "coverage")])
    assert "readstarts" in str(raised.value)


def test_a_readtype_is_checked_against_the_objects_own_count_columns():
    """--diffrts/--pcareadtypes take BARE readtypes, which are neither obs columns nor var
    values -- they name the 'nreads_<rt>_norm' family the object actually carries, which is
    the same thing resolve_readtype() consults when deciding whether a unique basis exists."""
    adata = _make_adata()
    toolsTG.validate_labels(adata, [("diffrts", "total", "readtype")])
    toolsTG.validate_labels(adata, [("diffrts", "wholecounts", "readtype")])

    with pytest.raises(toolsTG.UnknownLabelError) as raised:
        toolsTG.validate_labels(adata, [("diffrts", "fiveprim", "readtype")])
    assert "fiveprime" in str(raised.value)


def test_every_bad_label_is_reported_in_one_error_across_domains():
    """A graph command carries a dozen label-valued options. Aborting on the first bad one
    turns fixing a command line into as many round trips as there are typos."""
    adata = _make_adata()
    with pytest.raises(toolsTG.UnknownLabelError) as raised:
        toolsTG.validate_labels(adata, [
            ("covgrp", "goup", "obs"),
            ("volgrp", "group", "obs"),
            ("covtype", "readstart", "coverage"),
            ("diffrts", "fiveprim", "readtype"),
        ])
    message = str(raised.value)
    assert message.startswith("3 unknown labels")
    for expected in ("--covgrp", "--covtype", "--diffrts"):
        assert expected in message
    assert "--volgrp" not in message, "a valid label must not be reported as a problem"


def _make_grapher(args, adata):
    """A grapher without __init__, which would need a real .h5ad on disk."""
    grapher = object.__new__(adataGraph.anndataGrapher)
    grapher.args = args
    grapher.adata = adata
    grapher.logger = logging.getLogger("test")
    return grapher


def _valid_args(**overrides):
    args = dict(
        covgrp="group", covobs="trna", comparegrp1="group", comparegrp2="amino",
        heatgrp="group", volgrp="group", clustergrp="amino", clusterlabels=None,
        corrgroup="sample", logogrp="amino", pcamarkers="sample", pcacolors="group",
        radargrp="group", covtype="uniquecoverage",
        diffrts=["total"], pcareadtypes=["total"], graphtypes=["compare"], corrmask="none", heatorient="vertical",
    )
    args.update(overrides)
    return SimpleNamespace(**args)


def test_the_grapher_accepts_a_fully_valid_set_of_label_parameters():
    grapher = _make_grapher(_valid_args(), _make_adata())
    grapher._validate_label_args()


def test_the_grapher_reports_every_bad_label_parameter_at_once():
    """--clustergrp, --logogrp, --corrgroup, --pcamarkers, --pcacolors, --radargrp and
    --clusterlabels went through no validation at all and failed at first use inside pandas.
    All of them are checked here, in one pass, before any plotting begins."""
    adata = _make_adata()
    args = _valid_args(covgrp="goup", clustergrp="amin", pcacolors="grp")
    grapher = _make_grapher(args, adata)

    with pytest.raises(toolsTG.UnknownLabelError) as raised:
        grapher._validate_label_args()
    message = str(raised.value)
    for expected in ("--covgrp", "--clustergrp", "--pcacolors"):
        assert expected in message
    assert args.covgrp == "goup", "a rejected label must not be silently rewritten"


def test_an_unset_optional_label_parameter_is_not_validated():
    """--clusterlabels defaults to None, which means 'no labels', not a column named None."""
    grapher = _make_grapher(_valid_args(clusterlabels=None), _make_adata())
    grapher._validate_label_args()


def test_the_per_module_resolver_no_longer_substitutes_a_fallback_column():
    """plotsCoverage/plotsHeatmap/plotsCompare each resolve their own grouping argument so
    they stay usable standalone, and the Python API reaches them without going through
    anndataGrapher at all. If that path still fell back to 'sample', the silent substitution
    would simply move rather than be removed."""
    adata = _make_adata()
    assert toolsTG.resolve_grp_column(adata, "group", "covgrp") == "group"

    with pytest.raises(toolsTG.UnknownLabelError) as raised:
        toolsTG.resolve_grp_column(adata, "goup", "covgrp")
    assert "group" in str(raised.value)


def test_labels_are_validated_before_the_log2fc_precompute():
    """anndataGrapher.__init__ runs a PyDESeq2 precompute over --diffrts before any plotting.
    Validating after it would make a user wait through the most expensive part of the command
    to be told they mistyped a flag -- and the precompute reads heatgrp/volgrp itself, so a
    bad label would surface there first, as a pandas error naming neither the flag nor the
    alternatives."""
    import ast
    import inspect

    tree = ast.parse(inspect.getsource(adataGraph.anndataGrapher.__init__).lstrip())
    order = [(node.lineno, node.func.attr) for node in ast.walk(tree)
             if isinstance(node, ast.Call) and getattr(node.func, "attr", None)
             in ("_validate_label_args", "_precompute_and_persist_log2fc")]
    names = [name for _, name in sorted(order)]

    assert names == ["_validate_label_args", "_precompute_and_persist_log2fc"], (
        f"expected validation before the precompute, got {names}")


def test_a_weak_match_is_not_offered_as_a_suggestion():
    """Real case from vibrChol1: --covgrp 'treatement' against an object with no 'treatment'
    column produced "did you mean: fragment?" -- a confident-sounding suggestion for an
    unrelated column, which is worse than none. Every genuine typo measured on real labels
    scores >= 0.75 and every misleading match <= 0.667, so the threshold sits between them
    (and matches the --style colormap validator's own)."""
    adata = _make_adata()
    adata.obs["fragment"] = "whole"

    with pytest.raises(toolsTG.UnknownLabelError) as raised:
        toolsTG.validate_labels(adata, [("covgrp", "treatement", "obs")])
    assert "did you mean" not in str(raised.value)


def test_a_label_with_no_plausible_match_lists_what_is_available_instead():
    """Without a suggestion the message would be 'that is wrong' with no way forward, so the
    fallback is the vocabulary itself."""
    adata = _make_adata()
    with pytest.raises(toolsTG.UnknownLabelError) as raised:
        toolsTG.validate_labels(adata, [("covgrp", "zzzzzzzzzz", "obs")])
    message = str(raised.value)
    assert "group" in message and "amino" in message


def test_unknown_labels_and_invalid_parameters_share_a_usage_error_base():
    """Both are usage mistakes the CLI renders as a message rather than a traceback, so the
    guard catches one base type. They stay distinct types because they are distinct problems:
    an unknown label names something absent, an invalid parameter names things that exist in
    a combination that cannot work."""
    assert issubclass(toolsTG.UnknownLabelError, toolsTG.UsageError)
    assert issubclass(toolsTG.InvalidParameterError, toolsTG.UsageError)


def test_the_two_compare_columns_must_name_different_columns():
    """--comparegrp1 and --comparegrp2 both default to 'group'. log2fc_compare_df then pivots
    on a duplicated column and pandas raises `Grouper for 'group' not 1-dimensional`, which
    names neither flag. The fold change is taken BETWEEN comparegrp2 values WITHIN each
    comparegrp1 value, so the two naming one column is not a comparison at all."""
    grapher = _make_grapher(_valid_args(comparegrp1="group", comparegrp2="group"), _make_adata())

    with pytest.raises(toolsTG.InvalidParameterError) as raised:
        grapher._validate_label_args()
    message = str(raised.value)
    assert "--comparegrp1" in message and "--comparegrp2" in message
    assert "group" in message


def test_a_typo_and_the_compare_collision_are_reported_in_one_error():
    """Both are found in the same up-front pass, so fixing a command line takes one round trip
    rather than one per problem."""
    args = _valid_args(comparegrp1="group", comparegrp2="group", volgrp="grp")
    grapher = _make_grapher(args, _make_adata())

    with pytest.raises(toolsTG.UsageError) as raised:
        grapher._validate_label_args()
    message = str(raised.value)
    assert "2 problems" in message
    assert "--volgrp" in message
    assert "--comparegrp1" in message


def test_the_compare_collision_is_only_a_problem_when_compare_was_requested():
    """--comparegrp1 and --comparegrp2 BOTH DEFAULT to 'group', so treating the collision as
    unconditional would abort every ordinary run of a command that never asked for a compare
    plot. compare is excluded from `-g all` by design, so it is present in graphtypes only
    when named explicitly -- which is exactly when the collision matters."""
    adata = _make_adata()

    ordinary = _make_grapher(_valid_args(comparegrp1="group", comparegrp2="group",
                                         graphtypes=["all"]), adata)
    ordinary._validate_label_args()

    requested = _make_grapher(_valid_args(comparegrp1="group", comparegrp2="group",
                                          graphtypes=["compare", "pca"]), adata)
    with pytest.raises(toolsTG.InvalidParameterError):
        requested._validate_label_args()


def test_an_unrecognised_corrmask_is_reported_up_front_with_everything_else():
    """--corrmask is checked in the same pass as the labels rather than when the correlation
    module first draws: a run that plots six graph types before rejecting a typo has already
    spent the expensive part of the command."""
    grapher = _make_grapher(_valid_args(corrmask="top", volgrp="grp"), _make_adata())

    with pytest.raises(toolsTG.InvalidParameterError) as raised:
        grapher._validate_label_args()
    message = str(raised.value)
    assert "--corrmask" in message
    assert "upper" in message and "lower" in message
    assert "--volgrp" in message, "reported in the same batch as the label problems"


def test_a_valid_corrmask_passes():
    for value in ("none", "upper", "lower"):
        _make_grapher(_valid_args(corrmask=value, graphtypes=["all"]), _make_adata())._validate_label_args()


def test_an_unrecognised_heatorient_is_reported_up_front():
    """Same reasoning as --corrmask: heatmap is one of ten graph types, and a run should not
    reach it before rejecting a typo."""
    grapher = _make_grapher(_valid_args(heatorient="sideways", graphtypes=["all"]), _make_adata())

    with pytest.raises(toolsTG.InvalidParameterError) as raised:
        grapher._validate_label_args()
    message = str(raised.value)
    assert "--heatorient" in message
    assert "vertical" in message and "horizontal" in message
