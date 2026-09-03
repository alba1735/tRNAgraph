"""`-g agreement` reaching the plot module.

Same gate as `-g venn`: writing the `multivariate` block is the deliberate act, so once it is
there the type joins `-g all` rather than having to be named again. Without the block it is left
out and the run LOGS why -- a silently missing figure is indistinguishable from an empty result.
"""
from types import SimpleNamespace

from trnagraph.modules import adataGraph


def _status(multivariate=None, grp1='group', grp2='group'):
    suite = adataGraph.anndataGrapher.__new__(adataGraph.anndataGrapher)
    suite.args = SimpleNamespace(comparegrp1=grp1, comparegrp2=grp2)
    suite.config = SimpleNamespace(multivariate=multivariate)
    return suite._optional_graphtype_status()


def test_agreement_is_optional_not_part_of_the_fixed_all_list():
    assert 'agreement' in adataGraph.OPTIONAL_GRAPH_TYPES
    assert 'agreement' not in adataGraph.GRAPH_TYPES_ALL


def test_a_multivariate_block_lets_agreement_join_all():
    resolved, skipped = adataGraph.resolve_graphtypes(
        ['all'], _status(multivariate=SimpleNamespace()))

    assert 'agreement' in resolved
    assert not [name for name, _ in skipped if name == 'agreement']


def test_without_the_block_it_is_left_out_with_a_reason():
    resolved, skipped = adataGraph.resolve_graphtypes(['all'], _status())
    reasons = dict(skipped)

    assert 'agreement' not in resolved
    assert 'multivariate' in reasons['agreement']


def test_naming_it_explicitly_always_includes_it():
    """So the gate downstream can refuse BY NAME and say what to add, rather than the user
    getting no figure and no explanation."""
    resolved, _ = adataGraph.resolve_graphtypes(['agreement'], _status())

    assert 'agreement' in resolved


def test_the_precompute_covers_the_multivariate_grouping():
    """Agreement groups by multivariate.grouping, which need not equal volgrp. Fitting inside
    the worker pool deadlocks, so the grouping has to reach the pre-pool precompute."""
    suite = adataGraph.anndataGrapher.__new__(adataGraph.anndataGrapher)
    suite.args = SimpleNamespace(heatgrp='group', volgrp='group', graphtypes=['agreement'])
    suite.config = SimpleNamespace(multivariate=SimpleNamespace(grouping='timepoint'))

    assert 'timepoint' in suite._log2fc_groupings()


def test_the_precompute_does_not_add_a_grouping_no_plot_asked_for():
    suite = adataGraph.anndataGrapher.__new__(adataGraph.anndataGrapher)
    suite.args = SimpleNamespace(heatgrp='group', volgrp='group', graphtypes=['volcano'])
    suite.config = SimpleNamespace(multivariate=SimpleNamespace(grouping='timepoint'))

    assert 'timepoint' not in suite._log2fc_groupings()


def test_agreement_reads_the_grouping_column_s_colours():
    """Its points are levels of `multivariate.grouping`, so a style file naming Day 0 blue
    reaches it. Without this entry colormap resolved to None, the derivation never ran, and the
    figure silently used the built-in ramps instead -- found by rendering real data."""
    suite = adataGraph.anndataGrapher.__new__(adataGraph.anndataGrapher)
    suite.args = SimpleNamespace(covgrp=None, comparegrp1='group', comparegrp2='group')
    suite.config = SimpleNamespace(multivariate=SimpleNamespace(grouping='timepoint'))
    suite.cmap_dict = {'venn': 'venn'}

    assert suite._colormap_key('agreement') == 'timepoint'


def test_the_key_falls_back_when_there_is_no_block():
    suite = adataGraph.anndataGrapher.__new__(adataGraph.anndataGrapher)
    suite.args = SimpleNamespace(covgrp=None, comparegrp1='group', comparegrp2='group')
    suite.config = SimpleNamespace(multivariate=None)
    suite.cmap_dict = {'venn': 'venn'}

    assert suite._colormap_key('agreement') is None
