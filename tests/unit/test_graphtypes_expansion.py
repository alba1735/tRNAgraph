"""`-g all` unions with what else was asked for, and picks up opt-in types that are ready.

Two behaviours, one function.

The first is a defect: the expansion REPLACED the requested list, so `-g all -g venn` produced
exactly `all` and the venn was silently dropped. `-g all -g compare` had the same problem. A user
asking for two things and receiving one, with nothing said, is the worst shape a flag can have.

The second is the design change it enables. `compare` and `venn` are excluded from `all` because
they need choices nobody makes by typing `-g all` -- but when those choices HAVE been made, in a
config block or in --comparegrp1/2, requiring the type to be named again is ceremony rather than
a safeguard. So `all` now includes them when their prerequisites are satisfied.

What it does not do is skip silently. A type left out of `all` reports why, because "no Venn
appeared" is otherwise indistinguishable from "there was nothing to draw" -- the same confusion
the per-diagram skip messages exist to prevent.

Naming a type EXPLICITLY always includes it, prerequisites or not: the gate downstream then
refuses by name, which tells the user what to fix. Quietly dropping an explicit request would
not.
"""

from trnagraph.modules import adataGraph


def test_all_no_longer_swallows_an_explicitly_named_type():
    """The defect: `-g all -g venn` used to produce `all` alone."""
    types, _ = adataGraph.resolve_graphtypes(['all', 'venn'], {'venn': None, 'compare': 'no groups'})

    assert 'venn' in types
    assert 'volcano' in types, 'and still everything all covers'


def test_a_ready_optional_type_joins_all():
    types, skipped = adataGraph.resolve_graphtypes(
        ['all'], {'venn': None, 'agreement': None, 'compare': 'no groups'})

    assert 'venn' in types
    assert 'agreement' in types
    assert skipped == [('compare', 'no groups')]


def test_an_unready_optional_type_is_left_out_with_its_reason():
    types, skipped = adataGraph.resolve_graphtypes(
        ['all'], {'venn': 'no multivariate block', 'compare': 'no groups'})

    assert 'venn' not in types
    assert dict(skipped)['venn'] == 'no multivariate block'


def test_naming_a_type_explicitly_includes_it_even_when_unready():
    """So the downstream gate can refuse BY NAME, which tells the user what to add."""
    types, skipped = adataGraph.resolve_graphtypes(['venn'], {'venn': 'no multivariate block'})

    assert types == ['venn']
    assert skipped == [], 'an explicit request is never reported as skipped'


def test_asking_for_one_type_does_not_pull_in_the_others():
    types, _ = adataGraph.resolve_graphtypes(['volcano'], {'venn': None, 'compare': None})

    assert types == ['volcano']


def test_types_are_not_duplicated():
    types, _ = adataGraph.resolve_graphtypes(['all', 'volcano', 'venn'], {'venn': None})

    assert len(types) == len(set(types))


def test_all_still_expands_to_the_descriptive_set_when_nothing_is_ready():
    types, _ = adataGraph.resolve_graphtypes(['all'], {'venn': 'x', 'compare': 'y'})

    assert set(types) == set(adataGraph.GRAPH_TYPES_ALL)
