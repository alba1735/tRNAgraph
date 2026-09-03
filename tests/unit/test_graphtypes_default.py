"""The `-g` default is `all` alone, not `all` plus a list of everything.

The option's default used to enumerate every type alongside `all`. That was harmless only
because the expansion REPLACED the requested list, so the enumeration was discarded and `all`
decided everything. Once expansion started unioning -- so that `-g all -g venn` stops dropping
the venn -- the enumeration became a set of EXPLICIT requests, and an explicit request is
deliberately honoured even when its prerequisites are missing so the gate can refuse it by name.

The result was that every ordinary `graph` run, with no `-g` at all, demanded `venn` and aborted
unless a `multivariate` config block happened to exist. Found by the end-to-end suite, whose
split-variant step runs graph with no --config.
"""
from trnagraph import cli


def _default_graphtypes():
    return cli.graph.__defaults__[
        list(cli.graph.__code__.co_varnames).index('graphtypes')
        - (cli.graph.__code__.co_argcount - len(cli.graph.__defaults__))
    ].default


def test_the_default_is_all_alone():
    assert _default_graphtypes() == ['all']


def test_the_default_does_not_name_an_opt_in_type():
    """Naming one explicitly is a request, and a request is honoured rather than dropped."""
    from trnagraph.modules.adataGraph import OPTIONAL_GRAPH_TYPES

    assert not set(_default_graphtypes()) & set(OPTIONAL_GRAPH_TYPES)
