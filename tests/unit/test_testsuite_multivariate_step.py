"""The demo suite exercises both multivariate graph types, not just the Venn.

The default `-g all` runs carry no --config, so the `multivariate` block that gates these two
is absent and both are correctly skipped. A second, gated invocation is the only thing that
reaches them -- and it asked for `-g venn` alone, which left `-g agreement` with no end-to-end
coverage at all: every test of it stubbed the DESeq2 boundary.
"""
import inspect

from trnagraph.modules import toolsTestSuite


def _graph_db_source():
    return inspect.getsource(toolsTestSuite.demoPipeline.graph_db)


def test_the_gated_invocation_asks_for_both_types():
    source = _graph_db_source()

    assert '-g venn' in source
    assert '-g agreement' in source


def test_both_ride_the_same_multivariate_config():
    """One block declares the grouping and thresholds for both, so a second config file would
    be two places to keep in step."""
    source = _graph_db_source()
    venn_line = [line for line in source.splitlines() if '-g venn' in line][0]

    assert 'agreement' in venn_line, 'one invocation, not two'
