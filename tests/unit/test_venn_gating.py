"""`-g venn` is opt-in, and says so when it is asked for without being configured.

The multivariate plots are deliberate analyses over a specific experimental design, not
descriptive summaries of whatever was built. Producing them automatically would hand a user
figures whose sets and thresholds nobody chose, which is precisely how a wrong conclusion gets
drawn from a right pipeline. So they stay out of `-g all` and require a config block.

`-g compare` set this precedent and relies on the same mechanism -- it appears in graphtypes
only when named explicitly, which adataGraph depends on when checking its group arguments.
"""
import pytest

from trnagraph.modules import adataGraph, plotsVenn
from trnagraph.modules.toolsSchemas import RunConfig


def test_venn_is_not_part_of_all():
    assert 'venn' not in adataGraph.GRAPH_TYPES_ALL


def test_compare_is_still_not_part_of_all():
    """The existing precedent, pinned alongside so a future edit cannot quietly fold either in."""
    assert 'compare' not in adataGraph.GRAPH_TYPES_ALL


def test_all_still_expands_to_the_descriptive_plots():
    assert set(adataGraph.GRAPH_TYPES_ALL) == {
        'cluster', 'correlation', 'count', 'coverage', 'heatmap',
        'logo', 'mismatch', 'pca', 'radar', 'volcano'}


def test_asking_for_venn_without_a_config_block_names_what_is_missing():
    with pytest.raises(Exception) as excinfo:
        plotsVenn.require_multivariate_config(None)

    message = str(excinfo.value)
    assert 'multivariate' in message, 'names the block to add'
    assert '--config' in message, 'names the file it belongs in'


def test_a_config_without_the_block_is_also_refused():
    config = RunConfig.model_validate({'name': 'plain'})

    with pytest.raises(Exception) as excinfo:
        plotsVenn.require_multivariate_config(config)

    assert 'multivariate' in str(excinfo.value)


def test_a_configured_run_is_allowed_through():
    config = RunConfig.model_validate({
        'name': 'organoid', 'multivariate': {'grouping': 'group'}})

    assert plotsVenn.require_multivariate_config(config).grouping == 'group'
