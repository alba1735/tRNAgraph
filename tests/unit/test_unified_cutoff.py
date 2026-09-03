"""One read-count cutoff, not three.

`--volcutoff` and `--heatcutoff` were the same knob spelled twice -- the same pre-fit DESeq2
filter, the same default -- and adataGraph already deduplicated them before fitting. The
multivariate block carried a third, `presence_cutoff`, at a different default, which is how a
tRNA could sit in a Venn circle at 20 reads and be absent from every volcano beside it, having
never been fitted at 80.

They are now one `--cutoff`, default 20, governing the volcano, the heatmap, the agreement
volcano and Venn presence. Cluster's `--readcutoff` is deliberately untouched: it filters
observations before UMAP, which is a different question about a different axis.
"""
import inspect

import pytest

from trnagraph import cli
from trnagraph.modules import toolsSchemas


def test_the_graph_command_takes_one_cutoff():
    params = inspect.signature(cli.graph).parameters

    assert 'cutoff' in params
    assert 'volcutoff' not in params
    assert 'heatcutoff' not in params


def test_the_default_is_the_presence_cutoff_number():
    """20, not 80: the unification keeps Venn output as it was and widens the DE plots."""
    assert inspect.signature(cli.graph).parameters['cutoff'].default.default == 20


def test_cluster_keeps_its_own_read_cutoff():
    """A different operation -- it drops observations, not features -- so it is not folded in."""
    assert 'readcutoff' in inspect.signature(cli.cluster).parameters


def test_the_multivariate_block_no_longer_carries_a_cutoff():
    assert 'presence_cutoff' not in toolsSchemas.MultivariateConfig.model_fields


def test_a_config_naming_a_removed_cutoff_key_is_rejected():
    """extra='forbid' means an old config fails at load with a message, rather than silently
    running at a default the user did not choose."""
    with pytest.raises(Exception):
        toolsSchemas.MultivariateConfig(grouping='timepoint', presence_cutoff=20)


def test_the_graph_flags_block_exposes_the_unified_key():
    fields = toolsSchemas.GraphFlags.model_fields

    assert 'cutoff' in fields
    assert 'volcutoff' not in fields
    assert 'heatcutoff' not in fields
