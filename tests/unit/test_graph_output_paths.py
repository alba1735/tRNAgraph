"""Where each graph type's files land.

The rule the path builder follows: a directory segment is added only when nothing already in
the path distinguishes that output. `--variant` always qualifies, because the same plot of a
different normalization really is a different picture. The read basis qualifies most types --
but not the ones whose filenames already carry it, which is why coverage has no `allreads/`
segment (`--covtype` names the category) and why PCA should not have one either.
"""
from types import SimpleNamespace

import pytest

from trnagraph.modules import adataGraph, toolsTG


def _grapher(read_basis, variant='norm:full', output='figures'):
    grapher = object.__new__(adataGraph.anndataGrapher)
    grapher.args = SimpleNamespace(output=output)
    grapher.variant_spec = SimpleNamespace(raw=variant)
    grapher.read_basis = read_basis
    return grapher


@pytest.mark.parametrize('graph_type', ['volcano', 'heatmap', 'correlation', 'count'])
def test_the_read_basis_qualifies_an_ordinary_graph_type(graph_type):
    unique = _grapher(toolsTG.READ_BASIS_UNIQUE)
    all_reads = _grapher(toolsTG.READ_BASIS_ALL)

    assert unique._output_dir_for(graph_type) == f'figures/{graph_type}/'
    assert all_reads._output_dir_for(graph_type) == f'figures/{graph_type}/{toolsTG.READ_BASIS_ALL}/'


@pytest.mark.parametrize('graph_type', ['coverage', 'pca'])
def test_types_whose_filenames_already_name_the_basis_get_no_segment(graph_type):
    """
    coverage was already exempt because --covtype names the category. PCA is exempt for the
    same reason by a different route: every PCA filename carries _readtype_label() ('total'
    against 'total_unique'), so the basis is in the filename already.

    Without this, a default run and an --allreads run wrote exactly the same two files -- one
    into pca/, one into pca/allreads/ -- because OVERVIEW_TRNA_READTYPES pins the both-bases
    comparison so it survives whatever --pcareadtypes asks for.
    """
    unique = _grapher(toolsTG.READ_BASIS_UNIQUE)
    all_reads = _grapher(toolsTG.READ_BASIS_ALL)

    assert unique._output_dir_for(graph_type) == f'figures/{graph_type}/'
    assert all_reads._output_dir_for(graph_type) == f'figures/{graph_type}/'


def test_a_non_default_variant_always_qualifies_every_type():
    """Unlike the basis, a variant is never recoverable from a filename, so it always adds a
    segment -- including for the two types exempt from the basis segment."""
    grapher = _grapher(toolsTG.READ_BASIS_UNIQUE, variant='norm:u60')

    assert grapher._output_dir_for('pca') == 'figures/pca/norm_u60/'
    assert grapher._output_dir_for('coverage') == 'figures/coverage/norm_u60/'
