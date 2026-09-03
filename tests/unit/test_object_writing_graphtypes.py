"""The graph types that write the .h5ad back must not run concurrently.

`venn` and `agreement` both persist their membership onto the object, and graph types are
dispatched through a multiprocessing.Pool. Two workers writing the same file at once fails with
`BlockingIOError: unable to lock file` part-way through a run -- which is what happened the
first time the demo suite asked for both in one invocation. It was latent while venn was the
only writer.

They join `coverage` in non_pooled_graphs, which already exists for types that cannot share the
pool, so they run sequentially in the parent after it drains.
"""
from trnagraph.modules import adataGraph


def test_both_object_writers_are_declared():
    assert set(adataGraph.OBJECT_WRITING_GRAPH_TYPES) == {'venn', 'agreement'}


def test_every_declared_writer_actually_writes_the_object():
    """A type in this list that does not write is a needless serialization; one that writes and
    is missing from it is the BlockingIOError."""
    import inspect
    source = inspect.getsource(adataGraph.anndataGrapher.dispatch_plot)

    for gt in adataGraph.OBJECT_WRITING_GRAPH_TYPES:
        block = source.split(f"if gt == '{gt}':")[1].split('if gt ==')[0]
        assert 'adata_original.write' in block, f'{gt} is listed but does not write'


def test_no_other_graph_type_writes_the_object():
    import inspect
    source = inspect.getsource(adataGraph.anndataGrapher.dispatch_plot)

    writers = {block.split("'")[0] for block in source.split("if gt == '")[1:]
               if 'adata_original.write' in block.split('if gt ==')[0]}
    assert writers == set(adataGraph.OBJECT_WRITING_GRAPH_TYPES)


def test_they_are_pulled_out_of_the_pool():
    suite = adataGraph.anndataGrapher.__new__(adataGraph.anndataGrapher)
    requested = ['pca', 'venn', 'coverage', 'agreement']

    pooled, non_pooled = suite._partition_pooled_graphtypes(requested)

    assert 'venn' not in pooled and 'agreement' not in pooled
    assert set(non_pooled) == {'coverage', 'venn', 'agreement'}
    assert pooled == ['pca']


def test_a_run_without_them_is_unaffected():
    suite = adataGraph.anndataGrapher.__new__(adataGraph.anndataGrapher)

    pooled, non_pooled = suite._partition_pooled_graphtypes(['pca', 'volcano'])

    assert pooled == ['pca', 'volcano']
    assert non_pooled == []
