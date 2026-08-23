"""Regression test for roadmap.md's --vollabels default change. Labeling every significant point
by default (the previous behavior) has an unbounded cost via adjustText's per-plot label
repulsion -- a real comparison on the actual hg38 dataset produced 70-77 significant points on a
single pair. Capping the default at 100 bounds the worst case while leaving the flag's other
behavior (0 disables all labels, any explicit N means exactly N) unchanged."""
import inspect

from trnagraph.modules.plotsVolcano import visualizer


def test_cli_vollabels_default_is_100():
    from trnagraph import cli
    param = inspect.signature(cli.graph).parameters['vollabels']
    assert param.default.default == 100


def test_visualizer_toplabels_default_is_100():
    param = inspect.signature(visualizer).parameters['toplabels']
    assert param.default == 100
