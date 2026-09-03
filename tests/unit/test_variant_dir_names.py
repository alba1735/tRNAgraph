"""Where a build's results and graphs go, split or not.

A plain build writes straight to `results/`; a split build nests every variant, including the
unsplit one, under `results/complete/`, `results/u60/`, `results/o60/`. The rule is two lines in
variant_dir_names() and is shared by adataBuild and toolsSplit -- one writes what the other later
reads -- so it is pinned here rather than by running a whole second demo build to observe it.
"""
from types import SimpleNamespace

from trnagraph.modules import toolsTG


def _args(readlengthsplit=None):
    return SimpleNamespace(readlengthsplit=readlengthsplit)


def test_a_plain_build_writes_straight_to_results():
    """The common case: no split, so there is nothing to disambiguate and nothing to nest."""
    assert toolsTG.variant_dir_names(_args()) == ('results', 'graphs')


def test_a_split_build_nests_its_unsplit_variant_under_complete():
    assert toolsTG.variant_dir_names(_args(60)) == ('results/complete', 'graphs/complete')


def test_a_tagged_variant_always_nests():
    """A tag names one side of a split, so it nests whether or not the flag is on this call."""
    assert toolsTG.variant_dir_names(_args(60), 'u60') == ('results/u60', 'graphs/u60')
    assert toolsTG.variant_dir_names(_args(), 'o60') == ('results/o60', 'graphs/o60')


def test_results_and_graphs_stay_in_step():
    """adataBuild writes one and toolsSplit reads the other; a divergence here is a silent bug."""
    for args in (_args(), _args(60)):
        for tag in (None, 'u60'):
            results, graphs = toolsTG.variant_dir_names(args, tag)
            assert results.split('/')[1:] == graphs.split('/')[1:]


def test_a_workspace_left_over_from_the_two_build_era_is_still_ours(tmp_path):
    """`vibrChol1_nosplit` was created by an older version of the suite, so a workspace holding
    one must still be wipeable -- otherwise upgrading makes `--all` refuse to run until the user
    deletes a directory this suite put there."""
    for name in ('config', 'references', 'vibrChol1', 'vibrChol1_nosplit'):
        (tmp_path / name).mkdir()

    from trnagraph.modules import toolsTestSuite
    assert toolsTestSuite.unexpected_workspace_entries(str(tmp_path)) == []
