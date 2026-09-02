"""Read-specificity coverage, per group rather than aggregated over all of them.

`generate_partition_overview` was the only consumer of the four-way specificity partition, and
`_partition_frame` collapsed each category with `__coverage_transform__(..., singlecol=True)`
-- a --covmethod reduction over EVERY --covgrp value at once. So the plot ignored --covgrp
entirely and no per-group specificity view existed at any granularity: a treated and an
untreated sample were averaged into one trace with nothing saying so.

The replacement is a grid, tRNA down and group across, plus one individual plot per tRNA per
group.
"""
import anndata as ad
import numpy as np
import pandas as pd
import pytest

from trnagraph.modules import plotsCoverage, toolsTG


TRNAS = ['tRNA-Ala-AGC-1', 'tRNA-Gly-GCC-2']
GROUPS = ['ctrl', 'drug']
POSITIONS = 8
#: Per-group signal levels, chosen so a group-aware frame and an all-group average are
#: numerically distinguishable -- averaging these two gives 15, which is neither.
GROUP_LEVEL = {'ctrl': 10.0, 'drug': 20.0}


def _adata(groups=GROUPS, replicates=2):
    rows, values = [], []
    for trna in TRNAS:
        for group in groups:
            for replicate in range(replicates):
                rows.append({'trna': trna, 'group': group,
                             'sample': f'{group}_{replicate}'})
                values.append(np.full(POSITIONS * len(toolsTG.COVERAGE_PARTITION),
                                      GROUP_LEVEL.get(group, 10.0) / len(toolsTG.COVERAGE_PARTITION)))
    obs = pd.DataFrame(rows)
    obs.index = [f'obs{i}' for i in range(len(obs))]
    var = pd.DataFrame({
        'coverage': [c for c in toolsTG.COVERAGE_PARTITION for _ in range(POSITIONS)],
        'gap': [False] * (POSITIONS * len(toolsTG.COVERAGE_PARTITION)),
    })
    var.index = [f'v{i}' for i in range(len(var))]
    return ad.AnnData(X=np.array(values), obs=obs, var=var)


def _visualizer(**overrides):
    import matplotlib
    matplotlib.use('Agg')

    args = dict(adata=_adata(), threads=1, coverage_grp='group', coverage_obs='trna',
                coverage_type='uniquecoverage', coverage_gap=[False], coverage_method='mean',
                colormap=None, output='out/', settings=None)
    args.update(overrides)
    return plotsCoverage.visualizer(**args)


def test_the_partition_frame_is_keyed_by_covobs_and_by_group():
    frames = _visualizer()._partition_frame()

    assert set(frames) == set(TRNAS)
    assert set(frames[TRNAS[0]]) == set(GROUPS)


def test_each_group_carries_its_own_signal_rather_than_the_average_of_all_groups():
    """The defect in one assertion: with ctrl at 10 and drug at 20, the old frame reported
    15 for both -- a value neither group has."""
    frames = _visualizer()._partition_frame()

    for trna in TRNAS:
        for group in GROUPS:
            total = frames[trna][group].sum(axis=1)
            assert np.allclose(total, GROUP_LEVEL[group]), (
                f'{trna}/{group} should carry {GROUP_LEVEL[group]}, got {total.unique()}')


def test_every_partition_category_survives_into_each_groups_frame():
    """The four categories are what gets stacked, so losing one silently would understate
    coverage rather than fail."""
    frames = _visualizer()._partition_frame()

    for group_frames in frames.values():
        for frame in group_frames.values():
            assert list(frame.columns) == list(toolsTG.COVERAGE_PARTITION)


def _pages(visualizer):
    # _partition_pages streams its pages so the caller can save and release each one; these
    # tests inspect them all, so they materialize the stream.
    frames = visualizer._partition_frame()
    return list(visualizer._partition_pages(frames))


def test_the_grid_puts_a_trna_on_each_row_and_a_group_in_each_column():
    """Rows are tRNAs and columns are groups, so the same tRNA's groups sit side by side and
    can be read against each other -- which is the comparison the aggregated page destroyed."""
    pages = _pages(_visualizer())

    assert len(pages) == 1
    panels = pages[0].axes
    assert len(panels) == len(TRNAS) * len(GROUPS)


def test_every_group_is_named_on_the_page():
    """A grid whose columns are unlabelled is worse than the aggregate: it implies a
    distinction without saying what it is."""
    page = _pages(_visualizer())[0]

    rendered = ' '.join(text.get_text() for text in page.findobj(match=lambda o: hasattr(o, 'get_text')))
    for group in GROUPS:
        assert group in rendered
    for trna in TRNAS:
        assert trna in rendered


def test_more_groups_than_can_be_read_aborts_rather_than_capping_silently():
    """--covgrp trna on a human build is 436 columns. Plotting the first 24 would be a figure
    that misrepresents itself as complete, so the cardinality is refused outright."""
    # 30 groups, comfortably past the 24 a reader can take in across a page.
    wide = _adata(groups=[f'g{i}' for i in range(30)], replicates=1)
    visualizer = _visualizer(adata=wide, coverage_grp='group')

    with pytest.raises(toolsTG.InvalidParameterError) as raised:
        _pages(visualizer)
    message = str(raised.value)
    assert 'group' in message
    assert '30' in message


def test_panels_in_one_row_share_a_y_axis():
    """Two groups drawn to different scales look equal when they are not -- the specific way a
    small-multiples grid lies. The ctrl and drug fixtures differ 10 to 20, so an unshared axis
    would render them identically."""
    page = _pages(_visualizer())[0]
    panels = page.axes[:len(GROUPS)]

    tops = {ax.get_ylim()[1] for ax in panels}
    assert len(tops) == 1, f'panels in a row must share a ceiling, got {tops}'


def test_a_rows_scale_covers_every_group_not_only_those_on_this_page():
    """With more groups than fit on one page, a tRNA's second column-page must be drawn to the
    same scale as its first, or the two pages cannot be compared."""
    visualizer = _visualizer(adata=_adata(groups=[f'g{i}' for i in range(10)], replicates=1))
    pages = _pages(visualizer)

    tops = {ax.get_ylim()[1] for page in pages for ax in page.axes}
    assert len(tops) == 1


def test_the_grid_paginates_in_both_axes():
    """Six tRNAs and ten groups is two row-pages by two column-pages."""
    visualizer = _visualizer(
        adata=_adata(groups=[f'g{i}' for i in range(10)], replicates=1))
    visualizer_trnas = _pages(visualizer)

    # 2 tRNAs (one row-page) x 10 groups (two column-pages, 8 then 2)
    assert len(visualizer_trnas) == 2
    assert [len(page.axes) for page in visualizer_trnas] == [
        len(TRNAS) * 8, len(TRNAS) * 2]


def test_individual_plots_are_written_one_per_trna_per_group(tmp_path):
    """The per-group view at full size, mirroring how the coverage plots emit one file per
    group with `_with_endstarts_`."""
    visualizer = _visualizer(output=f'{tmp_path}/')
    visualizer.build_output_dirs()
    visualizer.generate_partition_split()

    written = sorted(p.name for p in tmp_path.rglob('*.pdf'))
    assert len(written) == len(TRNAS) * len(GROUPS)
    for trna in TRNAS:
        for group in GROUPS:
            assert any(trna in name and group in name for name in written), (
                f'no file for {trna}/{group} in {written}')


def test_a_low_coverage_group_is_diverted_like_the_coverage_plots(tmp_path):
    """Coverage sends a plot whose ceiling is under 20 into low_coverage/; specificity follows
    the same rule so the two sort identically. The fixture's ctrl group totals 10 and its drug
    group totals 20, so exactly one of each lands on each side."""
    visualizer = _visualizer(output=f'{tmp_path}/')
    visualizer.build_output_dirs()
    visualizer.generate_partition_split()

    low = {p.name for p in tmp_path.rglob('low_coverage/*.pdf')}
    ordinary = {p.name for p in tmp_path.rglob('*.pdf')} - low
    assert low and ordinary
    assert all('ctrl' in name for name in low)
    assert all('drug' in name for name in ordinary)


def test_specificity_output_lives_under_its_own_directory(tmp_path):
    """It shows every coverage category at once, so it belongs to no --covtype directory; and
    the combined page now sits beside the individuals it summarises."""
    visualizer = _visualizer(output=f'{tmp_path}/')
    visualizer.build_output_dirs()
    visualizer.generate_partition_split()
    visualizer.generate_partition_overview()

    written = [str(p.relative_to(tmp_path)) for p in tmp_path.rglob('*.pdf')]
    assert written, 'nothing was written'
    assert all(path.startswith('specificity/') for path in written), written
    assert any(path.startswith('specificity/combined_') for path in written), written


def _dispatch_calls(combinedpdfonly):
    """Which specificity methods adataGraph asks for, without doing any plotting."""
    import ast
    import inspect

    from trnagraph.modules import adataGraph

    source = inspect.getsource(adataGraph.anndataGrapher.dispatch_plot)
    tree = ast.parse(source.lstrip())
    return {node.func.attr for node in ast.walk(tree)
            if isinstance(node, ast.Call) and getattr(node.func, 'attr', None)}


def test_the_grapher_asks_for_both_specificity_views():
    calls = _dispatch_calls(combinedpdfonly=False)

    assert 'generate_partition_overview' in calls
    assert 'generate_partition_split' in calls, (
        'the individual specificity plots are never produced')


def test_the_individual_specificity_plots_sit_behind_combinedpdfonly():
    """--combinedpdfonly exists to skip the expensive per-tRNA files; specificity multiplies
    that count by the number of groups, so it must honour the same flag."""
    import ast
    import inspect

    from trnagraph.modules import adataGraph

    tree = ast.parse(inspect.getsource(adataGraph.anndataGrapher.dispatch_plot).lstrip())
    guarded = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.If):
            continue
        names = {n.attr for n in ast.walk(node.test) if isinstance(n, ast.Attribute)}
        if 'combinedpdfonly' not in names:
            continue
        guarded += [c.func.attr for c in ast.walk(node)
                    if isinstance(c, ast.Call) and getattr(c.func, 'attr', None)]

    assert 'generate_partition_split' in guarded


def test_an_individual_specificity_plot_carries_its_own_legend(tmp_path):
    """Four stacked colours with no key is a picture nobody can read. The grid states them
    once per page; a standalone file has no page to inherit from."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    drawn = []
    real_save = plotsCoverage.toolsTG.save_current

    def spy(path, settings=None):
        axes = plt.gcf().axes
        drawn.append([text.get_text() for text in axes[0].get_legend().get_texts()]
                     if axes and axes[0].get_legend() else [])
        return real_save(path, settings)

    visualizer = _visualizer(output=f'{tmp_path}/')
    visualizer.build_output_dirs()
    plotsCoverage.toolsTG.save_current = spy
    try:
        visualizer.generate_partition_split()
    finally:
        plotsCoverage.toolsTG.save_current = real_save

    assert drawn, 'nothing was saved'
    for entries in drawn:
        assert set(entries) == set(toolsTG.COVERAGE_CATEGORY_LABELS.values())


def _grapher(covgrp, graphtypes, adata):
    import logging
    from types import SimpleNamespace

    from trnagraph.modules import adataGraph

    grapher = object.__new__(adataGraph.anndataGrapher)
    grapher.adata = adata
    grapher.logger = logging.getLogger('test')
    grapher.args = SimpleNamespace(
        covgrp=covgrp, covobs='trna', clustergrp=None, clusterlabels=None, comparegrp1=None,
        comparegrp2=None, corrgroup=None, heatgrp=None, logogrp=None, pcamarkers=None,
        pcacolors=None, radargrp=None, volgrp=None, covtype=None, diffrts=None,
        pcareadtypes=None, graphtypes=graphtypes, corrmask='none', heatorient='vertical')
    return grapher


def test_an_unreadably_wide_covgrp_is_refused_before_any_plotting():
    """The guard inside the grid runs after generate_split and generate_combine, so on a human
    build the user would receive 400-odd coverage plots and only then be told the grouping
    cannot be drawn. Checked up front instead, with the other usage problems."""
    wide = _adata(groups=[f'g{i}' for i in range(30)], replicates=1)

    with pytest.raises(toolsTG.InvalidParameterError) as raised:
        _grapher('group', ['coverage'], wide)._validate_label_args()
    assert '30' in str(raised.value)


def test_the_width_limit_only_applies_when_coverage_is_being_drawn():
    """--covgrp is coverage's option alone, so a wide column is nobody's problem on a run that
    never draws a specificity grid."""
    wide = _adata(groups=[f'g{i}' for i in range(30)], replicates=1)

    _grapher('group', ['pca'], wide)._validate_label_args()
