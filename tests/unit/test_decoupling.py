"""`-g decoupling`: does a feature move the same way in two ways of measuring it?

The plot replaces `-g compare`, which stratified by an obs column and drew fold changes per
amino acid. That plot asked a question about two sequencing methods and could not be made to
mean anything about one experiment. This one fits the SAME between-sample contrast separately in
each of two measurement channels -- a read-length variant, or a read type -- and plots one
against the other, so the diagonal is perfect coupling and distance from it is decoupling.

Two fits rather than `~channel * timepoint`: the interaction test belongs to the multi-factor
engine that does not exist yet, and comparing CHANGES within each channel avoids the trap a
within-sample test walks into -- the read types partition one pool exactly, so a DESeq2 contrast
between them would estimate a size factor per sample x channel and absorb the very difference
being measured.
"""
import anndata as ad
import numpy as np
import pandas as pd
import pytest

from trnagraph.modules import plotsDecoupling, plotsPalette
from trnagraph.modules.toolsSchemas import (ChannelDeclaration, DecouplingChannel,
                                            DecouplingPlan, MultivariateConfig)

READTYPES = ['wholecounts', 'fiveprime', 'threeprime', 'other', 'total']


def _adata(with_split=True, n_trna=4):
    rows = []
    for t in range(n_trna):
        for sample in ('A_1', 'A_2', 'B_1', 'B_2'):
            row = {'trna': f'tRNA-Ala-AGC-{t}', 'amino': 'Ala', 'iso': 'Ala-AGC',
                   'sample': sample, 'group': sample.split('_')[0]}
            for rt in READTYPES:
                row[f'nreads_{rt}_unique_raw'] = 100
                row[f'nreads_{rt}_unique_norm'] = 100.0
                row[f'nreads_{rt}_raw'] = 120
                row[f'nreads_{rt}_norm'] = 120.0
            rows.append(row)
    obs = pd.DataFrame(rows)
    obs.index = [f'{r.trna}_{r["sample"]}' for _, r in obs.iterrows()]
    adata = ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs)
    if with_split:
        adata.uns['size_splits'] = {'u60': {}, 'o60': {}}
        for tag in ('u60', 'o60'):
            adata.obsm[f'size_split_{tag}'] = obs[[c for c in obs.columns
                                                   if c.startswith('nreads_')]].copy()
            adata.layers[f'norm_{tag}'] = np.zeros_like(adata.X)
    return adata


def _block(pairs):
    return MultivariateConfig.model_validate({'grouping': 'group', 'decoupling': pairs})


FRAG_FULL = {'name': 'frag_vs_full',
             'channels': [{'label': 'Fragment', 'variant': 'norm:u60'},
                          {'label': 'Full length', 'variant': 'norm:o60'}]}
ENDS = {'name': 'ends',
        'channels': [{'label': "5'", 'readtype': 'fiveprime'},
                     {'label': "3'", 'readtype': 'threeprime'}]}


# --- schema ---------------------------------------------------------------------------

def test_a_comparison_needs_exactly_two_channels():
    """The figure plots one channel against the other, so a third has no axis to live on."""
    with pytest.raises(ValueError, match='exactly two channels'):
        MultivariateConfig.model_validate(
            {'decoupling': [{'name': 'x', 'channels': [{'variant': 'norm:u60'}]}]})


def test_a_channel_cannot_carry_a_grouping_level():
    """ChannelDeclaration is purpose-built rather than a VennSetDeclaration with one field
    forbidden: a level is what a decoupling CONTRAST varies, not what distinguishes channels."""
    with pytest.raises(ValueError):
        ChannelDeclaration.model_validate({'variant': 'norm:u60', 'level': 'Day 0'})


# --- resolution and graceful degradation ------------------------------------------------

def test_a_readtype_channel_resolves_on_any_object():
    """5'/3' is a read-type axis every object has, so it must work regardless of protocol."""
    channel = plotsDecoupling.resolve_channel(
        _adata(with_split=False), ChannelDeclaration(readtype='fiveprime'), 'unique', 0)

    assert channel is not None
    assert channel.tag == 'full'
    assert 'fiveprime' in channel.readtype


def test_a_missing_split_tag_yields_no_channel_rather_than_raising():
    """Whether a read-length split is meaningful is a property of the sequencing protocol.
    OTTR-seq carries both fragments and full-length reads so a cutoff informs; ARM-seq is
    fragment-biased and DM-tRNA-seq full-length-biased, and objects from either legitimately
    carry no split at all."""
    channel = plotsDecoupling.resolve_channel(
        _adata(with_split=False), ChannelDeclaration(variant='norm:u60'), 'unique', 0)

    assert channel is None


def test_an_unavailable_pair_is_skipped_without_costing_the_others():
    """The unit that survives or is skipped is the PAIR. Unlike a Venn -- whose circles are one
    joint claim, so drawing two of four declared would produce a figure whose title says four
    things and whose content shows two -- decoupling pairs are independent comparisons."""
    plans, skipped = plotsDecoupling.declared_plans(
        _adata(with_split=False), _block([FRAG_FULL, ENDS]), read_basis='unique')

    assert [p.name for p in plans] == ['ends'], 'the read-type pair must survive'
    assert len(skipped) == 1 and 'frag_vs_full' in skipped[0]


def test_both_pairs_survive_when_the_object_has_the_split():
    plans, skipped = plotsDecoupling.declared_plans(
        _adata(), _block([FRAG_FULL, ENDS]), read_basis='unique')

    assert [p.name for p in plans] == ['frag_vs_full', 'ends']
    assert skipped == []


def test_two_identical_channels_are_refused():
    """Naming one channel twice leaves nothing to compare, the same defect the deleted compare
    plot had when both of its grouping flags named one column."""
    same = {'name': 'same', 'channels': [{'readtype': 'fiveprime'}, {'readtype': 'fiveprime'}]}
    plans, skipped = plotsDecoupling.declared_plans(_adata(), _block([same]),
                                                    read_basis='unique')

    assert plans == []
    assert 'nothing to compare' in skipped[0]


# --- the table the figure draws from ------------------------------------------------------

def _plan():
    return DecouplingPlan(name='p', title='p', channels=[
        DecouplingChannel(label='X', readtype='nreads_total_unique_norm', tag='u60'),
        DecouplingChannel(label='Y', readtype='nreads_total_unique_norm', tag='o60')])


def _frames(x_rows, y_rows):
    def frame(rows):
        return pd.DataFrame({'log2_A-B': [r[0] for r in rows.values()],
                             'pval_A-B': [r[1] for r in rows.values()]},
                            index=list(rows))
    return {'X': frame(x_rows), 'Y': frame(y_rows)}


def test_only_features_present_in_both_channels_are_plotted():
    """A point needs an x and a y."""
    frames = _frames({'a': (1.0, 0.01), 'b': (2.0, 0.01)},
                     {'a': (1.0, 0.01), 'c': (3.0, 0.01)})

    table = plotsDecoupling.decoupling_table(frames, _plan(), 'A-B', 1.5, 0.001)

    assert list(table['feature']) == ['a']


def test_a_feature_called_in_one_channel_only_is_the_decoupled_case():
    """'Significant in the fragment channel but not the full-length one' IS decoupling, which is
    why it is what the colour encodes rather than a secondary annotation."""
    frames = _frames({'a': (3.0, 0.0001), 'b': (0.1, 0.9)},
                     {'a': (0.1, 0.9), 'b': (0.2, 0.9)})

    table = plotsDecoupling.decoupling_table(frames, _plan(), 'A-B', 1.5, 0.001)
    row = table.set_index('feature').loc['a']

    assert bool(row['called_x']) and not bool(row['called_y'])
    assert plotsDecoupling._tier_of(row) == plotsDecoupling.TIER_X_ONLY


def test_a_feature_moving_together_in_both_channels_is_coupled():
    frames = _frames({'a': (3.0, 0.0001)}, {'a': (2.9, 0.0001)})

    table = plotsDecoupling.decoupling_table(frames, _plan(), 'A-B', 1.5, 0.001)

    assert plotsDecoupling._tier_of(table.iloc[0]) == plotsDecoupling.TIER_BOTH


def test_a_feature_below_both_thresholds_is_neither():
    frames = _frames({'a': (0.2, 0.9)}, {'a': (0.1, 0.9)})

    table = plotsDecoupling.decoupling_table(frames, _plan(), 'A-B', 1.5, 0.001)

    assert plotsDecoupling._tier_of(table.iloc[0]) == plotsDecoupling.TIER_NEITHER


def test_a_missing_contrast_names_the_cause():
    frames = _frames({'a': (1.0, 0.01)}, {'a': (1.0, 0.01)})

    with pytest.raises(plotsDecoupling.InvalidDecouplingError, match='different groupings'):
        plotsDecoupling.decoupling_table(frames, _plan(), 'C-D', 1.5, 0.001)


# --- colour ---------------------------------------------------------------------------

def test_channels_take_the_colourblind_safe_pair_by_default():
    assert plotsDecoupling.channel_palette(_plan()) == plotsPalette.venn_colors(2)


def test_a_style_entry_names_a_channel_by_label():
    """Keyed by label rather than derived from a grouping level: a channel is not a level of any
    obs column, so there is no entry elsewhere in the style file that names it."""
    palette = plotsDecoupling.channel_palette(_plan(), {'Y': '#123456'})

    assert palette == [plotsPalette.venn_colors(2)[0], '#123456']


def test_the_single_channel_tiers_take_their_own_channels_colour():
    """A point 'significant in the fragment channel only' is drawn in whatever the rest of the
    figure set already uses for fragments -- the reader learns one encoding, not two."""
    tiers = plotsDecoupling._tier_colors(_plan(), {'X': '#111111', 'Y': '#222222'})

    assert tiers[plotsDecoupling.TIER_X_ONLY] == '#111111'
    assert tiers[plotsDecoupling.TIER_Y_ONLY] == '#222222'
    assert tiers[plotsDecoupling.TIER_NEITHER] == plotsPalette.DE_NONSIGNIFICANT


# --- gating -----------------------------------------------------------------------------

def _grapher(pairs, graphtypes):
    """A grapher carrying a real RunConfig -- the status check reads `config.multivariate`,
    so handing it a bare MultivariateConfig would test the fixture rather than the gate."""
    import logging
    from types import SimpleNamespace

    from trnagraph.modules import adataGraph
    from trnagraph.modules.toolsSchemas import RunConfig

    config = RunConfig.model_validate(
        {'name': 'test', 'multivariate': {'grouping': 'group', 'decoupling': pairs} if pairs
         else {'grouping': 'group'}})
    grapher = object.__new__(adataGraph.anndataGrapher)
    grapher.config = config
    grapher.args = SimpleNamespace(graphtypes=graphtypes)
    grapher.logger = logging.getLogger('test')
    return grapher


def test_decoupling_needs_its_own_block_not_merely_a_multivariate_one():
    """venn and agreement are satisfied by the block existing. decoupling needs the channel
    pairs enumerated: which two ways of measuring the same features are worth comparing is a
    claim about the experiment, and an object carrying several split tags gives no basis for
    picking one."""
    status = _grapher(None, ['all'])._optional_graphtype_status()

    assert status['venn'] is None, 'the block alone satisfies venn'
    assert status['decoupling'] is not None
    assert 'decoupling' in status['decoupling']


def test_decoupling_is_ready_once_a_pair_is_declared():
    status = _grapher([ENDS], ['all'])._optional_graphtype_status()

    assert status['decoupling'] is None


def test_a_declared_pair_folds_decoupling_into_all():
    """Writing the config block IS the deliberate opt-in, which is the argument the codebase
    already makes for venn -- so the gate is the prerequisite, not the name."""
    from trnagraph.modules import adataGraph

    types, skipped = adataGraph.resolve_graphtypes(
        ['all'], {'venn': None, 'agreement': None, 'decoupling': None})

    assert 'decoupling' in types
    assert skipped == []


def test_without_a_declared_pair_all_leaves_it_out_with_the_reason():
    from trnagraph.modules import adataGraph

    types, skipped = adataGraph.resolve_graphtypes(
        ['all'], {'venn': None, 'agreement': None, 'decoupling': 'no channels declared'})

    assert 'decoupling' not in types
    assert ('decoupling', 'no channels declared') in skipped


# --- labels -----------------------------------------------------------------------------

def _label_table(rows):
    """rows: {feature: (x, y, called_x, called_y)}."""
    return pd.DataFrame({
        'feature': list(rows),
        'x': [r[0] for r in rows.values()], 'y': [r[1] for r in rows.values()],
        'x_padj': 0.0, 'y_padj': 0.0,
        'called_x': [r[2] for r in rows.values()], 'called_y': [r[3] for r in rows.values()],
    })


def test_labels_rank_by_distance_from_the_diagonal():
    """The volcano ranks by |log2FC| x -log10(p) because significance is its axis. Here the axis
    of interest is decoupling itself, so |x - y| is what orders the labels."""
    table = _label_table({'near': (2.0, 1.9, True, True),
                          'far': (3.0, -1.0, True, False),
                          'mid': (2.0, 0.5, True, False)})

    assert list(label_order(table)) == ['far', 'mid', 'near']


def label_order(table, toplabels=None):
    return plotsDecoupling.label_rows(table, toplabels)['feature']


def test_an_uncalled_feature_is_never_labelled_however_far_from_the_diagonal():
    """Without the called-in-one-channel restriction the ranking fills with noise: a feature at
    x = 0.1, y = -0.9 is a full point off the diagonal and means nothing."""
    table = _label_table({'noise': (0.1, -0.9, False, False),
                          'real': (2.0, 1.8, True, True)})

    assert list(label_order(table)) == ['real']


def test_toplabels_caps_the_number_requested():
    table = _label_table({'a': (3.0, -1.0, True, False), 'b': (2.0, 0.5, True, False),
                          'c': (1.6, 1.5, True, True)})

    assert list(label_order(table, toplabels=2)) == ['a', 'b']


def test_zero_disables_labelling_entirely():
    """Matching --vollabels elsewhere, where 0 means "draw none"."""
    table = _label_table({'a': (3.0, -1.0, True, False)})

    assert plotsDecoupling.label_rows(table, 0).empty


def test_none_labels_every_called_feature():
    table = _label_table({'a': (3.0, -1.0, True, False), 'b': (0.1, 0.0, False, False)})

    assert list(label_order(table, toplabels=None)) == ['a']


def test_the_scatter_actually_draws_the_labels():
    """End to end through the renderer, since ranking correctly and drawing nothing would pass
    every test above."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    table = _label_table({'tRNA-Arg-TCT-4': (3.0, -1.0, True, False),
                          'tRNA-Ala-AGC-1': (2.0, 1.9, True, True),
                          'tRNA-Gly-GCC-2': (0.1, 0.0, False, False)})
    plan = _plan()
    fig, ax = plt.subplots()
    try:
        plotsDecoupling.draw_scatter(ax, table, plan, 'A-B',
                                     plotsDecoupling._tier_colors(plan), toplabels=100)
        drawn = {t.get_text() for t in ax.texts}
    finally:
        plt.close(fig)

    assert 'tRNA-Arg-TCT-4' in drawn
    assert 'tRNA-Gly-GCC-2' not in drawn, 'an uncalled feature must not be labelled'


def test_the_scatter_draws_no_labels_when_disabled():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    table = _label_table({'tRNA-Arg-TCT-4': (3.0, -1.0, True, False)})
    plan = _plan()
    fig, ax = plt.subplots()
    try:
        plotsDecoupling.draw_scatter(ax, table, plan, 'A-B',
                                     plotsDecoupling._tier_colors(plan), toplabels=0)
        drawn = [t.get_text() for t in ax.texts]
    finally:
        plt.close(fig)

    assert drawn == []


# --- the exported table -------------------------------------------------------------------

def test_the_results_directory_is_a_sister_of_the_figure_directory():
    """Tables under results/, figures under graphs/, so a reader finds the table by the same
    route they found the plot."""
    from trnagraph.modules import toolsTG

    assert toolsTG.results_mirror('/out/graphs') == '/out/results'


def test_the_sister_holds_whatever_the_output_directory_is_called():
    """A user who points --output at figures/ still gets a sister, not a results/ nested inside
    their figures -- and a user whose project lives in ~/graphs_2024/ does not have that
    rewritten out from under them, which plotsHeatmap's replace('graphs', 'results') does."""
    from trnagraph.modules import toolsTG

    assert toolsTG.results_mirror('/out/figures') == '/out/results'
    assert toolsTG.results_mirror('/home/u/graphs_2024/out/graphs') == \
        '/home/u/graphs_2024/out/results'


def _export(rows):
    return plotsDecoupling.ranked_table(_label_table(rows), _plan(), 'A-B', 'trna')


def test_the_table_is_ordered_the_way_the_labels_are():
    """One notion of "notable" across the figure and the file beside it: called features first,
    then by distance from the diagonal. A table that disagreed with its own plot about what
    mattered would be worse than no table."""
    out = _export({'coupled_big': (4.0, 3.9, True, True),
                   'decoupled': (3.0, 0.5, True, False),
                   'uncalled_far': (0.9, -0.9, False, False)})

    assert list(out['trna']) == ['decoupled', 'coupled_big', 'uncalled_far']
    assert list(out['rank']) == [1, 2, 3]


def test_a_feature_that_moved_a_lot_in_both_channels_is_not_decoupled():
    """Significance alone is the wrong sort for this figure: a strong hit sitting on the diagonal
    is the coupled case, which is the opposite of what it is looking for."""
    out = _export({'on_diagonal': (5.0, 5.0, True, True),
                   'off_diagonal': (2.0, 0.2, True, False)})

    assert out.iloc[0]['trna'] == 'off_diagonal'
    assert out.set_index('trna').loc['on_diagonal', 'decoupling'] == pytest.approx(0.0)


def test_the_table_carries_what_a_reader_needs_to_re_sort_it():
    out = _export({'a': (3.0, 0.5, True, False)})

    for column in ('rank', 'trna', 'contrast', 'decoupling', 'min_padj', 'tier'):
        assert column in out.columns
    assert 'log2_X' in out.columns and 'padj_Y' in out.columns
    assert out.iloc[0]['tier'] == 'X only'


def test_min_padj_is_the_stronger_of_the_two_channels():
    """Carried so a reader who does want to sort by raw significance can, without recomputing."""
    out = _export({'a': (3.0, 0.5, True, False)})

    assert out.iloc[0]['min_padj'] == pytest.approx(0.0)


def test_the_decoupling_column_is_distance_from_the_diagonal():
    out = _export({'a': (3.0, 0.5, True, False)})

    assert out.iloc[0]['decoupling'] == pytest.approx(2.5)
