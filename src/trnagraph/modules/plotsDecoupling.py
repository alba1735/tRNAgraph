#!/usr/bin/env python3
'''
The decoupling scatter: does a feature move the same way in two ways of measuring it?

A volcano answers "did this tRNA change between these conditions" for one measurement channel.
It cannot answer "did the fragment pool change the same way the full-length pool did", because
each figure sees one channel. This plot fits the SAME between-sample contrast separately in each
channel and plots one against the other: the diagonal is perfect coupling, and distance from it
is decoupling.

Two fits rather than one model. `~channel * timepoint` -- a real interaction term, tested by
likelihood ratio -- is the honest version of this question and belongs to the planned
multi-factor engine (roadmap.md's "trajectory / decoupling analysis"), which does not exist yet.
Fitting each channel separately answers it descriptively with machinery that already exists, and
nothing it reports has to be un-taught when the interaction test arrives.

It also sidesteps a normalization trap that a within-sample test walks straight into. The read
types partition one pool exactly -- wholecounts + fiveprime + threeprime + other == total, to
the read -- so 5' and 3' are not independent measurements, and a DESeq2 contrast between them
would estimate a size factor per sample x channel, absorbing the very difference being measured.
Comparing CHANGES within each channel never puts two channel levels side by side, so each
channel's size factors cancel inside its own fold change and the question never arises.
'''

import logging
import os

import numpy as np
import pandas as pd

from . import plotsPalette, toolsTG
from .toolsSchemas import DecouplingChannel, DecouplingPlan

import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


class InvalidDecouplingError(ValueError):
    '''A decoupling figure that cannot be built as specified.'''


#: The read type a channel takes when it names only a variant. `total` is the whole tRNA pool
#: for that read-length split, which is what "fragment versus full length" means.
DEFAULT_READTYPE = 'total'

#: Feature axes the figure is drawn over. `trna` carries the scatter: the decoupling question is
#: per-gene -- WHICH tRNAs fragment differently -- and aggregating to amino acid averages that
#: away. The coarser axes keep the readable bar summary the old `-g compare` plot provided.
SCATTER_AXIS = 'trna'
BAR_AXES = ('amino', 'iso')

#: Below this many shared features a scatter says nothing a table would not say better.
MIN_FEATURES_TO_DRAW = 3


def resolve_channel(adata, declaration, read_basis, index):
    '''
    One declared channel resolved to (tag, obs column), or None if this object cannot supply it.

    Returns None rather than raising for a missing split tag: whether a read-length split is
    meaningful is a property of the sequencing protocol, not of this pipeline. OTTR-seq carries
    both fragments and full-length reads so a cutoff informs; ARM-seq is fragment-biased and
    DM-tRNA-seq full-length-biased, and objects from either legitimately carry no split at all.
    A 5'/3' pair is a read-type axis every object has, so it must keep working there.
    '''
    tag = 'full'
    if declaration.variant:
        try:
            tag = toolsTG.parse_variant(adata, declaration.variant).tag
        except ValueError as unavailable:
            logger.warning(
                f'Skipping the {declaration.variant!r} decoupling channel: {unavailable}')
            return None
    bare = declaration.readtype or DEFAULT_READTYPE
    try:
        readtype = toolsTG.resolve_readtype(bare, read_basis, adata)
    except Exception as unavailable:
        logger.warning(f'Skipping the {bare!r} decoupling channel: {unavailable}')
        return None
    label = declaration.label or _default_label(declaration, tag, bare, index)
    return DecouplingChannel(label=label, readtype=readtype, tag=tag)


def _default_label(declaration, tag, bare, index):
    '''A label naming whatever actually distinguishes this channel, not "channel 2".'''
    if declaration.variant and tag != 'full':
        return f'Fragment ({tag})' if tag.startswith('u') else f'Full length ({tag})'
    if declaration.readtype:
        return plotsPalette.bare_readtype(bare).replace('_', ' ')
    return f'Channel {index + 1}'


def declared_plans(adata, block, read_basis=None):
    '''
    Resolve `multivariate.decoupling` into drawable plans, plus the reasons any were dropped.

    A pair is the unit that survives or is skipped: unlike a Venn -- whose circles are one joint
    claim, so drawing two of four declared would produce a figure whose title says four things
    and whose content shows two -- decoupling pairs are independent comparisons. Losing
    fragment/full-length on an object with no split therefore costs the 5'/3' pair nothing.
    '''
    basis = read_basis or toolsTG.READ_BASIS_UNIQUE
    plans, skipped = [], []
    for declaration in (block.decoupling or []):
        channels = [resolve_channel(adata, spec, basis, i)
                    for i, spec in enumerate(declaration.channels)]
        if any(channel is None for channel in channels):
            skipped.append(
                f"Skipping the {declaration.name!r} decoupling figure: this object cannot supply "
                f"both of its channels (see the warning above). Other declared pairs are "
                f"unaffected.")
            continue
        if channels[0].tag == channels[1].tag and channels[0].readtype == channels[1].readtype:
            skipped.append(
                f"Skipping the {declaration.name!r} decoupling figure: both channels resolve to "
                f"the same variant and read type, so there is nothing to compare.")
            continue
        plans.append(DecouplingPlan(
            name=declaration.name,
            title=declaration.title or declaration.name.replace('_', ' '),
            channels=channels))
    return plans, skipped


def channel_frames(adata, plan, grouping, cutoff, config_name, index_col,
                   overwrite=False, shrink='apeGLM'):
    '''
    Each channel's fold-change frame for `index_col`, keyed by channel label.

    The `trna` axis goes through toolsTG.variant_log2fc, the one path designed to span variants:
    it caches each tag's result where that tag's data lives (uns['log2FC'] for 'full',
    uns['size_splits'][tag]['log2FC'] otherwise), so no channel can serve or overwrite another's
    numbers.

    The coarser axes bypass that cache deliberately. Its key is
    [config][compare][readtype][cutoff] with no slot for a feature axis, so caching an amino-level
    frame under it would hand a later `trna` reader the wrong rows -- the same reason the plot
    this replaces called log2fc_df() directly. Those fits are cheap at 22-54 features.
    '''
    frames = {}
    for channel in plan.channels:
        if index_col == SCATTER_AXIS:
            frame, _ = toolsTG.variant_log2fc(
                adata, channel.tag, grouping, channel.readtype, readcount_cutoff=cutoff,
                config_name=config_name, overwrite=overwrite, shrink=shrink)
        else:
            frame, _ = _uncached_frame(adata, channel, grouping, cutoff, shrink, index_col)
        frames[channel.label] = frame
    return frames


def _uncached_frame(adata, channel, grouping, cutoff, shrink, index_col):
    '''One channel's fold changes over a coarse feature axis, fitted against a slim carrier.'''
    import anndata as ad

    obs = toolsTG.variant_obs(adata, channel.tag)
    carrier = ad.AnnData(X=np.zeros((len(obs), 1), dtype='float32'), obs=obs)
    engine = toolsTG.adataLog2FC(carrier, grouping, channel.readtype, readcount_cutoff=cutoff,
                                 shrink=shrink,
                                 dispfittype=toolsTG.recorded_dispfittype(adata))
    return engine.log2fc_df(index_col=index_col)


def decoupling_table(frames, plan, contrast, log2fc_threshold, padj):
    '''
    One row per feature shared by both channels: its fold change and call in each.

    Only features present in BOTH is the whole point -- a point needs an x and a y. `called_x`
    and `called_y` use the project-wide PVAL_SIG/log2FC pair, and `agreement` records whether
    the two channels say the same thing, which is what the colour encodes.
    '''
    x_label, y_label = (channel.label for channel in plan.channels)
    log2_col, pval_col = f'log2_{contrast}', f'pval_{contrast}'
    x_frame, y_frame = frames[x_label], frames[y_label]
    for frame in (x_frame, y_frame):
        if log2_col not in frame.columns:
            raise InvalidDecouplingError(
                f'Contrast {contrast!r} is missing from a channel frame; the two channels were '
                f'fitted over different groupings.')
    shared = x_frame.index.intersection(y_frame.index)
    table = pd.DataFrame({
        'feature': list(shared),
        'x': x_frame.loc[shared, log2_col].to_numpy(),
        'y': y_frame.loc[shared, log2_col].to_numpy(),
        'x_padj': x_frame.loc[shared, pval_col].to_numpy(),
        'y_padj': y_frame.loc[shared, pval_col].to_numpy(),
    })
    table = table.dropna(subset=['x', 'y'])
    called_x = (table['x'].abs() >= log2fc_threshold) & (table['x_padj'] <= padj)
    called_y = (table['y'].abs() >= log2fc_threshold) & (table['y_padj'] <= padj)
    table['called_x'], table['called_y'] = called_x.fillna(False), called_y.fillna(False)
    return table


#: What each colour means. Ordered weakest claim first so the legend reads as a progression.
TIER_NEITHER = 'Neither'
TIER_X_ONLY = 'x only'
TIER_Y_ONLY = 'y only'
TIER_BOTH = 'Both'


def channel_palette(plan, colormap=None):
    '''
    One colour per channel, honouring a `colors.decoupling` block keyed by channel label.

    Keyed by LABEL rather than derived from a grouping level, for the reason the Venn's block is:
    a channel is not a level of any obs column, so there is no entry elsewhere in the style file
    that names it. Unnamed channels take the positional colourblind-safe pair, whose first two
    entries are the blue/orange that stay distinct under every common deficiency.
    '''
    defaults = plotsPalette.venn_colors(len(plan.channels))
    return [(colormap or {}).get(channel.label) or default
            for channel, default in zip(plan.channels, defaults)]


def _tier_colors(plan, colormap=None):
    '''
    Colour per tier. The two single-channel tiers take their channel's own colour, so a point
    "significant in the fragment channel only" is drawn in whatever the rest of the figure set
    already uses for fragments -- the reader learns one encoding, not two.
    '''
    palette = channel_palette(plan, colormap)
    return {
        TIER_NEITHER: plotsPalette.DE_NONSIGNIFICANT,
        TIER_X_ONLY: palette[0],
        TIER_Y_ONLY: palette[1],
        TIER_BOTH: plotsPalette.DE_UP,
    }


def _tier_of(row):
    if row['called_x'] and row['called_y']:
        return TIER_BOTH
    if row['called_x']:
        return TIER_X_ONLY
    if row['called_y']:
        return TIER_Y_ONLY
    return TIER_NEITHER


def label_rows(table, toplabels):
    '''
    Which points get a text label: the most decoupled of the called features.

    Ranked by distance from the diagonal, |x - y|, because that IS the quantity this figure
    exists to show -- the volcano ranks by |log2FC| x -log10(p) since significance is its axis,
    and the same reasoning lands somewhere different here.

    Restricted to features called in at least ONE channel first. Without that, the top of the
    ranking fills with noise: a feature at x = 0.1, y = -0.9 is a full point from the diagonal
    and means nothing, while the labels a reader wants are on real movement that one channel saw
    and the other did not.

    `toplabels` is a CEILING rather than a promise -- _place_labels culls whatever cannot be
    drawn legibly -- and 0 disables labelling entirely, matching --vollabels elsewhere.
    '''
    called = table[table['called_x'] | table['called_y']]
    if called.empty or toplabels == 0:
        return called.iloc[0:0]
    ranked = called.assign(_decoupling=(called['x'] - called['y']).abs())
    ranked = ranked.sort_values('_decoupling', ascending=False)
    return ranked if toplabels is None else ranked.head(toplabels)


def draw_scatter(ax, table, plan, contrast, tier_colors, settings=None, show_legend=True,
                 toplabels=100):
    '''
    One channel per axis, one point per feature, with the y = x line of perfect coupling.

    The diagonal is drawn from the data's own square extent rather than the axes limits, so it
    stays a true 45 degrees whatever the spread -- a "diagonal" that is not at 45 degrees makes
    distance from it unreadable, which is the only thing this figure is for.
    '''
    settings = settings or {}
    x_label, y_label = (channel.label for channel in plan.channels)
    marker_size = settings.get('marker_size') or 18
    threshold = settings.get('rasterize_over')
    rasterized = bool(threshold) and len(table) > threshold

    tiers = table.apply(_tier_of, axis=1) if len(table) else pd.Series(dtype=object)
    limit = float(np.nanmax(np.abs(np.concatenate([table['x'], table['y']])))) * 1.08 if len(table) else 1.0
    limit = max(limit, 0.5)

    ax.axhline(0, color=plotsPalette.GRID_LINE, linewidth=0.6, zorder=0)
    ax.axvline(0, color=plotsPalette.GRID_LINE, linewidth=0.6, zorder=0)
    ax.plot([-limit, limit], [-limit, limit], color=plotsPalette.REFERENCE_LINE,
            linestyle='--', linewidth=0.9, zorder=1)

    # Drawn weakest tier first so a strong point is never hidden under a grey one.
    for tier in (TIER_NEITHER, TIER_X_ONLY, TIER_Y_ONLY, TIER_BOTH):
        mask = (tiers == tier) if len(table) else []
        if len(table) == 0 or not np.any(mask):
            continue
        ax.scatter(table.loc[mask, 'x'], table.loc[mask, 'y'],
                   c=[tier_colors[tier]] * int(np.sum(mask)), s=marker_size, linewidths=0,
                   alpha=0.85 if tier != TIER_NEITHER else 0.35, zorder=2,
                   rasterized=rasterized, label=_tier_label(tier, x_label, y_label))

    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect('equal', adjustable='box')

    # After the limits are final: adjust_text solves in display space, so labels placed against
    # provisional limits would be repositioned by the rescale and land back on their markers.
    labelled = label_rows(table, toplabels)
    if not labelled.empty:
        from . import plotsVolcano

        # The FULL point set is handed over for repulsion, with only the ranked subset named --
        # so a label is pushed off every marker on the figure, not merely off the labelled ones.
        positions_x = pd.Series(table['x'].to_numpy(), index=pd.Index(table['feature']))
        positions_y = pd.Series(table['y'].to_numpy(), index=pd.Index(table['feature']))
        plotsVolcano._place_labels(ax, pd.Index(labelled['feature']), positions_x, positions_y)
    ax.set_xlabel(f'log2FC — {x_label}')
    ax.set_ylabel(f'log2FC — {y_label}')
    ax.set_title(f'{plan.title}\n{contrast}', fontsize=10)
    if show_legend:
        ax.legend(loc='upper left', frameon=False, fontsize=7, markerscale=1.4)
    return ax


def _tier_label(tier, x_label, y_label):
    if tier == TIER_X_ONLY:
        return f'{x_label} only'
    if tier == TIER_Y_ONLY:
        return f'{y_label} only'
    return 'Significant in both' if tier == TIER_BOTH else 'Neither'


def draw_bars(ax, table, plan, contrast, tier_colors, settings=None):
    '''
    The coarse axes as paired bars: one row per amino acid or isoacceptor, one bar per channel.

    Kept from the plot this replaces. A scatter carries hundreds of points comfortably and the
    bar chart cannot, but at 22 amino acids the bars stay the more readable summary of the two.
    '''
    settings = settings or {}
    x_label, y_label = (channel.label for channel in plan.channels)
    ordered = table.reindex(table[['x', 'y']].abs().max(axis=1).sort_values().index)
    positions = np.arange(len(ordered))
    height = 0.38
    edgewidth = toolsTG.linewidth_for(settings, 0.0)

    ax.barh(positions + height / 2, ordered['x'], height=height, color=tier_colors[TIER_X_ONLY],
            linewidth=edgewidth, label=x_label)
    ax.barh(positions - height / 2, ordered['y'], height=height, color=tier_colors[TIER_Y_ONLY],
            linewidth=edgewidth, label=y_label)
    ax.axvline(0, color=plotsPalette.GRID_LINE, linewidth=0.8)
    ax.set_yticks(positions)
    ax.set_yticklabels(ordered['feature'])
    ax.set_xlabel('log2 Fold-Change')
    ax.set_title(f'{plan.title}\n{contrast}', fontsize=10)
    ax.legend(loc='lower right', frameon=False, fontsize=8)
    return ax


def _write_table(table, plan, contrast, index_col, results_dir, name):
    '''The figure's own numbers as a TSV, rendered from the same frame the scatter draws.'''
    if not results_dir:
        return None
    path = os.path.join(results_dir, f'{name}_{index_col}_decoupling.tsv')
    out = table.rename(columns={
        'feature': index_col,
        'x': f'log2_{plan.channels[0].label}', 'y': f'log2_{plan.channels[1].label}',
        'x_padj': f'padj_{plan.channels[0].label}', 'y_padj': f'padj_{plan.channels[1].label}',
        'called_x': f'called_{plan.channels[0].label}',
        'called_y': f'called_{plan.channels[1].label}',
    })
    out.insert(1, 'contrast', contrast)
    out.to_csv(path, sep='\t', index=False)
    return path


def visualizer(adata, block, output, config_name='default', settings=None, read_basis=None,
               cutoff=20, colormap=None, threaded=False, overwrite=False, shrink='apeGLM',
               toplabels=100):
    '''
    Draw every declared decoupling pair, over every reference-anchored contrast.

    A combined scatter overview per pair, the per-contrast scatters and the coarse-axis bar
    charts under individual/, and a TSV per drawn figure. Contrasts come from
    `multivariate.grouping`/`reference` -- the same universe `-g agreement` uses -- so every
    figure in a run shares one denominator rather than mixing pairwise combinations.
    '''
    from . import plotsAgreement, plotsVenn

    grouping = block.grouping
    if grouping not in adata.obs.columns:
        raise InvalidDecouplingError(
            f'`multivariate.grouping` names {grouping!r}, which is not a column of obs.')
    reference = plotsAgreement.resolve_reference(adata.obs[grouping], block)
    levels = (list(adata.obs[grouping].cat.categories)
              if hasattr(adata.obs[grouping], 'cat')
              else list(dict.fromkeys(adata.obs[grouping])))
    contrasts = plotsAgreement.contrast_universe(levels, reference)

    plans, skipped = declared_plans(adata, block, read_basis=read_basis)
    messages = list(skipped)
    if not plans:
        messages.append('No decoupling figures could be drawn from this object.')
        return _emit(messages, threaded)

    individual_output = f'{output}individual/'
    messages.append(toolsTG.builder(individual_output))
    results_dir, results_message = plotsVenn.resolve_results_dir(adata)
    if results_message:
        messages.append(results_message)
    elif results_dir:
        os.makedirs(results_dir, exist_ok=True)

    for plan in plans:
        tier_colors = _tier_colors(plan, colormap)
        scatter_frames = channel_frames(adata, plan, grouping, cutoff, config_name,
                                        SCATTER_AXIS, overwrite=overwrite, shrink=shrink)
        drawn = []
        for reference_level, level in contrasts:
            contrast = f'{reference_level}-{level}'
            table = decoupling_table(scatter_frames, plan, contrast, block.log2fc, block.padj)
            if len(table) < MIN_FEATURES_TO_DRAW:
                messages.append(
                    f'Skipping {plan.name} {contrast}: only {len(table)} feature(s) are '
                    f'present in both channels, which is too few to read as a scatter.')
                continue
            drawn.append((contrast, table))

            fig, ax = plt.subplots(figsize=toolsTG.figsize_for(settings, (4.2, 4.2)))
            draw_scatter(ax, table, plan, contrast, tier_colors, settings=settings,
                         toplabels=toplabels)
            path = f'{individual_output}{plan.name}_{contrast}_decoupling.pdf'
            toolsTG.save_current(path, settings)
            plt.close(fig)
            messages.append(f'Saving decoupling scatter: {path}')

            written = _write_table(table, plan, contrast, SCATTER_AXIS, results_dir,
                                   f'{plan.name}_{contrast}')
            if written:
                messages.append(f'Saving decoupling table: {written}')

        # Combined overview: every contrast for this pair on one page. Deliberately NOT wired to
        # --style figsize -- a multi-panel page computes its geometry from how many panels it is
        # laying out, so a fixed size would either clip panels or leave dead space.
        if drawn:
            width = min(len(drawn), 3)
            rows = int(np.ceil(len(drawn) / width))
            fig, axes = plt.subplots(rows, width, figsize=(4.4 * width, 4.4 * rows),
                                     squeeze=False)
            for axis, (contrast, table) in zip(axes.ravel(), drawn):
                draw_scatter(axis, table, plan, contrast, tier_colors, toplabels=toplabels,
                             show_legend=(contrast == drawn[0][0]))
            for axis in axes.ravel()[len(drawn):]:
                axis.axis('off')
            fig.tight_layout()
            combined = f'{output}{plan.name}_decoupling.pdf'
            fig.savefig(combined, bbox_inches='tight')
            plt.close(fig)
            messages.append(f'Saving combined decoupling overview: {combined}')

        # The coarse axes, as the bar summary the plot this replaces provided.
        for axis_name in BAR_AXES:
            if axis_name not in adata.obs.columns:
                continue
            bar_frames = channel_frames(adata, plan, grouping, cutoff, config_name, axis_name,
                                        overwrite=overwrite, shrink=shrink)
            for reference_level, level in contrasts:
                contrast = f'{reference_level}-{level}'
                try:
                    table = decoupling_table(bar_frames, plan, contrast, block.log2fc, block.padj)
                except InvalidDecouplingError:
                    continue
                if table.empty:
                    continue
                fig, ax = plt.subplots(figsize=toolsTG.figsize_for(
                    settings, (5.5, max(3.0, 0.26 * len(table)))))
                draw_bars(ax, table, plan, contrast, tier_colors, settings=settings)
                path = f'{individual_output}{plan.name}_{contrast}_{axis_name}_decoupling.pdf'
                toolsTG.save_current(path, settings)
                plt.close(fig)
                messages.append(f'Saving decoupling bars: {path}')

    return _emit(messages, threaded)


def _emit(messages, threaded):
    '''
    Accumulate into the pooled-run string, or log directly when running serially.

    Each message is logged EXACTLY once, at a level taken from what it says. Logging a skip at
    its point of creation and then again here is how the sibling plots end up printing the same
    sentence twice, once as WARNING and once as INFO, which reads as two separate events.
    '''
    if threaded:
        return (threaded if isinstance(threaded, str) else '') + '\n'.join(messages) + '\n'
    for message in messages:
        logger.warning(message) if message.startswith('Skipping') else logger.info(message)
    return threaded
