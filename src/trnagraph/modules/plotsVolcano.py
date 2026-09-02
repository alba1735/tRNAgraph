#!/usr/bin/env python3

import logging

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from adjustText import adjust_text

from . import toolsTG
from . import plotsPalette
from . import plotsThresholds

logger = logging.getLogger(__name__)

plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Significance boundaries. LOG2FC_THRESHOLD + PVAL_SIG_LINE together define which points count
# as "significant" for coloring/opacity/labeling; PVAL_STRONG_LINE is drawn as an additional
# (non-classifying) reference line, matching the pre-existing dashed-line behavior.
#
# Derived from plotsThresholds rather than written as -log10 literals. These were previously
# 1.3 and 3, and 1.3 is a ROUNDED -log10(0.05) -- so the line labelled p = 0.05 sat at
# p = 0.050119, and every feature in that sliver was classified significant while failing the
# threshold the figure drew. Deriving the axis position from the p-value makes the two agree
# by construction.
LOG2FC_THRESHOLD = plotsThresholds.LOG2FC_THRESHOLD
PVAL_SIG_LINE = plotsThresholds.neglog10(plotsThresholds.PVAL_SIG)
PVAL_STRONG_LINE = plotsThresholds.neglog10(plotsThresholds.PVAL_STRONG)
SIG_ALPHA = 0.9
NONSIG_ALPHA = 0.25
DEFAULT_UP_COLOR = plotsPalette.DE_UP
DEFAULT_DOWN_COLOR = plotsPalette.DE_DOWN
NONSIG_COLOR = plotsPalette.DE_NONSIGNIFICANT
MARKER_SHAPES = {'tRNA': 'o', 'nonTRNA': 's'}

# X-axis capping. The axis used to run to 1.1 x the single largest |log2FC|, so one extreme
# feature set the scale for every other. Measured on the hg38 object, 27 of 36
# comparison/readtype combinations left the median |log2FC| occupying under 20% of the
# half-axis; the worst had a median of 1.12 against a max of 8.78 (11.6%).
#
# p95 rather than p99: on that worst case p95 caps at 4.55 and takes the bulk's share to
# 24.6%, while p99 lands near the max and changes almost nothing.
VOLCANO_CAP_PERCENTILE = 95
# Only cap when there is a real tail. A distribution whose max is within this multiple of its
# percentile keeps the old max-based axis, so a tame plot doesn't gain a pointless edge
# marker (vibrChol1's worst comparison maxes at 3.13 against a p95 near 2).
VOLCANO_CAP_ENGAGE_RATIO = 1.5
# Floor, preserving the previous behaviour of squaring small plots up to +-3.
VOLCANO_MIN_HALF_WIDTH = 3.0
# Breathing room so a boundary marker is drawn inside the axes rather than clipped in half.
VOLCANO_AXIS_MARGIN = 1.08
# Off-scale points are pinned to the boundary as triangles pointing the way they ran off.
OFFSCALE_MARKERS = {'right': '>', 'left': '<'}


# Repulsion settings shared by both label passes. The library defaults are tuned for a handful
# of labels; --vollabels defaults to 100. `expand` pads each label's bounding box so neighbours
# push each other apart rather than merely off their own point, and the stronger text force
# lets them actually travel to the space that makes.
_ADJUST_KWARGS = dict(expand=(1.35, 1.6), force_text=(0.4, 0.6), max_move=(20, 20))
LABEL_FONTSIZE = 6
# Below this many pixels of travel a connector line is noise rather than guidance.
MIN_CONNECTOR_PX = 6


def _place_labels(ax, names, x_plot, y):
    """
    Label `names` (already ordered most-interesting first), dropping any that cannot be drawn
    without landing on a better one.

    Repulsion alone does not solve a dense plot: asked to fit 100 labels into a volcano whose
    significant points are packed against the axis limits, adjust_text has nowhere to put them
    and converges on a stack of overlapping text that is worse than no labels at all. So the
    solve runs ONCE, and the labels that still collide afterwards are removed rather than left
    piled up. The user still asks for N labels with --vollabels; what changes is that N is a
    ceiling rather than a promise, and what gets drawn is the subset that is legible.

    Culling has to come last, and the connector lines are therefore drawn here rather than by
    adjust_text: letting it re-solve the survivors would move them again and could reintroduce
    the very overlaps that were just culled. Greedy in priority order, so the features a reader
    most wants named are the ones that survive.
    """
    texts = [ax.text(x_plot[n], y[n], str(n), fontsize=LABEL_FONTSIZE) for n in names]
    # Repel from every point, not just the labelled ones, so a label never lands on a marker.
    adjust_text(texts, x=x_plot.values, y=y.values, ax=ax, **_ADJUST_KWARGS)

    fig = ax.get_figure()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    kept, boxes = [], []
    for name, text in zip(names, texts):
        # Padded, so survivors are separated rather than merely not intersecting.
        box = text.get_window_extent(renderer).expanded(1.1, 1.25)
        if any(box.overlaps(other) for other in boxes):
            text.remove()
            continue
        boxes.append(box)
        kept.append((name, text))

    for name, text in kept:
        anchor = (x_plot[name], y[name])
        # Only connect a label that actually had to travel; a line under a label sitting on its
        # own point is noise. Distance is measured in display space so the threshold means the
        # same thing whatever the axis ranges are.
        moved = np.hypot(*(ax.transData.transform(text.get_position())
                           - ax.transData.transform(anchor)))
        if moved > MIN_CONNECTOR_PX:
            ax.annotate('', xy=anchor, xytext=text.get_position(), textcoords='data',
                        arrowprops={'arrowstyle': '-', 'color': plotsPalette.REFERENCE_LINE,
                                    'lw': 0.5, 'shrinkA': 2, 'shrinkB': 2})

    if len(kept) < len(names):
        logger.info(f'Labelled {len(kept)} of {len(names)} requested markers; the rest could '
                    f'not be placed without overlapping a more significant one.')
    return [name for name, _ in kept]


def resolve_x_limit(values, explicit=None, percentile=VOLCANO_CAP_PERCENTILE,
                    minimum=VOLCANO_MIN_HALF_WIDTH, engage_ratio=VOLCANO_CAP_ENGAGE_RATIO):
    '''
    Half-width of the volcano x-axis, and whether anything falls outside it.

    Returns `(cap, capped)`. `capped` is True when at least one point lies beyond `cap` and
    must therefore be drawn at the boundary rather than at its own position. Nothing is ever
    dropped -- capping changes where a point is drawn, never whether it is.
    '''
    if explicit is not None:
        if explicit <= 0:
            raise ValueError(f'--volxlim must be a positive log2 fold change, got {explicit}.')
        magnitudes = pd.Series(values).abs().dropna()
        return float(explicit), bool((magnitudes > explicit).any())

    magnitudes = pd.Series(values).abs().dropna()
    if magnitudes.empty:
        return float(minimum), False

    largest = float(magnitudes.max())
    candidate = max(float(minimum), float(np.percentile(magnitudes, percentile)))
    if largest <= candidate * engage_ratio:
        # No meaningful tail: keep the historical max-based axis.
        return max(float(minimum), largest), False
    return candidate, True

# The combined overview page always shows BOTH read bases side by side, regardless of
# --diffrts and regardless of --allreads. That is deliberate and is not the inconsistency
# --allreads exists to fix: the defect being avoided is two plots resting on different
# denominators *without saying so*, whereas a labelled side-by-side comparison is the only
# way to see how much transcript-level multi-mapping actually moves a dataset. plotsPca.py
# carries the same constant for the same reason.
OVERVIEW_TRNA_READTYPES = ('nreads_total_unique_norm', 'nreads_total_norm')


def _prettify_readtype(readtype):
    label = readtype.replace('nreads_', '').replace('_norm', '')
    return label.replace('_', ' ').title()


def _resolve_pair_colors(colormap, pair):
    '''
    A point upregulated toward pair[1] takes pair[1]'s configured --colormap color (and
    vice versa for pair[0]), so volcano coloring reuses the same per-group colors as PCA
    instead of introducing a separate config schema.
    '''
    up_color = colormap.get(pair[1], DEFAULT_UP_COLOR) if colormap else DEFAULT_UP_COLOR
    down_color = colormap.get(pair[0], DEFAULT_DOWN_COLOR) if colormap else DEFAULT_DOWN_COLOR
    return up_color, down_color


def _keep_square(standalone, settings):
    """
    Whether this volcano's axes should be forced square.

    Square is the default and stays that way on a COMBINED page whatever the style file says:
    those pages compute their own geometry from the panel count, and a grid of stretched panels
    reads badly. The one exception is a STANDALONE plot the user gave an explicit figsize for --
    there a forced square box would discard the extra width of a deliberately wide figure,
    which bbox_inches='tight' then crops away, so asking for 9.75x6.5 produced something close
    to square. An explicit size is the stronger statement of intent about that one plot.

    Note this is not how PCA and correlation behave: those stay square always, and figsize sets
    how big the square is. A correlation matrix is square because it is symmetric.
    """
    return not (standalone and (settings or {}).get('figsize'))


def _draw_volcano(ax, df_pairs, pair, colormap, toplabels, feature_types, title, show_legend=True,
                  xlim=None, settings=None, standalone=False):
    '''
    Draw a single volcano comparison onto a provided axes. Shared by the standalone individual
    plots and the combined overview page so both use identical styling. `feature_types`, when
    provided, is a Series aligned to df_pairs.index labeling each feature 'tRNA' or 'nonTRNA',
    used to draw tRNA points as circles and non-tRNA points as squares (combined allRNA plots).
    '''
    pair_name = f'{pair[0]}-{pair[1]}'
    log2fc = df_pairs[f'log2_{pair_name}']
    pval = df_pairs[f'pval_{pair_name}']
    # DESeq2 assigns NaN log2FC/padj to genes excluded by independent filtering (low mean count)
    # or Cook's-outlier detection -- these were never given a real significance call (not
    # "insignificant"), so they're dropped rather than plotted. A padj of exactly 0.0 is the
    # opposite case: a float64 underflow on an extremely significant hit, not noise -- clip it
    # away from zero instead of dropping it, so adjust_text's KDTree (which needs every point
    # finite) doesn't lose real data along with a plottable, if extreme, point.
    finite_mask = log2fc.notna() & pval.notna()
    x = log2fc[finite_mask]
    y = -np.log10(pval[finite_mask].clip(lower=np.nextafter(0, 1)))

    settings = settings or {}
    marker_size = settings.get('marker_size') or 40
    # Rasterize only the point layer past the configured threshold: a vector PDF carrying
    # tens of thousands of markers is slow to open and is often rejected on submission,
    # while text and axes stay vector so the figure remains editable.
    threshold = settings.get('rasterize_over')
    rasterized = bool(threshold) and len(df_pairs) > threshold

    up_color, down_color = _resolve_pair_colors(colormap, pair)

    up_mask = (x > LOG2FC_THRESHOLD) & (y > PVAL_SIG_LINE)
    down_mask = (x < -LOG2FC_THRESHOLD) & (y > PVAL_SIG_LINE)
    sig_mask = up_mask | down_mask

    colors = pd.Series([mplcolors.to_rgba(NONSIG_COLOR, alpha=NONSIG_ALPHA)] * len(x), index=x.index)
    colors[up_mask] = pd.Series([mplcolors.to_rgba(up_color, alpha=SIG_ALPHA)] * int(up_mask.sum()), index=x.index[up_mask])
    colors[down_mask] = pd.Series([mplcolors.to_rgba(down_color, alpha=SIG_ALPHA)] * int(down_mask.sum()), index=x.index[down_mask])

    # Cap the axis before drawing, and pin anything beyond it to the boundary. A point that
    # ran off the scale is drawn as a triangle pointing the way it went, so it is visibly an
    # out-of-range marker rather than a real value sitting on the edge.
    cap, capped = resolve_x_limit(x, explicit=xlim)
    x_plot = x.clip(-cap, cap)
    offscale = x.abs() > cap

    if feature_types is not None:
        for ftype, marker in MARKER_SHAPES.items():
            idx = feature_types[feature_types == ftype].index.intersection(x.index)
            idx = idx.difference(x.index[offscale])
            if len(idx) == 0:
                continue
            ax.scatter(x_plot.loc[idx], y.loc[idx], c=list(colors.loc[idx]), marker=marker, s=marker_size, linewidths=0, rasterized=rasterized)
    else:
        keep = x.index[~offscale]
        ax.scatter(x_plot.loc[keep], y.loc[keep], c=list(colors.loc[keep]), marker='o', s=marker_size, linewidths=0, rasterized=rasterized)

    if capped:
        for side, marker in OFFSCALE_MARKERS.items():
            side_mask = offscale & ((x > 0) if side == 'right' else (x < 0))
            idx = x.index[side_mask]
            if len(idx) == 0:
                continue
            ax.scatter(x_plot.loc[idx], y.loc[idx], c=list(colors.loc[idx]), marker=marker,
                       s=marker_size * 1.4, linewidths=0, rasterized=rasterized)

    # Add a line at log2FC = 1.5 and -1.5 and pval = 0.05 and 0.001
    ax.axvline(x=-LOG2FC_THRESHOLD, color=plotsPalette.REFERENCE_LINE, linestyle='--')
    ax.axvline(x=LOG2FC_THRESHOLD, color=plotsPalette.REFERENCE_LINE, linestyle='--')
    ax.axhline(y=PVAL_SIG_LINE, color=plotsPalette.REFERENCE_LINE, linestyle='--')
    ax.axhline(y=PVAL_STRONG_LINE, color=plotsPalette.REFERENCE_LINE, linestyle='--')

    # Set axis limits so that the plot is square. The half-width comes from resolve_x_limit()
    # above rather than from max(abs(x)), so one extreme feature no longer sets the scale for
    # every other one.
    ax.set_xlim(-cap * VOLCANO_AXIS_MARGIN, cap * VOLCANO_AXIS_MARGIN)
    if ax.get_ylim()[1] < 10:
        ax.set_ylim(0, 10)

    # Label markers of interest: significant points, most extreme (|log2FC| * -log10(pval)) first.
    # toplabels=None labels all significant points (unbounded cost -- the CLI default is 100
    # instead, see --vollabels), 0 disables labels, N labels the top N.
    if toplabels != 0:
        score = (x.abs() * y)[sig_mask].sort_values(ascending=False)
        label_idx = score.index if toplabels is None else score.index[:toplabels]
        if len(label_idx) > 0:
            # Labels follow the drawn (clipped) position, not the raw one -- otherwise a
            # capped feature's label is placed off-canvas and adjust_text drags a connector
            # line out of the axes chasing it.
            _place_labels(ax, label_idx, x_plot, y)

    ax.set_xlabel('log2(fold change)')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title(title)
    if _keep_square(standalone, settings):
        ax.set_box_aspect(1)

    if show_legend:
        legend_elements = [
            Line2D([0], [0], marker='o', linestyle='None', markerfacecolor=up_color, markeredgecolor='none', markersize=8, label=f'Up in {pair[1]}'),
            Line2D([0], [0], marker='o', linestyle='None', markerfacecolor=down_color, markeredgecolor='none', markersize=8, label=f'Up in {pair[0]}'),
            Line2D([0], [0], marker='o', linestyle='None', markerfacecolor=NONSIG_COLOR, markeredgecolor='none', alpha=NONSIG_ALPHA, markersize=8, label='Not significant'),
        ]
        if feature_types is not None:
            legend_elements += [
                Line2D([0], [0], marker='o', linestyle='None', markerfacecolor=plotsPalette.LEGEND_SHAPE_SWATCH, markeredgecolor='none', markersize=8, label='tRNA'),
                Line2D([0], [0], marker='s', linestyle='None', markerfacecolor=plotsPalette.LEGEND_SHAPE_SWATCH, markeredgecolor='none', markersize=8, label='Non-tRNA'),
            ]
        ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)


def _save_individual_volcano(df_pairs, pair, cutoff, colormap, toplabels, output, basename, variant_label, feature_types, threaded, xlim=None, settings=None):
    '''
    Draw and save one comparison as its own standalone PDF in the individual/ subfolder.
    '''
    pair_name = f'{pair[0]}-{pair[1]}'
    title = f'{pair[0]} vs {pair[1]} ({variant_label})'

    fig = plt.figure(figsize=toolsTG.figsize_for(settings, (8, 8)))
    ax = fig.gca()
    _draw_volcano(ax, df_pairs, pair, colormap, toplabels, feature_types, title, xlim=xlim, settings=settings, standalone=True)

    filename = f'{output}{basename}_{pair_name}_{cutoff}_volcano.pdf'
    if threaded:
        threaded += f'Saving figure: {filename}\n'
    else:
        logger.info(f'Saving figure: {filename}')
    plt.savefig(filename, bbox_inches='tight')
    plt.close(fig)
    return threaded


def _save_combined_volcano_page(pdf, pair, colormap, toplabels, slots, xlim=None, settings=None):
    '''
    Draw one page of the combined overview PDF for a single comparison pair. `slots` is a list
    of (variant_label, df_pairs, feature_types) tuples: 2 slots (tRNA total_unique + total) are
    stacked vertically, and 4 slots (the above plus non-tRNA + combined) form a 2x2 grid.
    '''
    n = len(slots)
    if n <= 2:
        fig, axs = plt.subplots(2, 1, figsize=(8, 17))
    else:
        fig, axs = plt.subplots(2, 2, figsize=(18, 17))
    axs = axs.flatten()

    for ax, (variant_label, df_pairs, feature_types) in zip(axs, slots):
        _draw_volcano(ax, df_pairs, pair, colormap, toplabels, feature_types, title=variant_label, show_legend=True, xlim=xlim, settings=settings)
    for ax in axs[n:]:
        ax.set_visible(False)
    fig.subplots_adjust(wspace=0.5, hspace=0.3)

    fig.suptitle(f'{pair[0]} vs {pair[1]}', fontsize=14, y=0.925)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def visualizer(adata, grp, readtypes, cutoff, output, colormap=None, toplabels=100, threaded=True, config_name='default', overwrite=False, is_full_variant=True, xlim=None, settings=None):
    '''
    Generate volcano visualizations for each pairwise group comparison in an AnnData object.

    Every individual plot (one per --diffrts read type / non-tRNA / combined family, per pair)
    is saved to an `individual/` subfolder, mirroring plotsCoverage.py's split-vs-combined output
    layout. A combined overview PDF is also saved at the top level with one page per comparison
    pair: 2 stacked plots (tRNA at total_unique + total, mirroring PCA's default --pcareadtypes)
    or a 2x2 grid of 4 (the same two plus non-tRNA and combined) when non-tRNA data is available.

    tRNA volcano plots use tRNAgraph's default tRNA/tRX-controlled DESeq2 normalization
    (adata.obs 'nreads_*_norm', matching adata.X) via toolsTG.adataLog2FC. The non-tRNA volcano
    and the combined tRNA + non-tRNA volcano instead use the all-feature-controlled
    DESeq2 normalization (adata.uns['nontRNA_counts'] and adata.uns['deseq2_sizefactors_allfeatures']),
    since tRNA-controlled size factors are not representative of non-tRNA library composition.
    This mirrors the naming/normalization/gating pattern used by plotsPca.py.
    '''
    individual_output = f'{output}individual/'
    if threaded:
        threaded += toolsTG.builder(individual_output) + '\n'
    else:
        logger.info(toolsTG.builder(individual_output))

    if colormap:
        colormap = {k: v if v[0] != '#' else mplcolors.to_rgb(v) for k, v in colormap.items()}

    trna_dfs = {}
    pairs = None
    for rt in readtypes:
        # `readtypes` arrives as resolved obs column names from adataGraph.resolved_diffrts().
        df, log2fc_dict = toolsTG.adataLog2FC(adata, grp, rt, readcount_cutoff=cutoff, config_name=config_name, overwrite=overwrite).main()
        trna_dfs[rt] = df
        pairs = log2fc_dict[config_name][grp]['pairs']
        for pair in pairs:
            threaded = _save_individual_volcano(df, pair, cutoff, colormap, toplabels, individual_output, basename=f'tRNA_{rt}', variant_label=_prettify_readtype(rt), feature_types=None, threaded=threaded, xlim=xlim, settings=settings)

    # The combined overview always shows total_unique + total for tRNA, regardless of --diffrts.
    for rt in OVERVIEW_TRNA_READTYPES:
        if rt not in trna_dfs:
            df, log2fc_dict = toolsTG.adataLog2FC(adata, grp, rt, readcount_cutoff=cutoff, config_name=config_name, overwrite=overwrite).main()
            trna_dfs[rt] = df
            pairs = log2fc_dict[config_name][grp]['pairs']

    # Non-tRNA and combined tRNA + non-tRNA volcano plots. These run automatically alongside the
    # per-readtype plots above (no separate flag), but only when non-tRNA feature counts are
    # available -- adata.uns['nontRNA_counts'] is only populated with real rows when
    # `trnagraph analyze build` was given a --gtf; otherwise it's present but empty.
    nontrna_pairs_df = None
    combined_pairs_df = None
    combined_feature_types = None

    nontrna_df, skip_message = toolsTG.resolve_nontrna_counts(adata, is_full_variant, 'and combined volcano plots')
    if nontrna_df is None:
        if threaded:
            threaded += skip_message + '\n'
        else:
            logger.warning(skip_message)
    else:
        sample_group_map = dict(zip(adata.obs['sample'].astype(str), adata.obs[grp]))

        # Non-tRNA-only volcano. adata.uns['nontRNA_counts'] is normalized against the
        # all-feature-controlled DESeq2 size factors (see adataBuild.py), which is the
        # statistically appropriate normalization for non-tRNA feature counts.
        nontrna_pairs_df, nontrna_pairs = toolsTG.log2fc_from_wide_df(nontrna_df, sample_group_map, readcount_cutoff=cutoff)
        for pair in nontrna_pairs:
            threaded = _save_individual_volcano(nontrna_pairs_df, pair, cutoff, colormap, toplabels, individual_output, basename='nontRNA', variant_label='Non-tRNA RNAs', feature_types=None, threaded=threaded, xlim=xlim, settings=settings)

        # Combined volcano: all tRNA reads (total, not unique-only) + non-tRNA RNAs, both normalized
        # against the all-feature-controlled size factors so the two feature sets are on the same
        # scale. This requires re-deriving tRNA total counts from raw counts, since adata.obs only
        # stores the tRNA/tRX-controlled (default) normalization.
        allfeature_sizefactors = adata.uns.get('deseq2_sizefactors_allfeatures')
        if allfeature_sizefactors is None or 'nreads_total_raw' not in adata.obs.columns:
            msg = ('No all-feature-controlled DESeq2 size factors found in AnnData object '
                   '(uns[\'deseq2_sizefactors_allfeatures\'] missing, or obs[\'nreads_total_raw\'] '
                   'unavailable). Skipping combined tRNA + non-tRNA volcano plot. Rebuild with the '
                   'current `trnagraph analyze build` to enable this plot.')
            if threaded:
                threaded += msg + '\n'
            else:
                logger.warning(msg)
        else:
            sf_map = {k: float(v) for k, v in dict(allfeature_sizefactors).items()}
            total_df = pd.DataFrame(adata.obs, columns=['trna', 'sample', 'nreads_total_raw'])
            missing_samples = sorted(set(total_df['sample'].astype(str)) - set(sf_map.keys()))
            if missing_samples:
                msg = f'All-feature size factors not found for samples {missing_samples}; defaulting to 1.0 for the combined volcano plot.'
                if threaded:
                    threaded += msg + '\n'
                else:
                    logger.warning(msg)
            total_df['nreads_total_norm_allfeatures'] = total_df['nreads_total_raw'] / total_df['sample'].astype(str).map(lambda s: sf_map.get(s, 1.0))
            total_df = total_df.pivot_table(index='trna', columns='sample', values='nreads_total_norm_allfeatures', observed=True)

            shared_cols = total_df.columns.intersection(nontrna_df.columns)
            if len(shared_cols) < len(total_df.columns) or len(shared_cols) < len(nontrna_df.columns):
                missing_from_nontrna = set(total_df.columns) - set(nontrna_df.columns)
                missing_from_total = set(nontrna_df.columns) - set(total_df.columns)
                msg = f'Sample columns did not fully match between tRNA and non-tRNA counts for combined volcano; using intersection of {len(shared_cols)} shared samples.'
                if missing_from_nontrna:
                    msg += f' Missing from nontRNA_counts: {sorted(missing_from_nontrna)}.'
                if missing_from_total:
                    msg += f' Missing from tRNA total counts: {sorted(missing_from_total)}.'
                if threaded:
                    threaded += msg + '\n'
                else:
                    logger.warning(msg)

            if len(shared_cols) == 0:
                msg = 'No overlapping sample columns between tRNA and non-tRNA counts. Skipping combined volcano plot.'
                if threaded:
                    threaded += msg + '\n'
                else:
                    logger.warning(msg)
            else:
                combined_df = pd.concat([total_df[shared_cols], nontrna_df[shared_cols]], axis=0)
                combined_feature_types = pd.Series(
                    ['tRNA' if idx in total_df.index else 'nonTRNA' for idx in combined_df.index],
                    index=combined_df.index
                )
                combined_pairs_df, combined_pairs = toolsTG.log2fc_from_wide_df(combined_df, sample_group_map, readcount_cutoff=cutoff)
                for pair in combined_pairs:
                    threaded = _save_individual_volcano(combined_pairs_df, pair, cutoff, colormap, toplabels, individual_output, basename='allRNA', variant_label='All tRNA + Non-tRNA RNAs', feature_types=combined_feature_types, threaded=threaded, xlim=xlim, settings=settings)

    # Combined overview: one PDF, one page per comparison pair.
    if pairs:
        overview_path = f'{output}{grp}_combined_{cutoff}_volcano.pdf'
        with PdfPages(overview_path) as pdf:
            for pair in pairs:
                slots = [
                    ('tRNA Unique', trna_dfs[OVERVIEW_TRNA_READTYPES[0]], None),
                    ('tRNA', trna_dfs[OVERVIEW_TRNA_READTYPES[1]], None),
                ]
                if nontrna_pairs_df is not None:
                    slots.append(('Other RNA', nontrna_pairs_df, None))
                if combined_pairs_df is not None:
                    slots.append(('All RNA', combined_pairs_df, combined_feature_types))
                _save_combined_volcano_page(pdf, pair, colormap, toplabels, slots, xlim=xlim, settings=settings)
        if threaded:
            threaded += f'Saving combined volcano overview: {overview_path}\n'
        else:
            logger.info(f'Saving combined volcano overview: {overview_path}')

    if threaded:
        return threaded

if __name__ == '__main__':
    pass
