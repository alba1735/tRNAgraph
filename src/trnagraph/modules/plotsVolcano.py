#!/usr/bin/env python3

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from adjustText import adjust_text

from . import toolsTG

plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Significance boundaries. LOG2FC_THRESHOLD + PVAL_SIG_LINE together define which points count
# as "significant" for coloring/opacity/labeling; PVAL_STRONG_LINE is drawn as an additional
# (non-classifying) reference line, matching the pre-existing dashed-line behavior.
LOG2FC_THRESHOLD = 1.5
PVAL_SIG_LINE = 1.3
PVAL_STRONG_LINE = 3
SIG_ALPHA = 0.9
NONSIG_ALPHA = 0.25
DEFAULT_UP_COLOR = '#d62728'
DEFAULT_DOWN_COLOR = '#1f77b4'
NONSIG_COLOR = '#7f7f7f'
MARKER_SHAPES = {'tRNA': 'o', 'nonTRNA': 's'}

# The combined overview page always uses these two tRNA read types, mirroring PCA's default
# --pcareadtypes (total_unique + total), regardless of what --diffrts requests for the
# individual per-readtype plots.
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


def _draw_volcano(ax, df_pairs, pair, colormap, toplabels, feature_types, title, show_legend=True):
    '''
    Draw a single volcano comparison onto a provided axes. Shared by the standalone individual
    plots and the combined overview page so both use identical styling. `feature_types`, when
    provided, is a Series aligned to df_pairs.index labeling each feature 'tRNA' or 'nonTRNA',
    used to draw tRNA points as circles and non-tRNA points as squares (combined allRNA plots).
    '''
    pair_name = f'{pair[0]}-{pair[1]}'
    x = df_pairs[f'log2_{pair_name}']
    y = -np.log10(df_pairs[f'pval_{pair_name}'])

    up_color, down_color = _resolve_pair_colors(colormap, pair)

    up_mask = (x > LOG2FC_THRESHOLD) & (y > PVAL_SIG_LINE)
    down_mask = (x < -LOG2FC_THRESHOLD) & (y > PVAL_SIG_LINE)
    sig_mask = up_mask | down_mask

    colors = pd.Series([mplcolors.to_rgba(NONSIG_COLOR, alpha=NONSIG_ALPHA)] * len(x), index=x.index)
    colors[up_mask] = pd.Series([mplcolors.to_rgba(up_color, alpha=SIG_ALPHA)] * int(up_mask.sum()), index=x.index[up_mask])
    colors[down_mask] = pd.Series([mplcolors.to_rgba(down_color, alpha=SIG_ALPHA)] * int(down_mask.sum()), index=x.index[down_mask])

    if feature_types is not None:
        for ftype, marker in MARKER_SHAPES.items():
            idx = feature_types[feature_types == ftype].index.intersection(x.index)
            if len(idx) == 0:
                continue
            ax.scatter(x.loc[idx], y.loc[idx], c=list(colors.loc[idx]), marker=marker, s=40, linewidths=0)
    else:
        ax.scatter(x, y, c=list(colors), marker='o', s=40, linewidths=0)

    # Add a line at log2FC = 1.5 and -1.5 and pval = 0.05 and 0.001
    ax.axvline(x=-LOG2FC_THRESHOLD, color='black', linestyle='--')
    ax.axvline(x=LOG2FC_THRESHOLD, color='black', linestyle='--')
    ax.axhline(y=PVAL_SIG_LINE, color='black', linestyle='--')
    ax.axhline(y=PVAL_STRONG_LINE, color='black', linestyle='--')

    # Set axis limits so that the plot is square
    if ax.get_xlim()[1] < 3 and ax.get_xlim()[0] > -3:
        ax.set_xlim(-3, 3)
    else:
        ax.set_xlim(-1.1 * max(abs(x)), 1.1 * max(abs(x)))
    if ax.get_ylim()[1] < 10:
        ax.set_ylim(0, 10)

    # Label markers of interest: significant points, most extreme (|log2FC| * -log10(pval)) first.
    # toplabels=None labels all significant points, 0 disables labels, N labels the top N.
    if toplabels != 0:
        score = (x.abs() * y)[sig_mask].sort_values(ascending=False)
        label_idx = score.index if toplabels is None else score.index[:toplabels]
        if len(label_idx) > 0:
            texts = [ax.text(x[name], y[name], str(name), fontsize=6) for name in label_idx]
            # Repel labels from every point (not just labeled ones) and from each other; only
            # draws a connector line back to its point if the label actually had to move.
            adjust_text(texts, x=x.values, y=y.values, ax=ax, arrowprops={'arrowstyle': '-', 'color': 'black', 'lw': 0.5})

    ax.set_xlabel('log2(fold change)')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title(title)
    ax.set_box_aspect(1)

    if show_legend:
        legend_elements = [
            Line2D([0], [0], marker='o', linestyle='None', markerfacecolor=up_color, markeredgecolor='none', markersize=8, label=f'Up in {pair[1]}'),
            Line2D([0], [0], marker='o', linestyle='None', markerfacecolor=down_color, markeredgecolor='none', markersize=8, label=f'Up in {pair[0]}'),
            Line2D([0], [0], marker='o', linestyle='None', markerfacecolor=NONSIG_COLOR, markeredgecolor='none', alpha=NONSIG_ALPHA, markersize=8, label='Not significant'),
        ]
        if feature_types is not None:
            legend_elements += [
                Line2D([0], [0], marker='o', linestyle='None', markerfacecolor='grey', markeredgecolor='none', markersize=8, label='tRNA'),
                Line2D([0], [0], marker='s', linestyle='None', markerfacecolor='grey', markeredgecolor='none', markersize=8, label='Non-tRNA'),
            ]
        ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)


def _save_individual_volcano(df_pairs, pair, cutoff, colormap, toplabels, output, basename, variant_label, feature_types, threaded):
    '''
    Draw and save one comparison as its own standalone PDF in the individual/ subfolder.
    '''
    pair_name = f'{pair[0]}-{pair[1]}'
    title = f'{pair[0]} vs {pair[1]} ({variant_label})'

    fig = plt.figure(figsize=(8, 8))
    ax = fig.gca()
    _draw_volcano(ax, df_pairs, pair, colormap, toplabels, feature_types, title)

    filename = f'{output}{basename}_{pair_name}_{cutoff}_volcano.pdf'
    if threaded:
        threaded += f'Saving figure: {filename}\n'
    else:
        print(f'Saving figure: {filename}')
    plt.savefig(filename, bbox_inches='tight')
    plt.close(fig)
    return threaded


def _save_combined_volcano_page(pdf, pair, colormap, toplabels, slots):
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
        _draw_volcano(ax, df_pairs, pair, colormap, toplabels, feature_types, title=variant_label, show_legend=True)
    for ax in axs[n:]:
        ax.set_visible(False)
    fig.subplots_adjust(wspace=0.5, hspace=0.3)

    fig.suptitle(f'{pair[0]} vs {pair[1]}', fontsize=14, y=0.925)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def visualizer(adata, grp, readtypes, cutoff, output, colormap=None, toplabels=None, threaded=True, config_name='default', overwrite=False):
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
        print(toolsTG.builder(individual_output))

    if colormap:
        colormap = {k: v if v[0] != '#' else mplcolors.to_rgb(v) for k, v in colormap.items()}

    trna_dfs = {}
    pairs = None
    for readtype in readtypes:
        rt = f'nreads_{readtype}_norm'
        df, log2fc_dict = toolsTG.adataLog2FC(adata, grp, rt, readcount_cutoff=cutoff, config_name=config_name, overwrite=overwrite).main()
        trna_dfs[rt] = df
        pairs = log2fc_dict[config_name][grp]['pairs']
        for pair in pairs:
            threaded = _save_individual_volcano(df, pair, cutoff, colormap, toplabels, individual_output, basename=f'tRNA_{rt}', variant_label=_prettify_readtype(rt), feature_types=None, threaded=threaded)

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

    nontrna_df = adata.uns.get('nontRNA_counts')
    if nontrna_df is None or nontrna_df.empty:
        msg = ('No non-tRNA feature counts found in AnnData object '
               '(uns[\'nontRNA_counts\'] missing or empty). Skipping non-tRNA and combined volcano plots. '
               'Re-run `trnagraph analyze build` with --gtf to enable these plots.')
        if threaded:
            threaded += msg + '\n'
        else:
            print(msg)
    else:
        sample_group_map = dict(zip(adata.obs['sample'].astype(str), adata.obs[grp]))

        # Non-tRNA-only volcano. adata.uns['nontRNA_counts'] is normalized against the
        # all-feature-controlled DESeq2 size factors (see adataBuild.py), which is the
        # statistically appropriate normalization for non-tRNA feature counts.
        nontrna_pairs_df, nontrna_pairs = toolsTG.log2fc_from_wide_df(nontrna_df, sample_group_map, readcount_cutoff=cutoff)
        for pair in nontrna_pairs:
            threaded = _save_individual_volcano(nontrna_pairs_df, pair, cutoff, colormap, toplabels, individual_output, basename='nontRNA', variant_label='Non-tRNA RNAs', feature_types=None, threaded=threaded)

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
                print(msg)
        else:
            sf_map = {k: float(v) for k, v in dict(allfeature_sizefactors).items()}
            total_df = pd.DataFrame(adata.obs, columns=['trna', 'sample', 'nreads_total_raw'])
            missing_samples = sorted(set(total_df['sample'].astype(str)) - set(sf_map.keys()))
            if missing_samples:
                msg = f'All-feature size factors not found for samples {missing_samples}; defaulting to 1.0 for the combined volcano plot.'
                if threaded:
                    threaded += msg + '\n'
                else:
                    print(msg)
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
                    print(msg)

            if len(shared_cols) == 0:
                msg = 'No overlapping sample columns between tRNA and non-tRNA counts. Skipping combined volcano plot.'
                if threaded:
                    threaded += msg + '\n'
                else:
                    print(msg)
            else:
                combined_df = pd.concat([total_df[shared_cols], nontrna_df[shared_cols]], axis=0)
                combined_feature_types = pd.Series(
                    ['tRNA' if idx in total_df.index else 'nonTRNA' for idx in combined_df.index],
                    index=combined_df.index
                )
                combined_pairs_df, combined_pairs = toolsTG.log2fc_from_wide_df(combined_df, sample_group_map, readcount_cutoff=cutoff)
                for pair in combined_pairs:
                    threaded = _save_individual_volcano(combined_pairs_df, pair, cutoff, colormap, toplabels, individual_output, basename='allRNA', variant_label='All tRNA + Non-tRNA RNAs', feature_types=combined_feature_types, threaded=threaded)

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
                _save_combined_volcano_page(pdf, pair, colormap, toplabels, slots)
        if threaded:
            threaded += f'Saving combined volcano overview: {overview_path}\n'
        else:
            print(f'Saving combined volcano overview: {overview_path}')

    if threaded:
        return threaded

if __name__ == '__main__':
    pass
