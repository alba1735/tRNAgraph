#!/usr/bin/env python3

import logging

import seaborn as sns
import numpy as np
import pandas as pd
import os

from . import toolsTG
from . import plotsPalette

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

logger = logging.getLogger(__name__)


def visualizer(adata, grp, readtypes, cutoff, heatbound, heatsubplots, output, threaded=False, config_name='default', overwrite=False, settings=None):
    '''
    Generate heatmap visualizations for each group in an AnnData object. `readtypes` are
    fully-resolved obs column names (e.g. 'nreads_total_unique_norm'), already carrying the
    read basis chosen once for the whole `graph` command -- see toolsTG.resolve_readtype().

    The combined,
    multi-comparison heatmap PDFs are saved at the top level; when --heatsubplots is set, the
    individual per-comparison heatmaps are also saved to an `individual/` subfolder, mirroring
    plotsCoverage.py's split-vs-combined output layout.
    '''
    grp = toolsTG.resolve_grp_column(adata, grp, 'group')

    individual_output = f'{output}individual/'
    if heatsubplots:
        logger.info(toolsTG.builder(individual_output))

    # Create an empty df to store the log2FC values for each group so a combined heatmap can be generated as well
    df_combine = pd.DataFrame()
    # Create a heatmap for each group
    for readtype in readtypes:
        # `readtypes` arrives as resolved obs column names -- adataGraph.resolved_diffrts()
        # applies the command-wide read basis once, so no module decides its own denominator.
        # Create a color palette for the heatmap
        cmap = plotsPalette.gradient(settings, 'lfc')
        # The significance panel sits beside the log2FC panel in the same figure, so it
        # is resolved here and threaded through rather than read inside heatmap_plot.
        sig_cmap = plotsPalette.gradient(settings, 'significance')
        # Create a correlation matrix from reads stored in adata observations
        df, log2fc_dict = toolsTG.adataLog2FC(adata, grp, readtype, readcount_cutoff=cutoff, config_name=config_name, overwrite=overwrite).main()
        if df.empty:
            if threaded:
                threaded += f'No data for {readtype} heatmap.\n'
            else:
                logger.warning(f'No data for {readtype} heatmap.')
            continue
        df['readtype'] = readtype
        # combine df with df_combine by stacking them vertically if readtype is not total_unique or total
        if readtype not in ('nreads_total_unique_norm', 'nreads_total_norm'):
            df_combine = pd.concat([df_combine, df], axis=0)
        # save df to csv
        csv_output = output.replace('graphs', 'results')
        os.makedirs(csv_output, exist_ok=True)
        df.to_csv(f'{csv_output}{grp}_{readtype}_{cutoff}_{heatbound}_heatmap.csv')
        # Create a pdf with a heatmap for sorted by each group on each page
        with PdfPages(f'{output}{grp}_{readtype}_{cutoff}_{heatbound}_heatmap.pdf') as pdf:
            for col in [i for i in df.columns.tolist() if 'log2' in i]:
                plt = heatmap_plot(df, col, cmap, sig_cmap, heatbound, settings)
                # Save figure
                if threaded:
                    threaded += f'Saving heatmap for {readtype} {col}...\n'
                else:
                    logger.info(f'Saving heatmap for {readtype} {col}...')
                if heatsubplots:
                    plt.savefig(f'{individual_output}{grp}_{readtype}_{cutoff}_{heatbound}_{col}_heatmap.pdf', bbox_inches='tight')
                pdf.savefig(bbox_inches='tight')
                plt.close()
    # Create a heatmap for the combined groups
    if not df_combine.empty:
        with PdfPages(f'{output}{grp}_combine_{cutoff}_{heatbound}_heatmap.pdf') as pdf:
            for col in [i for i in df_combine.columns.tolist() if 'log2' in i]:
                plt = heatmap_plot(df_combine, col, cmap, sig_cmap, heatbound, settings)
                # Save figure
                if threaded:
                    threaded += f'Saving heatmap for combine {col}...\n'
                else:
                    logger.info(f'Saving heatmap for combine {col}...')
                if heatsubplots:
                    plt.savefig(f'{individual_output}{grp}_combine_{cutoff}_{heatbound}_{col}_heatmap.pdf', bbox_inches='tight')
                pdf.savefig(bbox_inches='tight')
                plt.close()
    if threaded:
        return threaded

def heatmap_plot(df, col, cmap, sig_cmap, heatbound, settings=None):
    # Sort df by the sum of the log2FC values for each row
    tdf = df.sort_values(by=col, ascending=False)
    # subset df to only include the top and bottom heatbound values only if the heatmap is larger than heatbound
    if len(tdf) > int(heatbound):
        tdf = pd.concat([tdf.iloc[:int(heatbound), :], tdf.iloc[-int(heatbound):, :]])  
    # Create a heatmap
    # Width follows the comparison count; the height is a default the style file may replace.
    default_size = (16, 12) if len(tdf.columns) >= 20 else (6, 12)
    fig, axs = plt.subplots(1, 2, figsize=toolsTG.figsize_for(settings, default_size))

    log_tdf = tdf[[i for i in tdf.columns if 'log2' in i]]
    pval_tdf = tdf[[i for i in tdf.columns if 'pval' in i]]
    sns.heatmap(log_tdf, ax=axs[0], cmap=cmap, center=0, vmax=4, vmin=-4, cbar=True, square=True, cbar_kws={'fraction':0.05, 'pad':0.05})
    sns.heatmap(-np.log10(pval_tdf), ax=axs[1], cmap=sig_cmap, vmin=0, vmax=3, cbar=True, yticklabels=False, square=True, cbar_kws={'fraction':0.05, 'pad':0.05})
    axs[0].tick_params(axis='x', labelrotation=90)
    axs[1].tick_params(axis='x', labelrotation=90)
    axs[1].set_ylabel('')

    # # Set a column of symbols to the left of the heatmap based on readtype
    # for i, row in enumerate(tdf.index.tolist()):
    #     if tdf.loc[row, 'readtype'] == 'nreads_fiveprime_unique_norm' or tdf.loc[row, 'readtype'] == 'nreads_fiveprime_norm':
    #         axs[0].text(-0.5, i+0.5, '\u25CF', fontsize=8, horizontalalignment='center', verticalalignment='center')
    #     elif tdf.loc[row, 'readtype'] == 'nreads_threeprime_unique_norm' or tdf.loc[row, 'readtype'] == 'nreads_threeprime_norm':
    #         axs[0].text(-0.5, i+0.5, '\u25A0', fontsize=8, horizontalalignment='center', verticalalignment='center')
    #     elif tdf.loc[row, 'readtype'] == 'nreads_other_unique_norm' or tdf.loc[row, 'readtype'] == 'nreads_other_norm':
    #         axs[0].text(-0.5, i+0.5, '\u25B2', fontsize=8, horizontalalignment='center', verticalalignment='center')
    #     elif tdf.loc[row, 'readtype'] == 'nreads_whole_unique_norm' or tdf.loc[row, 'readtype'] == 'nreads_whole_norm':
    #         axs[0].text(-0.5, i+0.5, '\u25C6', fontsize=8, horizontalalignment='center', verticalalignment='center')
    #     elif tdf.loc[row, 'readtype'] == 'nreads_precounts_norm':
    #         axs[0].text(-0.5, i+0.5, '\u2715', fontsize=8, horizontalalignment='center', verticalalignment='center')
    # # Adjust xticklabels to the left of the heatmap to make room for the symbols
    # axs[0].set_xticklabels(['']+tdf.columns.tolist(), fontsize=8)

    cbar = axs[0].collections[0].colorbar
    # cbar.set_ticks([-2, -1, 0, 1, 2])
    # cbar.set_ticklabels(['<=-2', '-1', '0', '1', '>=2'])
    cbar.set_ticks([-4,-3,-2,-1,0,1,2,3,4])
    cbar.set_ticklabels(['<=-4','-3','-2','-1','0','1','2','3','>=4'])
    cbar.ax.tick_params(labelsize=8) 
    cbar = axs[1].collections[0].colorbar
    cbar.set_ticks([0, 1.3, 3])
    cbar.set_ticklabels(['1', '0.05', '<=0.001'])
    cbar.ax.tick_params(labelsize=8)

    # Panel titles first: the figure title is positioned relative to them below.
    axs[0].set_title('log2FC', fontsize=10)
    axs[1].set_title('pval', fontsize=10)
    # Anchor the figure title to where the panels ACTUALLY end, not to a fixed height.
    # `square=True` sizes the axes from the row count, so in a fixed-height figure a short
    # heatmap leaves a tall band of interior white that bbox_inches='tight' cannot trim,
    # because the title pins the top edge. Measured before this change: the gap ran to 44% of
    # the saved image at 4 rows against 5% at 20 -- exactly the row-count dependence this
    # removes. get_position() already reflects the aspect adjustment; the tight bbox is
    # preferred where a renderer is available because it also covers the panel titles.
    try:
        renderer = fig.canvas.get_renderer()
        inv = fig.transFigure.inverted()
        top = max(a.get_tightbbox(renderer).transformed(inv).y1 for a in axs)
    except AttributeError:
        # Backends without get_renderer: fall back to the axes box plus room for the titles.
        top = max(a.get_position().y1 for a in axs) + 0.03
    plt.suptitle(f'Heatmap of log2FC of tRNA read counts between groups \nsorted by {col}',
                 fontsize=12, y=min(top + 0.01, 0.999), va='bottom')

    return plt

if __name__ == '__main__':
    pass
