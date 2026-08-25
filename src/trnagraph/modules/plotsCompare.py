#!/usr/bin/env python3

import logging

import numpy as np

from . import toolsTG

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns

logger = logging.getLogger(__name__)

def visualizer(adata, comparegrp1, comparegrp2, colormap, output, threaded=True, read_basis=toolsTG.READ_BASIS_UNIQUE):
    # Fall back to 'sample' (with a warning) if the specified columns aren't in the AnnData object
    comparegrp1 = toolsTG.resolve_grp_column(adata, comparegrp1, 'comparegrp1')
    comparegrp2 = toolsTG.resolve_grp_column(adata, comparegrp2, 'comparegrp2')
    # Create a color palette for the p
    if colormap != None:
        pal = {k:v if v[0]!='#' else mplcolors.to_rgb(v) for k,v in colormap.items()}
    else:
        pal = sns.husl_palette(len(adata.obs[comparegrp1].unique()))
        pal = dict(zip(sorted(adata.obs[comparegrp1].unique()), pal))

    # Was hardcoded to all reads while plotsCluster was hardcoded to unique, so two plots of
    # one dataset rested on different denominators with nothing on either saying so.
    readtype = toolsTG.resolve_readtype('total', read_basis, adata)

    for countgrp in ['amino','iso']:
        # Get log2 fold change dataframe from analysis_tools
        df = toolsTG.log2fc_compare_df(adata, countgrp, [comparegrp1, comparegrp2], readtype, 0)
        # log2fc_compare_df can legitimately produce zero valid comparegrp1/comparegrp2 pairs
        # (e.g. comparegrp1 has no comparegrp2 value shared across every one of its own values --
        # notably always true when comparegrp1 is a per-observation-unique column like 'sample',
        # the Parameter Fallback default) -- its 'log2'/'pval' column levels don't exist at all
        # in that case, so skip this countgrp instead of a KeyError on the .loc[:, ('log2')] below.
        if 'log2' not in df.columns.get_level_values('stats'):
            logger.warning(
                f'WARNING: no valid {comparegrp1}/{comparegrp2} pairs found for {countgrp}; skipping compare plot.'
            )
            continue
        # Sort the df by the mean of the log2 fold change
        df = df.loc[df.loc[:, ('log2')].abs().mean(axis=1).sort_values(ascending=True).index, :]
        # create a stacked horizontal bar plot for each group 
        fig, ax = plt.subplots(figsize=(8, 12))
        cgrp1list = df.columns.get_level_values('cgrp1').unique()
        cgrp2list = df.columns.get_level_values('cgrp2').unique()
        # Create a bar widths table for amount of values in cgrp1
        barwidths = np.linspace(0.1, 0.9, len(cgrp1list)+1)
        bardiff = np.diff(barwidths)[0]
        barwidths = {k:v for k,v in zip(cgrp1list, barwidths)}
        # Enumerate the index of the dataframe so that the bar heights can be adjusted
        en_dict = dict(enumerate(df.index.values))
        for cgrp2 in cgrp2list:
            xminmax = tuple([-np.abs(df.loc[:, ('log2')]).max().max()*1.1, np.abs(df.loc[:, ('log2')]).max().max()*1.1])
            for cgrp1 in cgrp1list:
                for y,posname in en_dict.items():
                    if abs(df.loc[posname, ('log2',cgrp1,cgrp2)]) >= 1:
                        ax.barh(y+barwidths[cgrp1], df.loc[posname, ('log2',cgrp1,cgrp2)], color=pal[cgrp1], align='edge',
                                height=bardiff, linewidth=bardiff, edgecolor=pal[cgrp1], label=cgrp1)
                        # ax.hlines(y+0.25, xmin=xminmax[0], xmax=xminmax[1], color='lightgray', linewidth=bardiff, zorder=-1)
                        # ax.hlines(y+0.75, xmin=xminmax[0], xmax=xminmax[1], color='lightgray', linewidth=bardiff, zorder=-1)
                        ax.barh(y+0.4, xminmax, color='lightgray', align='edge', height=0.2, linewidth=0, zorder=-2)
                    else:
                        ax.barh(y+barwidths[cgrp1], df.loc[posname, ('log2',cgrp1,cgrp2)], color='white', align='edge',
                                height=bardiff, linewidth=bardiff, edgecolor=pal[cgrp1], label=cgrp1)
                    ax.hlines(y+0.5, xmin=xminmax[0], xmax=xminmax[1], color='lightgray', linewidth=0.5, zorder=1)
            # Set the xlim to the xminmax
            ax.set_xlim(xminmax)
            ax.set_xlabel('Log2 Fold-Change')
            # Set the yticks to the index values + 0.5
            ax.set_yticks(np.arange(len(df.index))+0.5)
            # Set the yticklabels to the index values
            ax.set_yticklabels(df.index)
            ax.set_ylabel(countgrp.capitalize())
            # Set ymin and ymax to -0.5 and len(df.index)+0.5
            ax.set_ylim(-0.5, len(df.index)+0.5)
            # Add light gray vertical lines at each integer
            ax.vlines(np.arange(round(xminmax[0]), round(xminmax[1])), ymin=-0.5, ymax=len(df.index)+0.5, color='lightgray',linewidth=0.5, zorder=-2)
            # Add a legend made manually from the bar colors
            handles = [plt.Rectangle((0,0),1,1, color=pal[i]) for i in cgrp1list] + \
                [plt.Rectangle((0,0),1,1, color='white')] + \
                [plt.Rectangle((0,0),1,1, color='black', linewidth=bardiff)] + \
                [plt.Rectangle((0,0),1,1, facecolor='white', linewidth=bardiff, edgecolor='black')]
            labels = cgrp1list.to_list() + ['','log2FC >= 1','log2FC < 1']
            ax.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.25, 1), frameon=False)
            # Add a title
            ax.set_title(f'{cgrp2} by {comparegrp1} {countgrp.capitalize()} Log2 Fold-Change')
            # Save the figure
            plt.savefig(f'{output}{comparegrp2}_{cgrp2}_by_{comparegrp1}_{countgrp}_log2fc.pdf', bbox_inches='tight')
            if threaded:
                threaded += f'Saving figure to {output}...\n'
                return threaded
            else:
                logger.info(f'Saving figure to {output}...')

if __name__ == '__main__':
    pass