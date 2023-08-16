#!/usr/bin/env python3

import seaborn as sns
import numpy as np
import pandas as pd
import anndata as ad
import argparse

import directory_tools
import analysis_tools

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


def visualizer(adata, grp, readtypes, cutoff, heatbound, heatsubplots, output):
    '''
    Generate heatmap visualizations for each group in an AnnData object.
    '''
    if grp not in adata.obs.columns:
        raise ValueError('Specified group not found in AnnData object.')

    for readtype in readtypes:
        readtype = f'nreads_{readtype}_norm'
        # Create a color palette for the heatmap
        cmap = sns.diverging_palette(255, 85, s=255, l=70, sep=20, as_cmap=True)
        # Create a correlation matrix from reads stored in adata observations
        df = analysis_tools.log2fc_df(adata, grp, readtype, cutoff)
        # save df to csv
        df.to_csv(f'{output}/{grp}_{readtype}_{cutoff}_{heatbound}_sum_heatmap.csv')
        # Create a pdf with a heatmap for sorted by each group on each page
        with PdfPages(f'{output}/{grp}_{readtype}_{cutoff}_{heatbound}_sum_heatmap.pdf') as pdf:
            for col in [i for i in df.columns.tolist() if 'log2' in i]:
                # Sort df by the sum of the log2FC values for each row
                tdf = df.sort_values(by=col, ascending=False)
                # subset df to only include the top and bottom heatbound values only if the heatmap is larger than heatbound
                if len(tdf) > heatbound:
                    tdf = pd.concat([tdf.iloc[:heatbound, :], tdf.iloc[-heatbound:, :]])
                # Create a heatmap
                if len(tdf.columns)>=20:
                    fig, axs = plt.subplots(1, 2, figsize=(16,12))
                else:
                    fig, axs = plt.subplots(1, 2, figsize=(6,12))

                # fig, axs = plt.subplots(1, len(tdf.columns)+3, figsize=(int((len(tdf.columns)+1)/2), heatbound/2))
                # midpoint = len(tdf.columns) // 2
                # for i, pair in enumerate(tdf.columns):
                #     # Create a heatmap for each pair
                #     if i < midpoint:
                #         if i == 0:
                #             sns.heatmap(tdf[pair].to_frame(), ax=axs[i], cmap=cmap, center=0, vmax=2, vmin=-2, cbar=False, square=True)
                #         else:
                #             sns.heatmap(tdf[pair].to_frame(), ax=axs[i], cmap=cmap, center=0, vmax=2, vmin=-2, cbar=False, yticklabels=False, square=True)
                #         axs[i].set_ylabel('')
                #         axs[i].tick_params(axis='x', labelrotation=90)
                #     else:
                #         sns.heatmap(-np.log10(tdf[pair].to_frame()), ax=axs[i+3], cmap='Greens', vmin=0, vmax=3, cbar=False, yticklabels=False, square=True)
                #         axs[i+3].set_ylabel('')
                #         axs[i+3].tick_params(axis='x', labelrotation=90)

                log_tdf = tdf[[i for i in tdf.columns if 'log2' in i]]
                pval_tdf = tdf[[i for i in tdf.columns if 'pval' in i]]
                sns.heatmap(log_tdf, ax=axs[0], cmap=cmap, center=0, vmax=2, vmin=-2, cbar=True, square=True, cbar_kws={'fraction':0.05, 'pad':0.05})
                sns.heatmap(-np.log10(pval_tdf), ax=axs[1], cmap='Greens', vmin=0, vmax=3, cbar=True, yticklabels=False, square=True, cbar_kws={'fraction':0.05, 'pad':0.05})
                axs[0].tick_params(axis='x', labelrotation=90)
                axs[1].tick_params(axis='x', labelrotation=90)
                axs[1].set_ylabel('')

                cbar = axs[0].collections[0].colorbar
                cbar.set_ticks([-2, -1, 0, 1, 2])
                cbar.set_ticklabels(['<=-2', '-1', '0', '1', '>=2'])
                cbar.ax.tick_params(labelsize=8) 
                cbar = axs[1].collections[0].colorbar
                cbar.set_ticks([0, 1.3, 3])
                cbar.set_ticklabels(['1', '0.05', '<=0.001'])
                cbar.ax.tick_params(labelsize=8)

                # Add a colorbar to the left of the heatmap for log2FC
                # axs[midpoint].set_axis_off()
                # axs[midpoint+1].set_axis_off()
                # axs[midpoint+2].set_axis_off()
                # cbar_ax = fig.add_axes([0.465, 0.15, 0.025, 0.7])
                # fig.colorbar(axs[0].collections[0], cax=cbar_ax)
                # Add a colorbar to the right of the heatmap for pvals
                # cbar_ax = fig.add_axes([0.915, 0.15, 0.025, 0.7])
                # fig.colorbar(axs[len(tdf.columns)].collections[0], cax=cbar_ax)
                # Change the tick labels to be -log10(pval) instead of pval
                # cbar_ax.set_yticks([0, 1.3, 3])
                # cbar_ax.set_yticklabels(['1', '0.05', '<=0.001'])

                # Add a title to the figure
                plt.suptitle(f'Heatmap of log2FC of tRNA read counts between groups \nsorted by {col}', fontsize=12, y=1.15)
                axs[0].set_title('log2FC', fontsize=10)
                axs[1].set_title('pval', fontsize=10)

                # Save figure
                print(f'Saving heatmap for {readtype} {col}...')
                if heatsubplots:
                    plt.savefig(f'{output}/{grp}_{readtype}_{cutoff}_{heatbound}_{col}_heatmap.pdf', bbox_inches='tight')
                pdf.savefig(bbox_inches='tight')
                plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='pca_tools.py',
        description='Generate heatmap visualizations for each group in an AnnData object.'
    )

    parser.add_argument('-i', '--anndata', help='Specify AnnData input', required=True)
    parser.add_argument('-o', '--output', help='Specify output directory', default='heatmap', required=False)
    parser.add_argument('--heatgrp', help='Specify group to use for heatmap', default='group', required=False)
    parser.add_argument('--heatrts', choices=['whole_unique', 'fiveprime_unique', 'threeprime_unique', 'other_unique', 'total_unique', \
         'wholecounts', 'fiveprime', 'threeprime', 'other', 'total', 'antisense', 'wholeprecounts', 'partialprecounts', 'trailercounts', 'all'], \
            help='Specify readtypes to use for heatmap (default: whole_unique, fiveprime_unique, threeprime_unique, other_unique, total_unique) (optional)', \
                nargs='+', default=['whole_unique', 'fiveprime_unique', 'threeprime_unique', 'other_unique', 'total_unique'], required=False)
    parser.add_argument('--heatcutoff', help='Specify readcount cutoff to use for heatmap', default=80, required=False, type=int)
    parser.add_argument('--heatbound', help='Specify range to use for bounding the heatmap to top and bottom counts', default=25, required=False)
    parser.add_argument('--heatsubplots', help='Specify wether to generate subplots for each comparasion in addition to the sum (default: False)', action='store_true', required=False)

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    directory_tools.builder(args.output)

    adata = ad.read_h5ad(args.anndata)

    visualizer(adata, args.heatgrp, args.heatrts, args.heatcutoff, args.heatbound, args.heatsubplots, args.output)
