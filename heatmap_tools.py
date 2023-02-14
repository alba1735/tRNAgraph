#!/usr/bin/env python3

import seaborn as sns
import numpy as np
import anndata as ad
import argparse

import directory_tools
import analysis_tools

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


def visualizer(adata, grp, readtype, cutoff, heatbound, output):
    '''
    Generate heatmap visualizations for each group in an AnnData object.
    '''
    # Create a color palette for the heatmap
    cmap = sns.diverging_palette(255, 85, s=255, l=70, sep=20, as_cmap=True)
    # Create a correlation matrix from reads stored in adata observations
    df = analysis_tools.log2fc_df(adata, grp, readtype, cutoff)
    # Create a pdf with a heatmap for sorted by each group on each page
    with PdfPages(f'{output}/{grp}_{readtype}_{cutoff}_{heatbound}_sum_heatmap.pdf') as pdf:
        for col in [i for i in df.columns.tolist() if 'log2' in i]:
            print(col)
            # Sort df by the sum of the log2FC values for each row
            tdf = df.sort_values(by=col, ascending=False)
            # subset df to only include the top 25 and bottom 25 rows
            tdf = tdf.iloc[:heatbound, :].append(tdf.iloc[-heatbound:, :])

            # Create a heatmap
            fig, axs = plt.subplots(1, len(tdf.columns)+1, figsize=(6, 12))
            midpoint = len(tdf.columns) // 2
            for i, pair in enumerate(tdf.columns):
                # Create a heatmap for each pair
                if i < midpoint:
                    if i == 0:
                        sns.heatmap(tdf[pair].to_frame(), ax=axs[i], cmap=cmap, center=0, vmax=2, vmin=-2, cbar=False)
                    else:
                        sns.heatmap(tdf[pair].to_frame(), ax=axs[i], cmap=cmap, center=0, vmax=2, vmin=-2, cbar=False, yticklabels=False)
                    axs[i].set_ylabel('')
                    axs[i].tick_params(axis='x', labelrotation=90)
                else:
                    sns.heatmap(-np.log10(tdf[pair].to_frame()), ax=axs[i+1], cmap='Greens', vmin=0, vmax=3, cbar=False, yticklabels=False)
                    axs[i+1].set_ylabel('')
                    axs[i+1].tick_params(axis='x', labelrotation=90)

            # Add a colorbar to the left of the heatmap for log2FC
            axs[midpoint].set_axis_off()
            cbar_ax = fig.add_axes([0.465, 0.15, 0.025, 0.7])
            fig.colorbar(axs[0].collections[0], cax=cbar_ax)

            # Add a colorbar to the right of the heatmap for pvals
            cbar_ax = fig.add_axes([0.915, 0.15, 0.025, 0.7])
            fig.colorbar(axs[len(tdf.columns)].collections[0], cax=cbar_ax)
            # Change the tick labels to be -log10(pval) instead of pval
            cbar_ax.set_yticks([0, 1.3, 3])
            cbar_ax.set_yticklabels(['1', '0.05', '<=0.001'])

            # Add a title to the figure
            fig.suptitle(f'Heatmap of log2FC of tRNA read counts between groups\nsorted by {col}', fontsize=12, y=0.925)

            # Save figure
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
    parser.add_argument('--heatrt', help='Specify readtype to use for heatmap', default='nreads_total_unique_norm', required=False)
    parser.add_argument('--heatcutoff', help='Specify readcount cutoff to use for heatmap', default=80, required=False)
    parser.add_argument('--heatbound', help='Specify range to use for bounding the heatmap to top and bottom counts', default=25, required=False)

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    directory_tools.builder(args.output)

    adata = ad.read_h5ad(args.anndata)

    visualizer(adata, args.heatgrp, args.heatrt, args.heatcutoff, args.heatbound, args.output)
