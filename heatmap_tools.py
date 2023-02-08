#!/usr/bin/env python3

import seaborn as sns
import anndata as ad
import argparse

import directory_tools
import analysis_tools

import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


def visualizer(adata, grp, readtype, cutoff, output):
    '''
    Generate heatmap visualizations for each group in an AnnData object.
    '''
    # Create a correlation matrix from reads stored in adata observations
    df = analysis_tools.log2fc_df(adata, grp, readtype, cutoff) 
    # subset df to only include the top 25 and bottom 25 rows
    df = df.iloc[:25, :].append(df.iloc[-25:, :])

    # Create a color palette for the heatmap
    cmap = sns.diverging_palette(255, 85, s=255, l=70, sep=20, as_cmap=True)
    # Create a heatmap
    fig, axs = plt.subplots(1, len(df.columns), figsize=(6, 12))
    for i, pair in enumerate(df.columns):
        # Create a heatmap for each pair
        if i != 0:
            sns.heatmap(df[pair].to_frame(), ax=axs[i], cmap=cmap, center=0, vmax=2, vmin=-2, cbar=False, yticklabels=False)
        else:
            sns.heatmap(df[pair].to_frame(), ax=axs[i], cmap=cmap, center=0, vmax=2, vmin=-2, cbar=False)
        axs[i].set_ylabel('')
        axs[i].tick_params(axis='x', labelrotation=90)

    # Add a colorbar to the right of the heatmap
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(axs[0].collections[0], cax=cbar_ax)
    cbar.set_label('log2FC', rotation=270, labelpad=15)

    # Add a title to the figure
    fig.suptitle('Heatmap of log2FC of tRNA read counts between groups', fontsize=16, y=0.925)

    # Save figure
    plt.savefig(f'{output}/heatmap.pdf', bbox_inches='tight')
    plt.close()

    # print(df[df.index == 'tRNA-Arg-TCT-4'])
    # print(df[df.index == 'tRNA-Arg-TCT-4'])

    # print((np.log2(100)-np.log2(10)))
    # print((np.log2(10)-np.log2(10)))
    # print((np.log2(10)-np.log2(100)))

    # print(2**(np.log2(100)-np.log2(10)))
    # print(2**(np.log2(10)-np.log2(10)))
    # print(2**(np.log2(10)-np.log2(100)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='pca_tools.py',
        description='Generate heatmap visualizations for each group in an AnnData object.'
    )

    parser.add_argument('-i', '--anndata', help='Specify AnnData input', required=True)
    parser.add_argument('-o', '--output', help='Specify output directory', default='heatmap', required=False)
    parser.add_argument('--heatgrp', help='Specify group to use for heatmap', default='group', required=False)
    parser.add_argument('--heatrt', help='Specify readtype to use for heatmap', default='nreads_total_norm', required=False)
    parser.add_argument('--heatcutoff', help='Specify readcount cutoff to use for heatmap', default=80, required=False)

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    directory_tools.builder(args.output)

    adata = ad.read_h5ad(args.anndata)

    visualizer(adata, args.heatgrp, args.heatrt, args.heatcutoff, args.output)
