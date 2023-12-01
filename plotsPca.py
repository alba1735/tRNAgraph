#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
# import anndata as ad
# import argparse

# import toolsTG

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


def visualizer(adata, pcamarkers, pcacolors, pcareadtypes, colormap, output):
    '''
    Generate PCA visualizations for each sample in an AnnData object.
    '''
    if pcamarkers not in adata.obs.columns:
        raise ValueError('Specified pcamarkers not found in AnnData object.')
    if pcacolors not in adata.obs.columns:
        raise ValueError('Specified pcacolor not found in AnnData object.')
    # Create a list of readtypes to iterate over if 'all' is specified
    if 'all' in pcareadtypes:
        pcareadtypes = ['whole_unique', 'fiveprime_unique', 'threeprime_unique', 'other_unique', 'total_unique', 'wholecounts',
                        'fiveprime', 'threeprime', 'other', 'total', 'antisense', 'wholeprecounts', 'partialprecounts', 'trailercounts']

    for readtype in pcareadtypes:
        # Rename the readtype column to nreads_{readtype}_norm to match adata.obs
        rt = f'nreads_{readtype}_norm'
        # Create a dataframe with trna, pcamarkers parameter, and nreads from adata and create dictory of sample and pcamarkers parameter for use in seaborn
        if pcamarkers == pcacolors:
            df = pd.DataFrame(adata.obs, columns=['trna', pcamarkers, rt])
            hue_dict = dict(zip(df[pcamarkers], df[pcamarkers]))
        else:
            df = pd.DataFrame(adata.obs, columns=['trna', pcamarkers, rt, pcacolors])
            hue_dict = dict(zip(df[pcamarkers], df[pcacolors]))
        if colormap != None:
            colormap = {k:v if v[0]!='#' else mplcolors.to_rgb(v) for k,v in colormap.items()}
            for v in hue_dict.values():
                if v not in colormap:
                    print(f'Color {v} not found in colormap. Using default colors instead.')
                    colormap = None
                    break

        # Pivot the dataframe to have trna as the index, sample as the columns, and nreads as the values for dimensionality reduction
        df = df.pivot_table(index='trna', columns='sample', values=rt)
        # Scale the data
        df = pd.DataFrame(StandardScaler().fit_transform(df), columns=df.columns, index=df.index)

        # Create a PCA object
        pca = PCA(n_components=min(len(df.columns), 5))
        pca.fit_transform(df)
        evr = pca.explained_variance_ratio_
        print('Principal components: {}'.format([f'PC{x}' for x in range(1, len(evr)+1)]))
        print('Explained variance: {}'.format([f'{i:.4f}' for i in pca.explained_variance_]))
        print('Explained variance ratio: {}'.format([f'{i*100:.2f}%' for i in evr]))
        # Transform the data and create a new dataframe
        pca_index = ['PC{}'.format(x) for x in range(1, len(evr)+1)]
        df_pca = pd.DataFrame(pca.components_, columns=df.columns, index=pca_index).T

        # Plot the explained variance ratio
        plt.figure(figsize=(6, 6))
        ax = sns.barplot(x=['PC{}'.format(x) for x in range(1, len(evr)+1)], y=evr, palette=sns.husl_palette(len(evr)), hue=['PC{}'.format(x) for x in range(1, len(evr)+1)])
        # Set the x and y labels and title
        ax.set_xlabel('Principal Component')
        ax.set_ylabel('Explained Variance Ratio')
        ax.set_title('Explained Variance Ratio of Principal Components')
        # Set the box aspect ratio to 1 so the plot is square
        plt.gca().set_box_aspect(1)
        # Save the plot
        plt.savefig(f'{output}{pcamarkers}_by_{pcacolors}_{readtype}_evr.pdf', bbox_inches='tight')
        print(f'Explained variance ratio graph saved to {output}{pcamarkers}_by_{pcacolors}_{readtype}_evr.pdf')
        plt.close()

        # Plot the data with seaborn
        plt.figure(figsize=(8, 8))
        if colormap:
            ax = sns.scatterplot(data=df_pca, x='PC1', y='PC2', s=100, palette=colormap, hue=hue_dict, legend='full')
        else:
            ax = sns.scatterplot(data=df_pca, x='PC1', y='PC2', s=100, palette=sns.husl_palette(len(set(hue_dict.values()))), hue=hue_dict, legend='full')
        ax.set_xlabel('PC1 ({:.2f}%)'.format(evr[0]*100))
        ax.set_ylabel('PC2 ({:.2f}%)'.format(evr[1]*100))
        # Capatilize the legend and move the legend outside the plot and remove the border around it
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=handles, labels=[x.capitalize() for x in labels])
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
        ax.legend_.set_title(pcacolors.capitalize())
        # Give the plot a title
        ax.set_title(f'PCA of {pcamarkers} colored by {pcacolors}')
        # Remove the ticks and tick labels
        ax.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelbottom=False, labelleft=False)
        # Set the box aspect ratio to 1 so the plot is square
        plt.gca().set_box_aspect(1)
        # Save the plot
        plt.savefig(f'{output}{pcamarkers}_by_{pcacolors}_{readtype}_pca.pdf', bbox_inches='tight')
        print(f'PCA graph saved to {output}{pcamarkers}_by_{pcacolors}_{readtype}_pca.pdf')
        plt.close()

        # Plot pairplot of the data with seaborn
        plt.figure(figsize=(10, 10))
        # Rename columns to add explained variance ratio to each PC as well as add colors for use in seaborn
        df_pca.columns = ['PC{} ({:.2f}%)'.format(i, evr[i-1]*100) for i in range(1, len(df_pca.columns)+1)]
        df_pca[pcacolors.capitalize()] = [hue_dict[x] for x in df_pca.index]
        # Plot the pairplot
        if colormap:
            ax = sns.pairplot(df_pca, hue=pcacolors.capitalize(), palette=colormap, hue_order=sorted(set(hue_dict.values())))
        else:
            ax = sns.pairplot(df_pca, hue=pcacolors.capitalize(), palette=sns.husl_palette(len(set(hue_dict.values()))), hue_order=sorted(set(hue_dict.values())))
        # Remove the ticks and tick labels
        ax.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelbottom=False, labelleft=False)
        # Give the plot a title
        ax.fig.suptitle(f'PCA Pairplot of {pcamarkers} colored by {pcacolors}', y=1.02)
        # Set the box aspect ratio to 1 so the plot is square
        plt.gca().set_box_aspect(1)
        plt.savefig(f'{output}{pcamarkers}_by_{pcacolors}_{readtype}_pairplot.pdf', bbox_inches='tight')
        print(f'Pairplot graph saved to {output}{pcamarkers}_by_{pcacolors}_{readtype}_pairplot.pdf')
        plt.close()


if __name__ == '__main__':
    pass
    # parser = argparse.ArgumentParser(
    #     prog='pca_tools.py',
    #     description='Generate PCA visualizations for each sample in an AnnData object.'
    # )

    # parser.add_argument('-i', '--anndata',
    #                     help='Specify AnnData input', required=True)
    # parser.add_argument('-o', '--output', 
    #                     help='Specify output directory', default='pca', required=False)
    # parser.add_argument('--pcamarkers', 
    #                     help='Specify AnnData column to use for PCA markers (default: sample) (optional)', default='sample')
    # parser.add_argument('--pcacolor', 
    #                     help='Specify AnnData column to color PCA markers by (default: group) (optional)', default='group')
    # parser.add_argument('--pcareadtypes', 
    #                     choices=['whole_unique', 'fiveprime_unique', 'threeprime_unique', 'other_unique', 'total_unique', 
    #                     'wholecounts', 'fiveprime', 'threeprime', 'other', 'total', 'antisense', 'wholeprecounts', 'partialprecounts', 'trailercounts', 'all'],
    #                     help='Specify read types to use for PCA markers (default: total_unique, total) (optional)', nargs='+', default=['total_unique', 'total'])
    # parser.add_argument('--colormap', help='Specify a colormap for coverage plots (optional)', default=None)

    # args = parser.parse_args()

    # # Create output directory if it doesn't exist
    # toolsTG.builder(args.output)

    # adata = ad.read_h5ad(args.anndata)

    # visualizer(adata, args.pcamarkers, args.pcacolor, args.pcareadtypes, args.colormap, args.output)
