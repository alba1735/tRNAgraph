#!/usr/bin/env python3

import pandas as pd
import anndata as ad
import argparse

import directory_tools

import sklearn.preprocessing
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns

class grapher():
    def __init__(self, adata, output, pcamarkers, pcacolors):
        self.adata = adata
        self.output = output
        self.pcamarkers = pcamarkers
        self.pcacolors = pcacolors

        if self.pcamarkers not in self.adata.obs.columns:
            raise ValueError('Specified pcamarkers not found in AnnData object.')
        if self.pcacolors not in self.adata.obs.columns:
            raise ValueError('Specified pcacolor not found in AnnData object.')

    def create(self):
        # Create a dataframe with trna, pcamarkers parameter, and nreads from adata and create dictory of sample and pcamarkers parameter for use in seaborn
        if self.pcamarkers == self.pcacolors:
            df = pd.DataFrame(self.adata.obs, columns=['trna', self.pcamarkers, 'nreads'])
            hue_dict = dict(zip(df[self.pcamarkers], df[self.pcamarkers]))
        else:
            df = pd.DataFrame(self.adata.obs, columns=['trna', self.pcamarkers, 'nreads', self.pcacolors])
            hue_dict = dict(zip(df[self.pcamarkers], df[self.pcacolors]))

        # Pivot the dataframe to have trna as the index, sample as the columns, and nreads as the values for dimensionality reduction
        df = df.pivot_table(index='trna', columns='sample', values='nreads')
        # Scale the data
        df = pd.DataFrame(sklearn.preprocessing.scale(df), columns=df.columns, index=df.index)

        # Create a PCA object
        pca = PCA(n_components=2)
        pca.fit_transform(df)
        evr = pca.explained_variance_ratio_
        print('Explained variance ratio: {}'.format(evr))

        # Transform the data and create a new dataframe
        df_pca = pd.DataFrame(pca.components_, columns=df.columns, index=['PC1', 'PC2']).T

        # Plot the data with seaborn
        plt.figure(figsize=(8, 8))
        ax = sns.scatterplot(data=df_pca, x='PC1', y='PC2', s=100, palette=sns.husl_palette(len(set(hue_dict.values()))), hue=hue_dict, legend='full')
        # Rename the x and y axis giving the percentage of variance explained by each PC
        ## It might be correct to swap the x and y axis here based on the output I have been getting
        ax.set_xlabel('PC1 ({:.2f}%)'.format(evr[0]*100))
        ax.set_ylabel('PC2 ({:.2f}%)'.format(evr[1]*100))
        # Capatilize the legend and move the legend outside the plot and remove the border around it
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=handles, labels=[x.capitalize() for x in labels])
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
        ax.legend_.set_title(self.pcacolors.capitalize())
        # Give the plot a title
        ax.set_title('PCA of {} colored by {}'.format(self.pcamarkers, self.pcacolors))
        # Remove the ticks and tick labels
        ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)

        # Set the box aspect ratio to 1 so the plot is square
        plt.gca().set_box_aspect(1)

        # Save the figure as a pdf without cropping the legend
        # plt.tight_layout()
        plt.savefig('{}/{}_by_{}_pca.pdf'.format(self.output, self.pcamarkers, self.pcacolors), bbox_inches='tight')
        print('PCA graph saved to {}/{}_by_{}_pca.pdf'.format(self.output, self.pcamarkers, self.pcacolors))    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='graphs_pca',
        description='Generate PCA graphs for each sample in an AnnData object.'
    )

    parser.add_argument('-i', '--anndata', help='Specify AnnData input', required=True)
    parser.add_argument('-o', '--output', help='Specify output directory', default='pca', required=False)
    parser.add_argument('--pcamarkers', help='Specify AnnData column to use for PCA markers (default: sample) (optional)', default='sample')
    parser.add_argument('--pcacolor', help='Specify AnnData column to color PCA markers by (default: sample) (optional)', default='sample')

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    directory_tools.builder(args.output)

    adata = ad.read_h5ad(args.anndata)

    grapher(adata, args.output, args.pcamarkers, args.pcacolor).create()