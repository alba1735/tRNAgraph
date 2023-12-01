#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
# import anndata as ad
# import argparse

# import toolsTG

import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

class visualizer():
    def __init__(self, adata, output):
        self.adata = adata
        self.output = output

    def main(self):
        # Generate overview plot
        self.overviewPlot(self.adata, 'sample', self.output)
        self.overviewPlot(self.adata, 'group', self.output)

    def overviewPlot(self, adata, umapgroup, output):
        # Define varibles
        umap1 = '_'.join([umapgroup,'umap1'])
        umap2 = '_'.join([umapgroup,'umap2'])
        cluster = '_'.join([umapgroup,'cluster'])
        # Subset the AnnData object to the umapgroup where not NaN (i.e. clustered data that wasn't filtered out)
        adata = adata[~adata.obs[cluster].isna(), :]
        # Create a list of clusters greater than or equal to 0 in size to filter out non-clustered reads
        hdbscan_annotated = adata.obs[cluster] >= 0
        # Create a 3 x 3 subplot with the umap projection and the cluster labels as the last subplot
        fig, axs = plt.subplots(2, 2, figsize=(12,12))
        # Plot first through ninth subplots
        plot_list = [('Amino Acid','amino',0,0), ('Isotype','iso',0,1), ('Total Number of Unique Reads','nreads_total_unique_norm',1,0), ('HDBScan',cluster,1,1)]
        for i in plot_list:
            # Sort the adata object by the categorical variable for legend purposes
            adata = adata[adata.obs[i[1]].sort_values().index, :]
            # Create a color palette for the categorical variable
            if i[2] == 1 and i[3] == 0:
                pal = sns.color_palette("mako")
            else:
                pal_size = len(pd.unique(adata.obs[i[1]]))
                if i[0] == 'HDBScan': # Remove from length if HDBScan because -1 will be plotted independently in grey
                    pal_size -= 1
                pal = sns.color_palette("hls", pal_size)

            # Plot the data
            if i[1] == cluster:
                sns.scatterplot(x=adata.obs[umap1][~hdbscan_annotated], y=adata.obs[umap2][~hdbscan_annotated], s=8, ax=axs[i[2],i[3]], c=(0.5,0.5,0.5), alpha=0.5, legend=False)
                sns.scatterplot(x=adata.obs[umap1][hdbscan_annotated], y=adata.obs[umap2][hdbscan_annotated], s=8, ax=axs[i[2],i[3]], hue=adata.obs[i[1]][hdbscan_annotated], palette=pal, legend=False)
            else:
                sns.scatterplot(x=adata.obs[umap1], y=adata.obs[umap2], s=8, ax=axs[i[2],i[3]], hue=adata.obs[i[1]], palette=pal, legend=False)
            axs[i[2],i[3]].set_title(i[0])
            axs[i[2],i[3]].set_xlabel('UMAP 1')
            axs[i[2],i[3]].set_ylabel('UMAP 2')
            axs[i[2],i[3]].set_xticks([])
            axs[i[2],i[3]].set_yticks([])
        # Add title
        fig.suptitle(f'UMAP Projection of tRNAs sorted by tRAX {umapgroup}')
        # Save figure
        print(f'Saving figure: {output}umap_{umapgroup}_overview.pdf')
        plt.savefig(f'{output}umap_{umapgroup}_overview.pdf', bbox_inches='tight')
        plt.close()

if __name__ == '__main__':
    pass
    # parser = argparse.ArgumentParser(
    #     prog='plotsCluster.py',
    #     description='Generate plots from UMAP projection of tRNA dataset',
    # )

    # parser.add_argument('-i', '--anndata',
    #                     help='Specify AnnData input', required=True)
    # parser.add_argument('-o', '--output', 
    #                     help='Specify output directory', default='umap', required=False)

    # args = parser.parse_args()

    # # Create output directory if it doesn't exist
    # toolsTG.builder(args.output)

    # adata = ad.read_h5ad(args.anndata)

    # visualizer(adata, args.output)