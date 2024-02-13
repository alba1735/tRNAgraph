#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

class visualizer():
    def __init__(self, adata, clustgrp, clustover, clusternumeric, clusterlabels, masking, colormap, output):
        self.adata = adata
        self.output = output
        self.clustgrp = clustgrp
        self.overview = clustover
        self.numeric = clusternumeric
        self.clusterlabels = clusterlabels
        self.masking = masking
        self.colormap = colormap
        self.point_size = 20

    def generate_plots(self):
        # Generate overview plot
        if self.overview:
            self.colormap = None
            self.overviewPlot(self.adata, 'sample', self.output)
            self.clusterPlot(self.adata, 'sample', 'amino', self.output)
            self.clusterPlot(self.adata, 'sample', 'iso', self.output)
            self.clusterPlot(self.adata, 'sample', 'nreads_total_unique_norm', self.output, numeric=True)
            self.clusterPlot(self.adata, 'sample', 'sample_cluster', self.output)
            self.overviewPlot(self.adata, 'group', self.output)
            self.clusterPlot(self.adata, 'group', 'amino', self.output)
            self.clusterPlot(self.adata, 'group', 'iso', self.output)
            self.clusterPlot(self.adata, 'group', 'nreads_total_unique_norm', self.output, numeric=True)
            self.clusterPlot(self.adata, 'group', 'group_cluster', self.output)
        # Generate cluster plots if overview is not selected
        else:
            if self.clustgrp == 'sample_cluster':
                self.clusterPlot(self.adata, 'sample', self.clustgrp, self.output)
            elif self.clustgrp == 'group_cluster':
                self.clusterPlot(self.adata, 'group', self.clustgrp, self.output)
            else:
                self.clusterPlot(self.adata, 'sample', self.clustgrp, self.output, numeric=self.numeric)
                self.clusterPlot(self.adata, 'group', self.clustgrp, self.output, numeric=self.numeric)

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
        fig, axs = plt.subplots(2, 2, figsize=(16,16))
        # Plot first through ninth subplots
        plot_list = [('Amino Acid','amino',0,0), ('Isotype','iso',0,1), ('Total Number of Unique Reads','nreads_total_unique_norm',1,0), ('HDBScan',cluster,1,1)]
        for i in plot_list:
            if i[2] == 1 and i[3] == 0:
                pal = dict(zip(sorted(pd.unique(adata.obs[i[1]])), sns.color_palette("mako", len(pd.unique(adata.obs[i[1]])))))
            else:
                if i[1] == cluster:
                    # pal = dict(zip(sorted(pd.unique(adata.obs[i[1]])), sns.color_palette("hls", len(pd.unique(adata.obs[i[1]])-1))))
                    pal = dict(zip(sorted(pd.unique(adata.obs[i[1]][hdbscan_annotated])), sns.color_palette("hls", len(pd.unique(adata.obs[i[1]]))-1)))
                else:
                    pal = dict(zip(sorted(pd.unique(adata.obs[i[1]])), sns.color_palette("hls", len(pd.unique(adata.obs[i[1]])))))
            # Sort the adata object by the categorical variable for legend purposes
            adata = adata[adata.obs[i[1]].sort_values().index, :]
            # Plot the data
            if i[1] == cluster:
                sns.scatterplot(x=adata.obs[umap1][~hdbscan_annotated], y=adata.obs[umap2][~hdbscan_annotated], s=self.point_size, linewidth=0.25, ax=axs[i[2],i[3]], color=np.array([(0.5,0.5,0.5)]), alpha=0.5, legend=False)
                sns.scatterplot(x=adata.obs[umap1][hdbscan_annotated], y=adata.obs[umap2][hdbscan_annotated], s=self.point_size, linewidth=0.25, ax=axs[i[2],i[3]], hue=adata.obs[i[1]][hdbscan_annotated], palette=pal, legend=False)
            else:
                sns.scatterplot(x=adata.obs[umap1], y=adata.obs[umap2], s=self.point_size, ax=axs[i[2],i[3]], hue=adata.obs[i[1]], palette=pal, legend=False)
            axs[i[2],i[3]].set_title(i[0])
            axs[i[2],i[3]].set_xlabel('UMAP 1')
            axs[i[2],i[3]].set_ylabel('UMAP 2')
            axs[i[2],i[3]].set_xticks([])
            axs[i[2],i[3]].set_yticks([])
        # Add title
        fig.suptitle(f'UMAP Projection of tRNAs sorted by tRAX {umapgroup}', y=0.925)
        # Save figure
        print(f'Saving figure: {output}umap_{umapgroup}_overview.pdf')
        plt.savefig(f'{output}overview_{umapgroup}.pdf', bbox_inches='tight')
        plt.close()

    def clusterPlot(self, adata, umapgroup, clustgrp, output, numeric=False):
        # Define varibles
        umap1 = '_'.join([umapgroup,'umap1'])
        umap2 = '_'.join([umapgroup,'umap2'])
        cluster = '_'.join([umapgroup,'cluster'])
        # Subset the AnnData object to the umapgroup where not NaN (i.e. clustered data that wasn't filtered out)
        adata = adata[~adata.obs[cluster].isna(), :]
        # Set masking parameters
        masking = False
        mask = ~adata.obs[clustgrp].isna()
        # Create a palette for the categorical variable
        if numeric:
            pal = dict(zip(sorted(pd.unique(adata.obs[clustgrp])), sns.color_palette("mako_r", len(pd.unique(adata.obs[clustgrp])))))
        else:
            pal = dict(zip(sorted(pd.unique(adata.obs[clustgrp][mask])), sns.color_palette("hls", len(pd.unique(adata.obs[clustgrp][mask])))))
        # Check if the user has defined a colormap
        if self.colormap:
            # Replace the palette values with the user defined ones if available
            pal = {k:self.colormap.get(k,v) for k,v in pal.items()}
        # Determine wether to mask NaN values after making pal and mask this helps with NaN problems
        if adata.obs[clustgrp].isna().any() or clustgrp == cluster or self.masking:
            masking = True
            mask = adata.obs[cluster] >= 0
            mask = mask & ~adata.obs[clustgrp].isna()
        # Sort the adata object by the categorical variable for legend purposes
        adata = adata[adata.obs[clustgrp].sort_values().index, :]
        # Create the figure
        fig, axs = plt.subplots(figsize=(8,8))
        if masking:
            sns.scatterplot(x=adata.obs[umap1][~mask], y=adata.obs[umap2][~mask], s=self.point_size, linewidth=0.25, ax=axs, color=np.array([(0.5,0.5,0.5)]), alpha=0.5)
            sns.scatterplot(x=adata.obs[umap1][mask], y=adata.obs[umap2][mask], s=self.point_size, linewidth=0.25, ax=axs, hue=adata.obs[clustgrp][mask], palette=pal)
        else:
            sns.scatterplot(x=adata.obs[umap1], y=adata.obs[umap2], s=self.point_size, ax=axs, hue=adata.obs[clustgrp], palette=pal)
        # Add cluster labels
        if self.clusterlabels:
            for j in adata.obs[cluster][mask].unique():
                x = adata.obs[umap1][adata.obs[cluster] == j].mean()
                y = adata.obs[umap2][adata.obs[cluster] == j].mean()
                if umapgroup == 'group':
                    name = adata.obs[self.clusterlabels][adata.obs[cluster] == j].unique()[0]
                else:
                    name = adata.obs[cluster][adata.obs[cluster] == j].unique()[0]
                axs.text(x, y, name, fontsize=10, color='black', fontweight='bold')
        # Set title, x and y labels, and remove ticks
        axs.set_title(clustgrp)
        axs.set_xlabel('UMAP 1')
        axs.set_ylabel('UMAP 2')
        axs.set_xticks([])
        axs.set_yticks([])
        # Create legend from pal adding outside of plot and also reduce the size of the legend
        if numeric:
            norm = plt.Normalize(adata.obs[clustgrp].min(), adata.obs[clustgrp].max())
            sm = plt.cm.ScalarMappable(cmap="mako_r", norm=norm)
            sm.set_array([])
            # Remove the legend and add a colorbar
            axs.get_legend().remove()
            plt.colorbar(sm, ax=axs, pad=0.02)
        else:
            plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0., frameon=False, title=clustgrp)
            # Add an astrix to the legend to indicate the masked data in addition to the legend
            if masking:
                handles, labels = axs.get_legend_handles_labels()
                # if value counts is 0, then add a gray dot to the legend
                mask_dict = {}
                for k,v in adata.obs[clustgrp][mask].value_counts().items():
                    if v == 0:
                        mask_dict[k] = '**'
                    if v != adata.obs[clustgrp].value_counts()[k]:
                        mask_dict[k] = '*'
                # Change the color of the handles to gray for masked data
                if len(mask_dict.items()) > 0: # If there are any masked data, this prevents NaN from being added to the legend
                    for i in range(len(handles)):
                        labels[i] = labels[i] + mask_dict.get(labels[i],'')
                    # Add disclaimer to the legend
                    handles.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='w', markersize=5))
                    labels.append('')
                    handles.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='w', markersize=5))
                    labels.append('* Data Partially Masked')
                    handles.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='w', markersize=5))
                    labels.append('** Data Fully Masked')
                axs.legend(handles, labels, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0., frameon=False, title=clustgrp)
        # Add title
        plt.title(f'UMAP projection of {clustgrp} by tRAX {umapgroup}')
        # Set layout to equal with gca
        # plt.gca().set_aspect('equal', adjustable='box')
        # Save figure
        if masking:
            plt.savefig(output + f'{umapgroup}_by_{clustgrp}_masked.pdf', bbox_inches='tight')
        else:
            plt.savefig(output + f'{umapgroup}_by_{clustgrp}.pdf', bbox_inches='tight')
        plt.close()

if __name__ == '__main__':
    pass