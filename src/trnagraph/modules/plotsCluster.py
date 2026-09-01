#!/usr/bin/env python3

import logging

import seaborn as sns
import pandas as pd
import numpy as np

from . import toolsTG
from . import plotsPalette

import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

class visualizer():
    def __init__(self, adata, clustgrp, clustover, clusternumeric, clusterlabels, masking, colormap, output, threaded=True, read_basis=toolsTG.READ_BASIS_UNIQUE, settings=None):
        self.logger = logging.getLogger(__name__)
        self.adata = adata
        self.output = output
        self.threaded = threaded
        self.clustgrp = clustgrp
        self.overview = clustover
        self.numeric = clusternumeric
        self.clusterlabels = clusterlabels
        self.masking = masking
        self.colormap = colormap
        self.numericcolormap = plotsPalette.gradient(settings, 'ordered') # sns.diverging_palette(255, 85, s=255, l=70, sep=128, as_cmap=True)
        # --style's marker_size when set, otherwise the tuned default.
        self.settings = settings
        self.point_size = (settings or {}).get('marker_size') or 20
        # --allreads can only change how an EXISTING embedding is coloured: the UMAP
        # coordinates and HDBSCAN labels were computed by `analyze cluster` and written into
        # obs long before `graph` ran, so no graph-time flag can alter the projection itself.
        self.readtype = toolsTG.resolve_readtype('total', read_basis, adata)
        self.readtype_label = 'Unique Reads' if read_basis == toolsTG.READ_BASIS_UNIQUE else 'Reads'

    def generate_plots(self):
        # Generate overview plot
        if self.overview:
            self.colormap = None
            self.overviewPlot(self.adata, 'sample', self.output)
            self.clusterPlot(self.adata, 'sample', 'amino', self.output)
            self.clusterPlot(self.adata, 'sample', 'iso', self.output)
            self.clusterPlot(self.adata, 'sample', self.readtype, self.output, numeric=True)
            self.clusterPlot(self.adata, 'sample', 'fragment', self.output)
            self.clusterPlot(self.adata, 'sample', 'sample', self.output)
            self.clusterPlot(self.adata, 'sample', 'group', self.output)
            self.clusterPlot(self.adata, 'sample', 'sample_cluster', self.output)
            self.overviewPlot(self.adata, 'group', self.output)
            self.clusterPlot(self.adata, 'group', 'amino', self.output)
            self.clusterPlot(self.adata, 'group', 'iso', self.output)
            self.clusterPlot(self.adata, 'group', self.readtype, self.output, numeric=True)
            self.clusterPlot(self.adata, 'group', 'fragment', self.output)
            self.clusterPlot(self.adata, 'group', 'group', self.output)
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
        # Return the threaded string if threaded
        if self.threaded:
            return self.threaded

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
        fig, axs = plt.subplots(2, 3, figsize=(24,16))
        # Plot first through ninth subplots
        plot_list = [('Amino Acid','amino',0,0), ('Isotype','iso',0,1), (f'Total Number of {self.readtype_label}',self.readtype,0,2),
                     ('Fragment','fragment',1,0), (f'{umapgroup.capitalize()}',f'{umapgroup}',1,1), ('HDBScan',cluster,1,2)]
        for i in plot_list:
            if i[2] == 0 and i[3] == 2:
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
        fig.suptitle(f'UMAP Projection of tRNAs sorted by tRNAgraph {umapgroup}', y=0.925)
        # Save figure
        if self.threaded:
            self.threaded += f'Saving figure: {output}umap_{umapgroup}_overview.pdf\n'
        else:
            self.logger.info(f'Saving figure: {output}overview_{umapgroup}_overview.pdf')
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
            # discrete_colors, not sns.color_palette: numericcolormap is a resolved Colormap
            # object (so it can be handed straight to ScalarMappable below), and
            # color_palette takes a NAME. Sampling here also has to match color_palette's
            # own spacing, or these points shift colour against every other ordered plot.
            levels = sorted(pd.unique(adata.obs[clustgrp]))
            pal = dict(zip(levels, plotsPalette.discrete_colors(self.numericcolormap, len(levels))))
        else:
            pal = dict(zip(sorted(pd.unique(adata.obs[clustgrp][mask])), sns.color_palette("hls", len(pd.unique(adata.obs[clustgrp][mask])))))
        # Check if the user has defined a colormap
        if self.colormap:
            # Replace the palette values with the user defined ones if available
            pal = {k:self.colormap.get(k,v) for k,v in pal.items()}
        # Determine wether to mask NaN values after making pal and mask this helps with NaN problems
        if (adata.obs[clustgrp].isna().any() or clustgrp == cluster or self.masking) and not numeric:
            masking = True
            mask = adata.obs[cluster] >= 0
            mask = mask & ~adata.obs[clustgrp].isna()
        # Sort the adata object by the categorical variable for legend purposes
        adata = adata[adata.obs[clustgrp].sort_values().index, :]
        # Create the figure
        fig, axs = plt.subplots(figsize=toolsTG.figsize_for(self.settings, (8, 8)))
        if masking:
            sns.scatterplot(x=adata.obs[umap1][~mask], y=adata.obs[umap2][~mask], s=self.point_size, linewidth=0.25, ax=axs, color=np.array([(0.5,0.5,0.5)]), alpha=0.5)
            sns.scatterplot(x=adata.obs[umap1][mask], y=adata.obs[umap2][mask], s=self.point_size, linewidth=0.25, ax=axs, hue=adata.obs[clustgrp][mask], palette=pal)
        else:
            sns.scatterplot(x=adata.obs[umap1], y=adata.obs[umap2], s=self.point_size, linewidth=0.25, ax=axs, hue=adata.obs[clustgrp], palette=pal) #sns.diverging_palette(255, 85, s=255, l=70, sep=64, as_cmap=True))
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
            # sm = plt.cm.ScalarMappable(cmap=self.numericcolormap, norm=norm)
            sm = plt.cm.ScalarMappable(cmap=self.numericcolormap, norm=norm)
            sm.set_array([])
            # Remove the legend and add a colorbar
            if axs.get_legend():
                axs.get_legend().remove()
            fig.colorbar(sm, cax=fig.add_axes([1, 0.125, 0.03, 0.725]), label=clustgrp)
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
                    handles.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='darkgray', markersize=5))
                    labels.append('* Data Partially Masked')
                    handles.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='darkgray', markersize=5))
                    labels.append('** Data Fully Masked')
                axs.legend(handles, labels, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0., frameon=False, title=clustgrp)
        # Add title
        if numeric:
            plt.suptitle(f'UMAP projection of {clustgrp} by tRNAgraph {umapgroup}')
        else:
            plt.title(f'UMAP projection of {clustgrp} by tRNAgraph {umapgroup}')
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