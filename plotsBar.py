#!/usr/bin/env python3

import numpy as np
import pandas as pd
import seaborn as sns
# import anndata as ad
# import argparse

# import toolsTG

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


class visualizer():
    def __init__(self, adata, barcol, bargrp, colormap, output):
        self.barcol = barcol
        self.bargrp = bargrp
        df = self.gen_df(adata)
        # Create a dictionary of colors for the bargrp
        if colormap:
            self.palette = colormap
        else:
            self.palette = dict(zip(sorted(df[bargrp].unique()), sns.color_palette("hls", len(pd.unique(df[bargrp])))))
        self.output = output
        self.stacked_barplots(df)
    
    def gen_df(self, adata):
        df = adata.obs[[self.barcol,self.bargrp,'nreads_total_unique_norm']]
        # groupby the bargrp
        df = df.groupby([self.barcol,self.bargrp], observed=True).agg({'nreads_total_unique_norm':['mean','sum','count']})
        # remove the multiindex
        df = df.droplevel(0, axis=1)
        # reset the index
        df = df.reset_index()
        # Normalize the data to sum to 1 within each from the first column
        for i in df[self.barcol].unique():
            df.loc[df[self.barcol] == i, 'mean'] = df.loc[df[self.barcol] == i, 'mean'] / df.loc[df[self.barcol] == i, 'mean'].sum()
            df.loc[df[self.barcol] == i, 'sum'] = df.loc[df[self.barcol] == i, 'sum'] / df.loc[df[self.barcol] == i, 'sum'].sum()
            df.loc[df[self.barcol] == i, 'count'] = df.loc[df[self.barcol] == i, 'count'] / df.loc[df[self.barcol] == i, 'count'].sum()
        # If the combo of df[colgrp].unique() and df[bargrp].unique() does not exist in the df, add it with a value of 0
        for i in df[self.barcol].unique():
            for j in df[self.bargrp].unique():
                if len(df[(df[self.barcol] == i) & (df[self.bargrp] == j)]) == 0:
                    df = pd.concat([df, pd.DataFrame([[i,j,0,0,0]], columns=[self.barcol,self.bargrp,'mean','sum','count'])], ignore_index=True)
        # Sort the df by the colgrp and bargrp
        df = df.sort_values(by=[self.barcol,self.bargrp])
        # Reset the index
        df = df.reset_index(drop=True)
        # Enumerate the colgrp into a new column for positioning of the bars
        col_dict = dict(zip(sorted(df[self.barcol].unique()), range(len(df[self.barcol].unique()))))
        df['x_pos'] = df[self.barcol].map(col_dict)

        return df

    def stacked_barplots(self, df):
        '''
        Create stacked barplots for readtypes.
        '''
        # Create stacked barplots
        fig, axs = plt.subplots(figsize=(12,8))
        # Create a bar heights list
        bar_heights = [0] * len(df[self.barcol].unique())
        # Loop through the bargrp and create a bar for each
        legend_order = []
        for bar in sorted(df[self.bargrp].unique()):
            tdf = df[df[self.bargrp] == bar]
            tdf = tdf.sort_values(by=['x_pos'])
            tdf = tdf.reset_index(drop=True)
            axs.bar(tdf['x_pos'], tdf['count'], 0.9, bottom=bar_heights, color=self.palette[bar], linewidth=0.5, \
                    edgecolor='black', label=bar, clip_on=False)
            bar_heights += df[df[self.bargrp] == bar]['count'].values
            legend_order.append(bar)
        # Remove all spines
        sns.despine(left=True, bottom=True)
        # Set the y axis limits
        axs.set_ylabel('Percentage of Total Reads')
        axs.set_yticks([0, 0.25, 0.50, 0.75, 1])
        axs.set_yticklabels(['0%', '25%', '50%', '75%', '100%'])
        axs.set_title('Percentage of Total Reads')
        axs.set_xlabel('Group')
        # rotate x labels and set them to the barcol
        plt.xticks(rotation=90)
        axs.set_xticks(range(len(df[self.barcol].unique())))
        axs.set_xticklabels(df[self.barcol].unique())
        # Legend - Set to the pal dictionary in reverse order
        handles, labels = axs.get_legend_handles_labels()
        axs.legend(handles[::-1], labels[::-1], title=self.bargrp.capitalize() + ' Group', \
                   loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
        # Set the box aspect ratio to 1 so the plot is square
        # plt.gca().set_box_aspect(1)
        # Save figure
        plt.savefig(f'{self.output}{self.barcol}_percent_{self.bargrp}_stacked.pdf', bbox_inches='tight')
        plt.close()

if __name__ == '__main__':
    pass