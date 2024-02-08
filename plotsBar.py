#!/usr/bin/env python3

import pandas as pd
import seaborn as sns

from functools import partial
from multiprocessing import Pool
import os

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


class visualizer():
    def __init__(self, adata, threads, barcol, bargrp, barsubgrp, barsort, barlabel, colormap, output):
        self.adata = adata
        self.threads = threads
        self.barcol = barcol
        self.bargrp = bargrp
        self.barsort = barsort
        self.barlabel = barlabel
        self.barsubgrp = barsubgrp
        self.colormap = colormap
        self.output = output
        
    def plots(self):
        '''
        Default function to plot the barplots.
        '''
        fig, axs = plt.subplots(figsize=(12,8))
        # Create a dictionary of colors for the bargrp
        if not self.colormap:
            self.colormap = dict(zip(sorted(self.adata.obs[self.bargrp].unique()), sns.color_palette('hls', len(self.adata.obs[self.bargrp].unique()))))
        # Create a df for the bargrps
        df = self.gen_df(self.adata)
        self.generate_plot(df, axs, self.bargrp)
        # Save figure
        plt.savefig(f'{self.output}{self.barcol}_percent_{self.bargrp}.pdf', bbox_inches='tight')

    def subplots(self):
        '''
        Generate subplots for each unique value in the barsubgrp column.
        '''
        ulist = self.adata.obs[self.barsubgrp].unique()
        # Generate list of tRNAs to plot split by 16 for each page
        ulist = [ulist[i*16:(i+1)*16] for i in range((len(ulist)+15)//16)]  # Split list into n sublists for pdfPages
        # Use multiprocessing to generate plotsand return them as a list so they can be saved to a pdf in order
        with PdfPages(f'{self.output}{self.barcol}_percent_{self.bargrp}_sub_{self.barsubgrp}.pdf') as pdf:
            with Pool(self.threads) as p:
                pages = p.map(partial(self.generate_combine_page), ulist)
            for page in pages:
                pdf.savefig(page, bbox_inches='tight')

    def generate_combine_page(self, ulist):
        '''
        Generate combined coverage plots page for PdfPages via multiprocessing.
        '''
        # Generate figure
        fig_pdf = plt.figure(figsize=(32,24))
        for i, bg in enumerate(ulist):
            # Turn off xaxis for all but the last 4 plots of ulist
            xaxis = False if i < len(ulist)-4 else True
            # Create a new df for each bargrp
            ax = fig_pdf.add_subplot(4,4,i+1)
            # df = self.gen_df(self.adata[self.adata.obs[self.barsubgrp].isin(bg)])
            df = self.gen_df(self.adata[self.adata.obs[self.barsubgrp].isin([bg])])
            self.colormap = dict(zip(sorted(df[self.bargrp].unique()), sns.color_palette('hls', len(df[self.bargrp].unique()))))
            self.generate_plot(df, ax, bg, xaxis=xaxis)
        return fig_pdf
    
    def gen_df(self, adata):
        df = adata.obs[[self.barcol,self.bargrp,'nreads_total_unique_norm']]
        # groupby the bargrp
        df = df.groupby([self.barcol,self.bargrp], observed=True).agg({'nreads_total_unique_norm':['mean','median','max','sum']})
        # remove the multiindex
        df = df.droplevel(0, axis=1)
        # reset the index
        df = df.reset_index()
        # Normalize the data to sum to 1 within each from the first column
        for i in df[self.barcol].unique():
            df.loc[df[self.barcol] == i, 'mean'] = df.loc[df[self.barcol] == i, 'mean'] / df.loc[df[self.barcol] == i, 'mean'].sum()
            df.loc[df[self.barcol] == i, 'median'] = df.loc[df[self.barcol] == i, 'median'] / df.loc[df[self.barcol] == i, 'median'].sum()
            df.loc[df[self.barcol] == i, 'max'] = df.loc[df[self.barcol] == i, 'max'] / df.loc[df[self.barcol] == i, 'max'].sum()
            df.loc[df[self.barcol] == i, 'sum'] = df.loc[df[self.barcol] == i, 'sum'] / df.loc[df[self.barcol] == i, 'sum'].sum()
            
        # If the combo of df[colgrp].unique() and df[bargrp].unique() does not exist in the df, add it with a value of 0
        for i in df[self.barcol].unique():
            for j in df[self.bargrp].unique():
                if len(df[(df[self.barcol] == i) & (df[self.bargrp] == j)]) == 0:
                    df = pd.concat([df, pd.DataFrame([[i,j,0,0,0,0]], columns=[self.barcol,self.bargrp,'mean','median','max','sum'])], ignore_index=True)
        # Reset the index
        df = df.reset_index(drop=True)
        # Enumerate the colgrp into a new column for positioning of the bars
        if self.barsort:
            # Get sort order start by dropping the NaN values
            sort_df = adata.obs[[self.barsort,self.barcol]].dropna().groupby(self.barcol)
            # Get the mean of the sort order
            sort_df = sort_df.mean()
            # Create a dictionary of the sort order
            sort_df = sort_df.sort_values(by=self.barsort)
            # Create a dictionary of the sort order
            col_dict = dict(zip(sort_df.index, range(len(sort_df.index))))
        else:
            # Create a dictionary of the sort order
            col_dict = dict(zip(sorted(df[self.barcol].unique()), range(len(df[self.barcol].unique()))))
        # Map the col_dict to the df
        df['x_pos'] = df[self.barcol].map(col_dict)
        # if barlabel, make a dict of that column from adata with the barcol as the key, then map it to the df
        if self.barlabel:
            label_df = adata.obs[[self.barcol,self.barlabel]].dropna()
            # Check that each value in the barcol has only one value in the barlabel
            for i in label_df[self.barcol].unique():
                if len(label_df[label_df[self.barcol] == i][self.barlabel].unique()) > 1:
                    raise ValueError(f'Each value in the {self.barcol} column must have only one value in the {self.barlabel} column.'+ \
                                    'Having more than one value in the --barlabel column for a single value in the --barcol column will cause the bar labels to be incorrect.')
            # Using the first value of the label for each barcol based on the previous check
            label_df = label_df.groupby(self.barcol).agg({self.barlabel:'first'})
            label_dict = dict(zip(label_df.index, label_df[self.barlabel]))
            df['label'] = df[self.barcol].map(label_dict)
        else:
            df['label'] = df[self.barcol]
        # Sort the df by the x_pos and bargrp
        df = df.sort_values(by=['x_pos',self.bargrp])

        return df

    def generate_plot(self, df, axs, title, xaxis=True):
        '''
        Create stacked barplots for readtypes.
        '''
        # Create a bar heights list
        bar_heights = [0] * len(df[self.barcol].unique())
        # Loop through the bargrp and create a bar for each
        legend_order = []
        for bar in sorted(df[self.bargrp].unique()):
            tdf = df[df[self.bargrp] == bar]
            tdf = tdf.sort_values(by=['x_pos'])
            tdf = tdf.reset_index(drop=True)
            axs.bar(tdf['x_pos'], tdf['mean'], 0.9, bottom=bar_heights, color=self.colormap[bar], linewidth=0.5, \
                    edgecolor='black', label=bar, clip_on=False)
            bar_heights += df[df[self.bargrp] == bar]['mean'].values
            legend_order.append(bar)
        # Remove all spines
        sns.despine(left=True, bottom=True)
        # Set the y axis limits
        axs.set_ylabel('Percentage of Total Reads')
        axs.set_yticks([0, 0.25, 0.50, 0.75, 1])
        axs.set_yticklabels(['0%', '25%', '50%', '75%', '100%'])
        axs.set_xlabel('Group')
        # rotate x labels and set them to the barcol
        if xaxis:
            plt.xticks(rotation=90)
            axs.set_xticks(range(len(df[self.barcol].unique())))
            axs.set_xticklabels(df['label'].unique())
        else:
            axs.set_xticks([])
        # Set the title to the bargrp
        axs.set_title(f'Percentage of Total Reads of {title.capitalize()}')
        # Legend - Set to the pal dictionary in reverse order
        handles, labels = axs.get_legend_handles_labels()
        axs.legend(handles[::-1], labels[::-1], title=self.bargrp.capitalize() + ' Group', \
                   loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
        # Set the box aspect ratio to 1 so the plot is square
        if self.barsubgrp:
            plt.gca().set_box_aspect(1)

if __name__ == '__main__':
    pass