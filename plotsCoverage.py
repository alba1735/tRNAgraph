#!/usr/bin/env python3

import numpy as np
import pandas as pd
import anndata as ad

from functools import partial
from multiprocessing import Pool

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns

class visualizer():
    '''
    Generate coverage plots for each sample in an AnnData object.
    '''
    def __init__(self, adata, threads, coverage_grp, coverage_obs, coverage_type, coverage_gap, coverage_method, colormap, output):
        self.threads = threads
        if coverage_grp not in adata.obs.columns:
            raise ValueError(f'Specified coveragegrp: {coverage_grp} not found in AnnData object.')
        self.coverage_obs = coverage_obs
        self.coverage_grp = coverage_grp
        self.coverage_type = coverage_type
        self.coverage_gap = coverage_gap
        self.coverage_method = coverage_method
        # Clean AnnData object for plotting
        self.adata, self.readstarts, self.readends = self.clean_adata(adata)
        if colormap != None: #and self.coverage_combine_all == False:
            self.coverage_pal = {k:v if v[0]!='#' else mplcolors.to_rgb(v) for k,v in colormap.items()}
        else:
            coverage_pal = sns.husl_palette(len(self.adata.obs[self.coverage_grp].unique()))
            self.coverage_pal = dict(zip(sorted(self.adata.obs[self.coverage_grp].unique()), coverage_pal))
        self.output = output

    def clean_adata(self, adata):
        '''
        Clean AnnData object for plotting.
        '''
        # Subset gaps from AnnData variables
        adata = adata[:,np.isin(adata.var.gap, self.coverage_gap)]
        # Drop nan values from AnnData
        adata = adata[~adata.obs[self.coverage_grp].isna()]
        # Subset just the readstarts and readends from AnnData variables
        readstarts = adata[:,np.isin(adata.var.coverage, ['readstarts'])].copy()
        readends = adata[:,np.isin(adata.var.coverage, ['readends'])].copy()
        # Subset by coverage type from AnnData variables
        adata = adata[:,np.isin(adata.var.coverage, [self.coverage_type])]

        return adata, readstarts, readends
    
    def __coverage_transform__(self, df, singlecol=False):
        '''
        Transform coverage data for plotting.
        '''
        # If the coverage method is mean, median, max, or min, transform the df by using groupby on column names
        if self.coverage_method == 'mean':
            if singlecol:
                df = df.mean(axis=1)
            else:
                df = df.T.groupby(level=0, observed=False).mean().T
        elif self.coverage_method == 'median':
            if singlecol:
                df = df.median(axis=1)
            else:
                df = df.T.groupby(level=0, observed=False).median().T
        elif self.coverage_method == 'max':
            if singlecol:
                df = df.max(axis=1)
            else:
                df = df.T.groupby(level=0, observed=False).max().T
        elif self.coverage_method == 'min':
            if singlecol:
                df = df.min(axis=1)
            else:
                df = df.T.groupby(level=0, observed=False).min().T
        elif self.coverage_method == 'sum':
            if singlecol:
                df = df.sum(axis=1)
            else:
                df = df.T.groupby(level=0, observed=False).sum().T
        else:
            raise ValueError(f'Invalid coverage method: {self.coverage_method}')
        
        return df

    def generate_combine(self):
        '''
        Generate combined coverage plots for all tRNAs using multiprocessing.
        '''
        # Generate list of tRNAs to plot sorting by name alphabetically and numerically with the copy number since sorting tRNAs is annoying
        ulist = sorted(self.adata.obs[self.coverage_obs].unique())
        # Sort by anticodon and then by copy number if tRNA
        if self.coverage_obs == 'trna':
            ulist = sorted(ulist, key=lambda x: ('-'.join(x.split('-')[:-1]), int(x.split('-')[-1])))
        # Generate list of tRNAs to plot split by 16 for each page
        ulist = [ulist[i*16:(i+1)*16] for i in range((len(ulist)+15)//16)]
        # Use multiprocessing to generate plotsand return them as a list so they can be saved to a pdf in order
        # Generate plots with confidence intervals
        outend = f'combined_{self.coverage_obs}_{self.coverage_type}_by_{self.coverage_grp}_with_ci_{self.coverage_method}.pdf'
        with PdfPages(f'{self.output}{outend}') as pdf:
            with Pool(self.threads) as p:
                pages = p.map(partial(self.generate_combine_page, coverage_fill='ci'), ulist)
            for page in pages:
                pdf.savefig(page, bbox_inches='tight')
        # Generate plots with fill
        outend = f'combined_{self.coverage_obs}_{self.coverage_type}_by_{self.coverage_grp}_with_fill_{self.coverage_method}.pdf'
        with PdfPages(f'{self.output}{outend}') as pdf:
            with Pool(self.threads) as p:
                pages = p.map(partial(self.generate_combine_page, coverage_fill='fill'), ulist)
            for page in pages:
                pdf.savefig(page, bbox_inches='tight')

    def generate_combine_page(self, ulist, coverage_fill):
        '''
        Generate combined coverage plots page for PdfPages via multiprocessing.
        '''
        # Generate figure
        fig_pdf = plt.figure(figsize=(24,22))
        for i, covobs in enumerate(ulist):
            # Turn off xaxis for all but bottom row
            xaxis = True if ulist.index(covobs) > 11 else False
            # Turn off legend for all but right column
            lgnd = True if ulist.index(covobs) % 4 == 3 else False
            # Generate subplot
            ax = fig_pdf.add_subplot(4,4,i+1)
            df = pd.DataFrame(self.adata[self.adata.obs[self.coverage_obs] == covobs].X.T, columns=self.adata[self.adata.obs[self.coverage_obs] == covobs].obs[self.coverage_grp].values)
            self.generate_plot(df, ax, covobs, coverage_fill=coverage_fill, lgnd=lgnd, xaxis=xaxis)
        return fig_pdf

    def generate_split(self):
        ulist = self.adata.obs[self.coverage_obs].unique()
        # Use multiprocessing to generate plots
        with Pool(self.threads) as p:
            p.map(self.generate_split_single, ulist)

    def generate_split_single(self, covobs):
        # Create df for single tRNA
        lgnd = True
        if self.coverage_grp == self.coverage_obs:
            lgnd = False
        df = pd.DataFrame(self.adata[self.adata.obs[self.coverage_obs] == covobs].X.T, columns=self.adata[self.adata.obs[self.coverage_obs] == covobs].obs[self.coverage_grp].values)
        low_coverage, outstart = self._assemble_plot_(df,covobs,lgnd,'ci')
        self._assemble_plot_(df,covobs,lgnd,'fill')
        # Generate plot with readstarts and readends
        rslist = self.adata.obs[self.coverage_grp].unique()
        if self.coverage_grp == self.coverage_obs:
            rslist = [covobs]
        for i in rslist:
            if not low_coverage:
                # Generate plot with readstarts and readends
                fig, ax = plt.subplots(figsize=(6,5.5))
                # Subset the df to the current group
                df_grp = df[i]
                if df_grp.sum().sum() != 0:
                    # Get the readstarts/ends for the current group and get the method of each position
                    df_readstarts = pd.DataFrame(self.readstarts[self.readstarts.obs[self.coverage_obs] == covobs].X.T, columns=self.readstarts[self.readstarts.obs[self.coverage_obs] == covobs].obs[self.coverage_grp].values)[i]
                    df_readends = pd.DataFrame(self.readends[self.readends.obs[self.coverage_obs] == covobs].X.T, columns=self.readends[self.readends.obs[self.coverage_obs] == covobs].obs[self.coverage_grp].values)[i]
                    df_readstarts = self.__coverage_transform__(df_readstarts, singlecol=True)
                    df_readends = self.__coverage_transform__(df_readends, singlecol=True)
                    # Combine the readstarts and readends into a single df
                    df_readstartends = pd.concat([df_readstarts, df_readends], axis=1)
                    df_readstartends.columns = ['readstarts', 'readends']
                    # Scale the df so that sum of each column is 1
                    df_readstartends = df_readstartends/df_readstartends.sum()
                    # Add position column to df
                    df_readstartends['position'] = df_readstartends.index
                    # Melt the df so that readstarts and readends are in the same column for plotting
                    df_readstartends = df_readstartends.melt(id_vars='position', value_vars=['readstarts', 'readends'])
                    self.generate_plot(df_grp, ax, covobs, coverage_fill='endstarts', rse=df_readstartends)
                    # Save the plot
                    if self.coverage_grp == self.coverage_obs:
                        outend = f'{covobs}_{self.coverage_type}_with_endstarts_{self.coverage_method}.pdf'
                    else:
                        outend = f'{covobs}_{self.coverage_type}_by_{self.coverage_grp}_{i}_with_endstarts_{self.coverage_method}.pdf'
                    plt.savefig(outstart+outend, bbox_inches='tight')
            plt.close()

    def _assemble_plot_(self, df, covobs, lgnd, cov_fill):
        # Generate plot with confidence intervals
        fig, ax = plt.subplots(figsize=(6,5.5))
        self.generate_plot(df, ax, covobs, coverage_fill=cov_fill, lgnd=lgnd)
        # Get max y value for plot
        outend = f'{covobs}_{self.coverage_type}_by_{self.coverage_grp}_with_{cov_fill}_{self.coverage_method}.pdf'
        if self.coverage_grp == self.coverage_obs:
            outend = f'{covobs}_{self.coverage_type}_with_{cov_fill}_{self.coverage_method}.pdf'
        if cov_fill == 'ci':
            outend = '_'.join(outend.split('_')[:-1]) + '.pdf'
        outstart = f'{self.output}{self.coverage_obs}/'
        low_coverage = False
        if ax.get_ylim()[1] < 20:
            outstart += 'low_coverage/'
            low_coverage = True
        plt.savefig(outstart+outend, bbox_inches='tight')
        plt.close()

        return low_coverage, outstart

    def generate_plot(self, df, ax, trna, coverage_fill, lgnd=True, xaxis=True, rse=None):
        '''
        Generate coverage plots for a single tRNA.
        '''
        # Get the cov method of all the columns
        if coverage_fill == 'ci': 
            sns.lineplot(data=df, linewidth=2, dashes=False, palette=self.coverage_pal, errorbar=('se',2), zorder=2, ax=ax)
        elif coverage_fill == 'fill':
            df = self.__coverage_transform__(df)
            sns.lineplot(data=df, linewidth=2, dashes=False, palette=self.coverage_pal, errorbar=('se',False), zorder=2, ax=ax)
        elif coverage_fill == 'endstarts':
            # Get the cov method of all the columns
            df = self.__coverage_transform__(df, singlecol=True)
            # Scale the df so that sum of all values is 1
            df = df/df.sum()
            sns.lineplot(data=df, linewidth=1.5, dashes=False, color='dimgrey', zorder=3, ax=ax)
            sns.histplot(data=rse, x='position', weights='value', hue='variable', palette=['magenta','cyan'], alpha=0.5, discrete=True, zorder=2, \
                         stat='probability', common_norm=False, linewidth=0, legend=True, ax=ax)
        # Fill the area under the curve with the mean by going in order of the column with the highest mean to the lowest if fill/both is specified
        if coverage_fill == 'fill':
            df_mean = df.T.groupby(level=0, observed=False).mean().T # This is the mean of the columns with the same name
            for i in df_mean.mean().sort_values(ascending=False).index:
                ax.fill_between(df_mean.index, df_mean[i], color=self.coverage_pal.get(i), alpha=0.35, zorder=1)
        # Set plot parameters
        plt.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=True, labelleft=True)
        ax.set_title(f'{trna} {self.coverage_type}')
        ax.set_ylabel("Normalized Readcounts")
        if coverage_fill == 'endstarts':
            ylim = [0,1]
            # Set as 0 to 100%
            ax.set_yticks([0,0.25,0.5,0.75,1])
            ax.set_yticklabels(['0%','25%','50%','75%','100%'])
        else:
            ylim = ax.get_ylim()
        ax.set_ylim(0, ylim[1])
        # Add dashed lines to plot for common tRNA modifications
        if coverage_fill == 'endstarts':
            # Create a df with the same columns as the rse df but empty
            names_df = pd.DataFrame({'position_start':[18.01,35.01,57.01], 'position_end':[18.01,35.01,57.01], 'name':['\nD-Arm','\nA-Arm','\nT-Arm']},\
                                     columns=['position_start','position_end','name'])
            # Add vertical dashed lines to plot for histobars that are above 10%
            for rt in rse['variable'].unique():
                current_start, current_bound = 0, 0
                endstart_switch = False
                for i, row in rse[rse['variable']==rt].iterrows():
                    if row['value'] >= 0.1:
                        if row['position'] - current_start > 1:
                            plt.plot([row['position']-0.5,row['position']-0.5],[0,ylim[1]],linewidth=1,ls='--',color='black', zorder=3)
                            current_bound = row['position']
                        current_start = row['position']
                        endstart_switch = True
                    else:
                        if endstart_switch:
                            plt.plot([row['position']-0.5,row['position']-0.5],[0,ylim[1]],linewidth=1,ls='--',color='black', zorder=3)
                            # Set name as the avg of the start and end positions
                            if current_bound+1 == row['position']:
                                names_df = pd.concat([names_df, pd.DataFrame({'position_start':[current_bound-1], 'position_end':[row['position']], 'name':[str(current_bound)]},\
                                                                              columns=['position_start','position_end','name'])])
                            else:
                                names_df = pd.concat([names_df, pd.DataFrame({'position_start':[current_bound-1], 'position_end':row['position'], 'name':[str(current_bound+1)+'-'+str(row['position'])]},\
                                                                            columns=['position_start','position_end','name'])])
                            endstart_switch = False
                if endstart_switch:
                    if current_bound == 75:
                        names_df = pd.concat([names_df, pd.DataFrame({'position_start':[current_bound], 'position_end':[75], 'name':[str(current_bound)]},\
                                                                      columns=['position_start','position_end','name'])])
                    else:
                        names_df = pd.concat([names_df, pd.DataFrame({'position_start':[current_bound], 'position_end':[75], 'name':[str(current_bound+1)+'-'+str(75)]},\
                                                                      columns=['position_start','position_end','name'])])
        else:
            for i in [37, 58]:
                plt.plot([i,i],[0,ylim[1]],linewidth=1,ls='--',color='black', zorder=3)
        # Set xaxis parameters
        if xaxis:
            ax.set_xlabel("Positions on tRNA")
        if coverage_fill == 'endstarts':
            ax.set_xlim(-0.5, 75.5)
            # Add horizontal dashed line to plot for histo above 10%
            plt.plot([-0.5,75.5],[0.1,0.1],linewidth=1,ls='--',color='black', zorder=3)
            # Add the xticks and xticklabels
            names_df['avg'] = (names_df['position_end'] + names_df['position_start'])/2
            names_df = names_df.sort_values(by=['avg'])
            ax.set_xticks(names_df['avg'])
            ax.set_xticklabels(names_df['name'])
        else:
            ax.set_xlim(0, 75) # ax.set_xlim(0, 73)
            ax.set_xticks([18.01,35.01,37,57.01,58])
            ax.set_xticklabels(['\nD-Arm','\nA-Arm','37','\nT-Arm','58'])
        # Add coverage regions to background
        for i in [[14,21],[32,38],[54,60],[10,25],[27,43],[49,65]]:
            ax.fill_between(i,[ylim[1],ylim[1]], color='#cacaca', alpha=0.35, zorder=0)
        # Capatilize the legend and move the legend outside the plot and remove the border around it
        if lgnd:
            if coverage_fill == 'endstarts':
                classes = ['Readstarts', 'Readends']
                class_colours = ['magenta', 'cyan']
                recs = []
                for i in range(0, len(class_colours)):
                    recs.append(mpatches.Rectangle((0, 0), 1, 1, fc=class_colours[i]))
                ax.legend(recs, classes, loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
                ax.legend_.set_title('Coverage')
            else:
                handles, labels = ax.get_legend_handles_labels()
                ax.legend(handles=handles, labels=[x.capitalize() for x in labels])
                ax.legend(loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
                ax.legend_.set_title(self.coverage_grp.capitalize())
        else:
            ax.legend_.remove()
        # Set the box aspect ratio to 1 so the plot is square
        plt.gca().set_box_aspect(1)

if __name__ == '__main__':
    pass