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
    def __init__(self, adata, threads, coverage_grp, coverage_combine, coverage_type, coverage_gap, colormap, output):
        self.threads = threads
        if coverage_grp not in adata.obs.columns:
            raise ValueError('Specified coveragegrp not found in AnnData object.')
        self.coverage_grp = coverage_grp
        self.coverage_combine = coverage_combine
        self.coverage_type = coverage_type
        self.coverage_gap = coverage_gap
        # Verify adata is valid for chosen coverage group or obs
        #if adata.obs[self.coverage_grp].isna().any():
        #    raise ValueError('Coverage group contains NaN values.\nThis most likely means that you forgot to include samples in your metadata file that are present in your trax directory.\n' \
        #                     'Try adding the samples to your metadata file and rebuilding the AnnData object or selecting a different coverage group.')
        # if self.coverage_obs:
        #     if adata.obs[self.coverage_obs].isna().any().any():
        #         raise ValueError('Coverage obs contains NaN values.\nThis most likely means that you forgot to include samples in your metadata file that are present in your trax directory.\n' \
        #                         'Try adding the samples to your metadata file and rebuilding the AnnData object or selecting a different coverage obs.')
        # Clean AnnData object for plotting
        self.adata, self.readstarts, self.readends = self.clean_adata(adata)
        if colormap != None and self.coverage_combine == None:
            self.coverage_pal = {k:v if v[0]!='#' else mplcolors.to_rgb(v) for k,v in colormap.items()}
        else:
            if self.coverage_combine:
                coverage_pal = sns.husl_palette(len(self.adata.obs[self.coverage_combine].unique()))
                self.coverage_pal = dict(zip(sorted(self.adata.obs[self.coverage_combine].unique()), coverage_pal))
            else:
                coverage_pal = sns.husl_palette(len(self.adata.obs[self.coverage_grp].unique()))
                self.coverage_pal = dict(zip(sorted(self.adata.obs[self.coverage_grp].unique()), coverage_pal))
        self.output = output

    def clean_adata(self, adata):
        '''
        Clean AnnData object for plotting.
        '''
        # Subset AnnData observations if specified
        # if self.coverage_obs:
        #     for k,v in self.coverage_obs.items():
        #         adata = adata[np.isin(adata.obs[k].values, v),:]
        # Subset gaps from AnnData variables
        adata = adata[:,np.isin(adata.var.gap, self.coverage_gap)]
        # If coverage combine is specified, take the mean of the coverage_combine obs column
        if self.coverage_combine:
            xdf = pd.DataFrame()
            obs = pd.DataFrame()
            for i in adata.obs[adata.obs[self.coverage_combine].notna()][self.coverage_combine].unique():
                # Take the mean of the coverage_combine obs column and add as row to obs
                obs = pd.concat([obs, adata[adata.obs[self.coverage_combine]==i,:].obs.iloc[0,:]], axis=1)
                # Take the mean of the coverage_combine X column
                xdf = pd.concat([xdf, pd.DataFrame(adata[adata.obs[self.coverage_combine]==i,:].X.mean(axis=0))], axis=1)
            # Transpose the dataframes into a new AnnData object
            xdf = xdf.T.reset_index(drop=True)
            obs = obs.T.reset_index(drop=True)
            adata = ad.AnnData(xdf, obs=obs, var=adata.var)
        # Subset just the readstarts and readends from AnnData variables
        readstarts = adata[:,np.isin(adata.var.coverage, ['readstarts'])].copy()
        readends = adata[:,np.isin(adata.var.coverage, ['readends'])].copy()
        # Subset by coverage type from AnnData variables
        adata = adata[:,np.isin(adata.var.coverage, [self.coverage_type])]

        return adata, readstarts, readends

    def generate_combine(self):
        '''
        Generate combined coverage plots for all tRNAs using multiprocessing.
        '''
        # Generate list of tRNAs to plot sorting by name alphabetically and numerically with the copy number since sorting tRNAs is annoying
        if self.coverage_combine:
            ulist = sorted(self.adata.obs[self.coverage_combine].unique())
        else:
            ulist = self.adata.obs.trna.unique()
            ulist = sorted(ulist, key=lambda x: ('-'.join(x.split('-')[:-1]), int(x.split('-')[-1])))
        # Generate list of tRNAs to plot split by 16 for each page
        ulist = [ulist[i*16:(i+1)*16] for i in range((len(ulist)+15)//16)]  # Split list into n sublists for pdfPages
        # Use multiprocessing to generate plotsand return them as a list so they can be saved to a pdf in order
        # Generate plots with confidence intervals
        outend = f'{self.coverage_type}_by_{self.coverage_grp}_with_ci.pdf'
        with PdfPages(f'{self.output}{outend}') as pdf:
            with Pool(self.threads) as p:
                pages = p.map(partial(self.generate_combine_page, coverage_fill='ci'), ulist)
            for page in pages:
                pdf.savefig(page, bbox_inches='tight')
        # Generate plots with fill
        outend = f'{self.coverage_type}_by_{self.coverage_grp}_with_fill.pdf'
        with PdfPages(f'{self.output}{outend}') as pdf:
            with Pool(self.threads) as p:
                pages = p.map(partial(self.generate_combine_page, coverage_fill='fill'), ulist)
            for page in pages:
                pdf.savefig(page, bbox_inches='tight')

    def generate_combine_page(self, trna_list, coverage_fill):
        '''
        Generate combined coverage plots page for PdfPages via multiprocessing.
        '''
        # Generate figure
        fig_pdf = plt.figure(figsize=(24,22))
        for i, trna in enumerate(trna_list):
            # Turn off xaxis for all but bottom row
            xaxis = True if trna_list.index(trna) > 11 else False
            # Turn off legend for all but right column
            lgnd = True if trna_list.index(trna) % 4 == 3 else False
            # Generate subplot
            ax = fig_pdf.add_subplot(4,4,i+1)
            # Create df for single tRNA
            if self.coverage_combine:
                df = pd.DataFrame(self.adata[self.adata.obs[self.coverage_combine] == trna].X.T, columns=[trna])
                lgnd = False
            else:
                df = pd.DataFrame(self.adata[self.adata.obs.trna == trna].X.T, columns=self.adata[self.adata.obs.trna == trna].obs[self.coverage_grp].values)
            self.generate_plot(df, ax, trna, coverage_fill=coverage_fill, lgnd=lgnd, xaxis=xaxis)
        return fig_pdf

    def generate_split(self):
        if self.coverage_combine:
            ulist = self.adata.obs[self.coverage_combine].unique()
        else:
            ulist = self.adata.obs.trna.unique()
        # Use multiprocessing to generate plots
        with Pool(self.threads) as p:
            p.map(self.generate_split_single, ulist)

    def generate_split_single(self, trna):
        # Create df for single tRNA
        lgnd = True
        if self.coverage_combine:
            df = pd.DataFrame(self.adata[self.adata.obs[self.coverage_combine] == trna].X.T, columns=[trna])
            lgnd = False
        else:
            df = pd.DataFrame(self.adata[self.adata.obs.trna == trna].X.T, columns=self.adata[self.adata.obs.trna == trna].obs[self.coverage_grp].values)
        # Generate plot with confidence intervals
        fig, ax = plt.subplots(figsize=(6,5.5))
        self.generate_plot(df, ax, trna, coverage_fill='ci', lgnd=lgnd)
        # Get max y value for plot
        outend = f'{trna}_{self.coverage_type}_by_{self.coverage_grp}_with_ci.pdf'
        if self.coverage_combine:
            outend = f'{self.coverage_combine}_{trna}_with_ci.pdf'
        if ax.get_ylim()[1] >= 20:
            outstart = f'{self.output}/single/'
            if self.coverage_combine:
                outstart += 'combined/'
            low_coverage = False
        else:
            outstart = f'{self.output}/single/low_coverage/'
            if self.coverage_combine:
                outstart = f'{self.output}/single/combined/low_coverage/'
            low_coverage = True
        plt.savefig(outstart+outend, bbox_inches='tight')
        plt.close()
        # Generate plot with fill
        fig, ax = plt.subplots(figsize=(6,5.5))
        self.generate_plot(df, ax, trna, coverage_fill='fill', lgnd=lgnd)
        # Get max y value for plot
        outend = f'{trna}_{self.coverage_type}_by_{self.coverage_grp}_with_fill.pdf'
        if self.coverage_combine:
            outend = f'{self.coverage_combine}_{trna}_with_fill.pdf'
        plt.savefig(outstart+outend, bbox_inches='tight')
        plt.close()
        # Generate plot with readstarts and readends
        rslist = self.adata.obs[self.coverage_grp].unique()
        if self.coverage_combine:
            rslist = [trna]
        for i in rslist:
            if not low_coverage:
                # Generate plot with readstarts and readends
                fig, ax = plt.subplots(figsize=(6,5.5))
                # Subset the df to the current group
                df_grp = df[i]
                if df_grp.sum().sum() != 0:
                    # Get the readstarts/ends for the current group and get the mean of each position
                    if self.coverage_combine:
                        df_readstarts = pd.DataFrame(self.readstarts[self.readstarts.obs[self.coverage_combine] == trna].X.T, columns=[trna])[i]
                        df_readends = pd.DataFrame(self.readends[self.readends.obs[self.coverage_combine] == trna].X.T, columns=[trna])[i]
                    else:
                        df_readstarts = pd.DataFrame(self.readstarts[self.readstarts.obs.trna == trna].X.T, columns=self.readstarts[self.readstarts.obs.trna == trna].obs[self.coverage_grp].values)[i]
                        df_readends = pd.DataFrame(self.readends[self.readends.obs.trna == trna].X.T, columns=self.readends[self.readends.obs.trna == trna].obs[self.coverage_grp].values)[i]
                        df_readstarts = df_readstarts.mean(axis=1)
                        df_readends = df_readends.mean(axis=1)
                    # Combine the readstarts and readends into a single df
                    df_readstartends = pd.concat([df_readstarts, df_readends], axis=1)
                    df_readstartends.columns = ['readstarts', 'readends']
                    # Scale the df so that sum of each column is 1
                    df_readstartends = df_readstartends/df_readstartends.sum()
                    # Add position column to df
                    df_readstartends['position'] = df_readstartends.index
                    # Melt the df so that readstarts and readends are in the same column for plotting
                    df_readstartends = df_readstartends.melt(id_vars='position', value_vars=['readstarts', 'readends'])
                    self.generate_plot(df_grp, ax, trna, coverage_fill='endstarts', rse=df_readstartends)
                    # Save the plot
                    outend = f'{trna}_{self.coverage_type}_by_{self.coverage_grp}_{i}_with_endstarts.pdf'
                    if self.coverage_combine:
                        outend = f'{self.coverage_combine}_{trna}_with_endstarts.pdf'
                    plt.savefig(outstart+outend, bbox_inches='tight')
            plt.close()

    def generate_plot(self, df, ax, trna, coverage_fill, lgnd=True, xaxis=True, rse=None):
        '''
        Generate coverage plots for a single tRNA.
        '''
        if coverage_fill == 'ci': 
            sns.lineplot(data=df, linewidth=2, dashes=False, palette=self.coverage_pal, errorbar=('ci'), zorder=2, ax=ax)
        elif coverage_fill == 'fill':
            sns.lineplot(data=df, linewidth=2, dashes=False, palette=self.coverage_pal, errorbar=('ci', False), zorder=2, ax=ax)
        elif coverage_fill == 'endstarts':
            # Get the mean of all the columns
            if not self.coverage_combine:
                df = df.mean(axis=1)
            # Scale the df so that sum of all values is 1
            df = df/df.sum()
            sns.lineplot(data=df, linewidth=1.5, dashes=False, color='dimgrey', zorder=3, ax=ax)
            sns.histplot(data=rse, x='position', weights='value', hue='variable', palette=['magenta','cyan'], alpha=0.5, discrete=True, zorder=2, \
                         stat='probability', common_norm=False, linewidth=0, legend=True, ax=ax)
        # Fill the area under the curve with the mean by going in order of the column with the highest mean to the lowest if fill/both is specified
        if coverage_fill == 'fill':
            # df_mean = df.groupby(df.columns, axis=1).mean() 
            # The above is deprecated because it doesn't work with the new version of pandas use DataFrame.groupby with axis=1 is deprecated. Do `frame.T.groupby(...)` without axis instead.
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