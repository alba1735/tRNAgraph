#!/usr/bin/env python3

import numpy as np
import pandas as pd
import anndata as ad
import argparse

from functools import partial
from itertools import repeat

from multiprocessing import Pool

import directory_tools

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns

class visualizer():
    '''
    Generate coverage plots for each sample in an AnnData object.
    '''
    def __init__(self, adata, threads, coverage_grp, coverage_obs, coverage_type, coverage_gap, colormap, output):
        self.threads = threads
        if coverage_grp not in adata.obs.columns:
            raise ValueError('Specified coveragegrp not found in AnnData object.')
        self.coverage_grp = coverage_grp
        self.coverage_obs = coverage_obs
        self.coverage_type = coverage_type
        self.coverage_gap = coverage_gap
        # Verify adata is valid for chosen coverage group or obs
        if adata.obs[self.coverage_grp].isna().any():
            raise ValueError('Coverage group contains NaN values.\nThis most likely means that you forgot to include samples in your metadata file that are present in your trax directory.\n' \
                             'Try adding the samples to your metadata file and rebuilding the AnnData object or selecting a different coverage group.')
        if self.coverage_obs:
            if adata.obs[self.coverage_obs].isna().any().any():
                raise ValueError('Coverage obs contains NaN values.\nThis most likely means that you forgot to include samples in your metadata file that are present in your trax directory.\n' \
                                'Try adding the samples to your metadata file and rebuilding the AnnData object or selecting a different coverage obs.')
        # Clean AnnData object for plotting
        self.adata = self.clean_adata(adata)
        if colormap != None:
            self.coverage_pal = {k:v if v[0]!='#' else mplcolors.to_rgb(v) for k,v in colormap.items()}
        else:
            coverage_pal = sns.husl_palette(len(self.adata.obs[self.coverage_grp].unique()))
            self.coverage_pal = dict(zip(sorted(self.adata.obs[self.coverage_grp].unique()), coverage_pal))
        self.output = output

    def clean_adata(self, adata):
        '''
        Clean AnnData object for plotting.
        '''
        # Subset AnnData observations if specified
        if self.coverage_obs:
            for k,v in self.coverage_obs.items():
                adata = adata[np.isin(adata.obs[k].values, v),:]
        # Subset by coverage type from AnnData variables
        adata = adata[:,np.isin(adata.var.coverage, [self.coverage_type])]
        # Subset gaps from AnnData variables
        adata = adata[:,np.isin(adata.var.gap, self.coverage_gap)]
        
        return adata

    def generate_combine(self):
        '''
        Generate combined coverage plots for all tRNAs using multiprocessing.
        '''
        # Generate list of tRNAs to plot sorting by name alphabetically and numerically with the copy number since sorting tRNAs is annoying
        trna_lists = self.adata.obs.trna.unique()
        trna_lists = sorted(trna_lists, key=lambda x: ('-'.join(x.split('-')[:-1]), int(x.split('-')[-1])))
        # Generate list of tRNAs to plot split by 16 for each page
        trna_lists = [trna_lists[i*16:(i+1)*16] for i in range((len(trna_lists)+15)//16)]  # Split list into n sublists for pdfPages
        # Use multiprocessing to generate plotsand return them as a list so they can be saved to a pdf in order
        # Generate plots with confidence intervals
        with PdfPages('{}/{}_by_{}_with_ci.pdf'.format(self.output, self.coverage_type, self.coverage_grp, self.coverage_obs)) as pdf:
            with Pool(self.threads) as p:
                pages = p.map(partial(self.generate_combine_page, coverage_fill='ci'), trna_lists)
            for page in pages:
                pdf.savefig(page, bbox_inches='tight')
        # Generate plots with fill
        with PdfPages('{}/{}_by_{}_with_fill.pdf'.format(self.output, self.coverage_type, self.coverage_grp, self.coverage_obs)) as pdf:
            with Pool(self.threads) as p:
                pages = p.map(partial(self.generate_combine_page, coverage_fill='fill'), trna_lists)
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
            df = pd.DataFrame(self.adata[self.adata.obs.trna == trna].X.T, columns=self.adata[self.adata.obs.trna == trna].obs[self.coverage_grp].values)
            self.generate_plot(df, ax, trna, coverage_fill=coverage_fill, lgnd=lgnd, xaxis=xaxis)
        return fig_pdf

    def generate_split(self):
        # Use multiprocessing to generate plots
        with Pool(self.threads) as p:
            p.map(self.generate_split_single, self.adata.obs.trna.unique())
    
    def generate_split_single(self, trna):
        # Creat df for single tRNA
        df = pd.DataFrame(self.adata[self.adata.obs.trna == trna].X.T, columns=self.adata[self.adata.obs.trna == trna].obs[self.coverage_grp].values)
        # Generate plot with confidence intervals
        fig, ax = plt.subplots(figsize=(6,5.5))
        self.generate_plot(df, ax, trna, coverage_fill='ci')
        # Get max y value for plot
        if ax.get_ylim()[1] >= 20:
            plt.savefig('{}/single/{}_{}_by_{}_with_ci.pdf'.format(self.output, trna, self.coverage_type, self.coverage_grp, self.coverage_obs), bbox_inches='tight')
            low_coverage = False
        else:
            plt.savefig('{}/single/low_coverage/{}_{}_by_{}_with_ci.pdf'.format(self.output, trna, self.coverage_type, self.coverage_grp, self.coverage_obs), bbox_inches='tight')
            low_coverage = True
        plt.close()
        # Generate plot with fill
        fig, ax = plt.subplots(figsize=(6,5.5))
        self.generate_plot(df, ax, trna, coverage_fill='fill')
        # Get max y value for plot
        if not low_coverage:
            plt.savefig('{}/single/{}_{}_by_{}_with_fill.pdf'.format(self.output, trna, self.coverage_type, self.coverage_grp, self.coverage_obs), bbox_inches='tight')
        else:
            plt.savefig('{}/single/low_coverage/{}_{}_by_{}_with_fill.pdf'.format(self.output, trna, self.coverage_type, self.coverage_grp, self.coverage_obs), bbox_inches='tight')
        plt.close()

    def generate_plot(self, df, ax, trna, coverage_fill, lgnd=True, xaxis=True):
        '''
        Generate coverage plots for a single tRNA.
        '''
        if coverage_fill == 'ci': 
            sns.lineplot(data=df, linewidth=2, dashes=False, palette=self.coverage_pal, errorbar=('ci'), zorder=2, ax=ax)
        else:
            sns.lineplot(data=df, linewidth=2, dashes=False, palette=self.coverage_pal, errorbar=('ci', False), zorder=2, ax=ax)

        # Fill the area under the curve with the mean by going in order of the column with the highest mean to the lowest if fill/both is specified
        if coverage_fill == 'fill':
            df_mean = df.groupby(df.columns, axis=1).mean() # This is the mean of the columns with the same name
            for i in df_mean.mean().sort_values(ascending=False).index:
                ax.fill_between(df_mean.index, df_mean[i], color=self.coverage_pal.get(i), alpha=0.35, zorder=1)

        # Set plot parameters
        ax.set_ylabel("Normalized Readcounts")
        ylim = ax.get_ylim()
        ax.set_ylim(0, ylim[1])

        if xaxis:
            ax.set_xlabel("Positions on tRNA")
        ax.set_xlim(0, 73)
        ax.set_xticks([18.01,35.01,37,57.01,58])
        ax.set_xticklabels(['\nD-Arm','\nA-Arm','37','\nT-Arm','58'])

        plt.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=True, labelleft=True)
        ax.set_title('{} {}'.format(trna, self.coverage_type))

        # Add dashed lines to plot for common tRNA modifications
        for i in [37, 58]:
            plt.plot([i,i],[0,ylim[1]],linewidth=1,ls='--',color='black', zorder=3)

        # Add coverage regions to background
        for i in [[14,21],[32,38],[54,60],[10,25],[27,43],[49,65]]:
            ax.fill_between(i,[ylim[1],ylim[1]], color='#cacaca', alpha=0.35, zorder=0)

        # Capatilize the legend and move the legend outside the plot and remove the border around it
        if lgnd:
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles=handles, labels=[x.capitalize() for x in labels])
            ax.legend(loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
            ax.legend_.set_title(self.coverage_grp.capitalize())
        else:
            ax.legend_.remove()

        # Set the box aspect ratio to 1 so the plot is square
        plt.gca().set_box_aspect(1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='coverage_tools.py',
        description='Generate tRNA coverage plots'
    )

    parser.add_argument('-i', '--anndata', help='Specify AnnData input', required=True)
    parser.add_argument('-o', '--output', help='Specify output directory', default='coverage', required=False)
    parser.add_argument('-n', '--threads', help='Specify number of threads to use (default: 1) (optional)', default=1, type=int)
    parser.add_argument('--coveragegrp', help='Specify a grouping variable to generate coverage plots for (default: group) (optional)', default='group')
    parser.add_argument('--coverageobs', help='Specify a observation subsetting for coverage plots (optional)', nargs='+', default=None)
    parser.add_argument('--coveragetype', help='Specify a coverage type for coverage plots corresponding to trax coverage file outputs (default: uniquecoverage) (optional)', \
        choices=['coverage', 'readstarts', 'readends', 'uniquecoverage', 'multitrnacoverage', 'multianticodoncoverage', 'multiaminocoverage','tRNAreadstotal', 'mismatchedbases', \
            'deletedbases', 'adenines', 'thymines', 'cytosines', 'guanines', 'deletions'], default='uniquecoverage')
    parser.add_argument('--coveragegap', help='Specify wether to include gaps in coverage plots (default: False) (optional)', default=False)
    parser.add_argument('--combineonly', help='Do not generate single tRNA coverage plot PDFs for every tRNA, only keep the combined output (optional)', action='store_false', required=False)
    parser.add_argument('--colormap', help='Specify a colormap for coverage plots (optional)', default=None)

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    directory_tools.builder(args.output)
    if args.splitpdfs:
        directory_tools.builder(args.output+'/single')
        directory_tools.builder(args.output+'/single/low_coverage')

    adata = ad.read_h5ad(args.anndata)

    if args.combineonly:
        visualizer(adata, args.threads, args.coveragegrp, args.coverageobs, args.coveragetype, args.coveragegap, args.colormap, args.output).generate_combined()
    else:
        visualizer(adata, args.threads, args.coveragegrp, args.coverageobs, args.coveragetype, args.coveragegap, args.colormap, args.output).generate_split()
        visualizer(adata, args.threads, args.coveragegrp, args.coverageobs, args.coveragetype, args.coveragegap, args.colormap, args.output).generate_combined()