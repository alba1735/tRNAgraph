#!/usr/bin/env python3

import numpy as np
import pandas as pd
import anndata as ad
import argparse

import directory_tools

import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns

class visualizer():
    '''
    Generate coverage plots for each sample in an AnnData object.
    '''
    def __init__(self, adata, coverage_grp, coverage_obs, coverage_type, coverage_gap, coverage_fill, output):
        if coverage_grp not in adata.obs.columns:
            raise ValueError('Specified coveragegrp not found in AnnData object.')
        self.coverage_grp = coverage_grp
        self.coverage_obs = coverage_obs
        self.coverage_type = coverage_type
        self.coverage_gap = coverage_gap
        self.coverage_fill = coverage_fill
        self.adata = self.clean_adata(adata)
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

    def generate_all(self):
        '''
        Generate coverage plots for all tRNAs.
        '''
        for trna in sorted(self.adata.obs.trna.unique()):
            self.generate_single(trna)

    def generate_single(self, trna):
        '''
        Generate coverage plots for a single tRNA.
        '''
        # Subset AnnData object to a single tRNA
        adata_t = self.adata[self.adata.obs.trna == trna].copy()
        df = pd.DataFrame(adata_t.X.T, columns=adata_t.obs[self.coverage_grp].values)

        # Generate plot
        plt.figure(figsize=(6,5.5))
        if self.coverage_fill == 'ci': 
            ax = sns.lineplot(data=df, linewidth=2, dashes=False, palette=self.coverage_pal, errorbar=('ci'), zorder=2)
        else:
            ax = sns.lineplot(data=df, linewidth=2, dashes=False, palette=self.coverage_pal, errorbar=('ci', False), zorder=2)

        # Fill the area under the curve with the mean by going in order of the column with the highest mean to the lowest if fill/both is specified
        if self.coverage_fill == 'fill':
            df_mean = df.groupby(df.columns, axis=1).mean() # This is the mean of the columns with the same name
            for i in df_mean.mean().sort_values(ascending=False).index:
                ax.fill_between(df_mean.index, df_mean[i], color=self.coverage_pal.get(i), alpha=0.35, zorder=1)

        # Set plot parameters
        ax.set_xlabel("Positions on tRNA")
        ax.set_ylabel("Normalized Readcounts")
        ylim = ax.get_ylim()
        ax.set_ylim(0, ylim[1])
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
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=handles, labels=[x.capitalize() for x in labels])
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
        ax.legend_.set_title(self.coverage_grp.capitalize())

         # Set the box aspect ratio to 1 so the plot is square
        plt.gca().set_box_aspect(1)
        # Save the plot
        if self.coverage_fill != 'none':
            plt.savefig('{}/{}_by_{}_{}_with_{}.pdf'.format(self.output, self.coverage_type, self.coverage_grp, trna, self.coverage_fill), bbox_inches='tight')
        else:
            plt.savefig('{}/{}_by_{}_{}.pdf'.format(self.output, self.coverage_type, self.coverage_grp, trna), bbox_inches='tight')
        plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='coverage_tools.py',
        description='Generate tRNA coverage plots'
    )

    parser.add_argument('-i', '--anndata', help='Specify AnnData input', required=True)
    parser.add_argument('-o', '--output', help='Specify output directory', default='coverage', required=False)
    parser.add_argument('--coveragegrp', help='Specify a grouping variable to generate coverage plots for (default: sample) (optional)', default='group')
    parser.add_argument('--coverageobs', help='Specify a observation subsetting for coverage plots (optional)', nargs='+', default=None)
    parser.add_argument('--coveragetype', help='Specify a coverage type for coverage plots corresponding to trax coverage file outputs (default: uniquecoverage) (optional)', \
        choices=['coverage', 'readstarts', 'readends', 'uniquecoverage', 'multitrnacoverage', 'multianticodoncoverage', 'multiaminocoverage','tRNAreadstotal', 'mismatchedbases', \
            'deletedbases', 'adenines', 'thymines', 'cytosines', 'guanines', 'deletions'], default='uniquecoverage')
    parser.add_argument('--coveragegap', help='Specify wether to include gaps in coverage plots (default: False) (optional)', default=False)
    parser.add_argument('--coveragefill', choices=['fill', 'ci', 'none'], help='Specify wether to fill area under coverage plots or use confidence intervals (default: ci) (optional)', default='ci')

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    directory_tools.builder(args.output)

    adata = ad.read_h5ad(args.anndata)

    visualizer(adata, args.coveragegrp, args.coverageobs, args.coveragetype, args.coveragegap, args.coveragefill, args.output).generate_all()