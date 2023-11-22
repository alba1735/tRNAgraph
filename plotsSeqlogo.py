
#!/usr/bin/env python3

import numpy as np
import pandas as pd
import anndata as ad
import argparse

import toolsDirectory

import logomaker

import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

class visualizer():
    def __init__(self, adata, output):
        self.adata = adata
        self.output = output
        self.psuedocount = 20
        self.coverage_grp = 'amino'
        self.pal_dict = {'T':'#ea4335','A':'#34a853','C':'#4285f4', 'G':'#fbbc05'}

        self.cutoff = 0
        self.cov_mask = False
        self.divergence_plot = False

    def generate_plots(self):
        for i in sorted(self.adata.obs[self.coverage_grp].values.unique()):
            seq_list = [i for i in self.adata.obs[self.adata.obs[self.coverage_grp] == i].refseq.values]
            # Create a matrix of base counts at each position
            df_readcounts, df_seqreads = self.norm_counts(seq_list, self.coverage_grp, i)
            # Plot the logo - information
            title_base = f'{self.coverage_grp}_{i}_normalized_bits'
            # if self.cov_mask:
            #     for x in range(len(df_readcounts)):
            #         if self.cov_mask[x] == False: 
            #             df_readcounts.iloc[x] = [0,0,0,0]
            #     title = title_base + f'_frequency_cutoff_{self.cutoff}_cov_testcovmask'
            # Convert the counts to information
            df_info_seqreads = logomaker.transform_matrix(df_seqreads, from_type='counts', to_type='information')
            df_info_readcounts = logomaker.transform_matrix(df_readcounts, from_type='counts', to_type='information')
            # df_counts = self.freq_matrix(df)
            self.logo_plot(df_info_seqreads, 'information', title=title_base+'_seqreads')
            self.logo_plot(df_info_readcounts, 'information', title=title_base+'_readcounts')
    
    def freq_matrix(self, df):
        fsize = df.max().max()
        for k,v in df.iterrows():
            tl = []
            for i in v:
                if i/fsize <= self.cutoff:
                    tl.append(0)
                else:
                    tl.append(i*2/fsize)
            df.at[k,['A','C','G','T']] = tl
        return df
    
    def norm_counts(self, seq_list, category, unit):
        # Create a matrix of 0s and add the pseudocount
        df_readcounts = pd.DataFrame(np.full((len(seq_list[0]),4), self.psuedocount), columns=['A','C','G','T'])
        df_seqreads = pd.DataFrame(np.full((len(seq_list[0]),4), 0), columns=['A','C','G','T'])
        # Subset the adata matrix to match the sequence list
        mx = self.adata[np.isin(self.adata.obs[category], [unit]),:]
        mx = mx[:,np.isin(mx.var.gap.values, [False])]
        # Generate a matrix of base counts at each position as a factor of the total counts
        base_list = ['adenines','cytosines','guanines','thymines']
        for i in range(len(base_list)):
            mxx = mx[:,np.isin(mx.var.coverage.values, [base_list[i]])]
            df_readcounts.iloc[:,i] += mxx.layers['raw'].sum(axis=0)
        # Iterate through the sequence list and add the counts to the matrix to generate a count matrix
        for seq in seq_list:
            for i in range(len(seq)):
                if seq[i] != '-':
                    df_seqreads.at[i,seq[i]] += 1
        # Remove the last 3 positions from the matrix to remove the CCA tail
        df_readcounts = df_readcounts.iloc[:-3,:]
        df_seqreads = df_seqreads.iloc[:-3,:]
        # Create a mask for the coverage
        if self.cov_mask:
            cov_df = pd.DataFrame(mx.X).mean(axis=0)
            self.cov_mask = [True if i >= self.cov_mask else False for i in cov_df/cov_df.max()]

        return df_readcounts, df_seqreads

    def logo_plot(self, df, seq_type, title):
        fig, ax = plt.subplots(figsize=(9,3))
        logomaker.Logo(df, color_scheme=self.pal_dict, ax=ax)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xlim([-1, 73])
        ax.set_xticks([])
        ax.set_title(title, fontsize='x-large')
        
        if self.divergence_plot == True:
            ax.spines['left'].set_bounds(-2, 2)
            y = -2.25
            yys = 0.935
            ax.set_ylim([-2.5, 2])
            ax.set_yticks([-2, -1, 0, 1, 2])
            ax.set_yticklabels(['-2','-1','0','1','2'])
            ax.set_ylabel("       Bits", fontsize=12, verticalalignment='center')
        else:
            ax.spines['left'].set_bounds(0, 2)
            y = -.16
            yys = 0.6
            ax.set_ylim([-0.25, 2])
            if seq_type == 'information':
                ax.set_yticks([0, 1, 2])
                ax.set_yticklabels(['0','1','2'])
                ax.set_ylabel("       Bits", fontsize=12, verticalalignment='center')
            elif seq_type == 'counts':
                ax.set_yticks([0,2])
                ax.set_yticklabels(['0%','100%'])
                ax.set_ylabel("       Frequency", fontsize=12, verticalalignment='center')    
                if self.cutoff > 0:
                    ax.plot([-1,73],[self.cutoff*2,self.cutoff*2],linewidth=1,ls='--',color='black')
        
        # Highlight the tRNA loops and arms
        ax.text(16.5, yys*y, 'D-Arm/Loop', fontsize=10, verticalalignment='top', horizontalalignment='center', color='white')
        ax.text(34, yys*y, 'A-Arm/Loop', fontsize=10, verticalalignment='top', horizontalalignment='center', color='white')
        ax.text(56, yys*y, 'T-Arm/Loop', fontsize=10, verticalalignment='top', horizontalalignment='center', color='white')
        ax.axvspan(8.5, 24.5, color='#f0f0f0', zorder=-2)
        ax.axvspan(25.5, 42.5, color='#f0f0f0', zorder=-2)
        ax.axvspan(47.5, 64.5, color='#f0f0f0', zorder=-2)
        ax.axvspan(12.5, 20.5, color='#cacaca', zorder=-1)
        ax.axvspan(30.5, 37.5, color='#cacaca', zorder=-1)
        ax.axvspan(52.5, 59.5, color='#cacaca', zorder=-1)
        # Plot dashed lines
        ax.plot([8.5, 8.5], [0, 2], color='k', linewidth=1, linestyle='--', zorder=0)
        ax.plot([24.5, 24.5], [0, 2], color='k', linewidth=1, linestyle='--', zorder=0)
        ax.plot([25.5, 25.5], [0, 2], color='k', linewidth=1, linestyle='--', zorder=0)
        ax.plot([42.5, 42.5], [0, 2], color='k', linewidth=1, linestyle='--', zorder=0)
        ax.plot([47.5, 47.5], [0, 2], color='k', linewidth=1, linestyle='--', zorder=0)
        ax.plot([64.5, 64.5], [0, 2], color='k', linewidth=1, linestyle='--', zorder=0)
        # Draw tRNA scheme at the bottom
        ax.axhline(y, color='k', linewidth=1)
        xs = np.arange(-3, len(df), 10)
        ys = y*np.ones(len(xs))
        ax.plot(xs, ys, marker='4', linewidth=0, markersize=7, color='k')
        ax.plot([8.5, 24.5], [y, y], color='k', linewidth=12, solid_capstyle='butt')
        ax.plot([25.5, 42.5], [y, y], color='k', linewidth=12, solid_capstyle='butt')
        ax.plot([47.5, 64.5], [y, y], color='k', linewidth=12, solid_capstyle='butt')

        # title_save = title.replace(' ','_').replace('(','').replace(')','').lower()
        plt.savefig('{}/{}.pdf'.format(self.output, title), bbox_inches='tight')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='logo_tools.py',
        description='Generate seqlogos from anndata.',
    )

    parser.add_argument('-i', '--anndata', help='Specify AnnData input', required=True)
    parser.add_argument('-o', '--output', help='Specify output directory', default='seqlogo', required=False)

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    toolsDirectory.builder(args.output)

    adata = ad.read_h5ad(args.anndata)

    visualizer(adata, args.output)