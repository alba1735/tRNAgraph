#!/usr/bin/env python3

import numpy as np
import pandas as pd
# import anndata as ad
# import argparse

# import toolsTG

import logomaker
# from sklearn.preprocessing import RobustScaler
# from sklearn.preprocessing import MinMaxScaler
# from sklearn.preprocessing import StandardScaler
# from sklearn.preprocessing import Normalizer

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

class visualizer():
    def __init__(self, adata, output):
        self.adata = adata
        self.output = output
        self.psuedocount = 0
        self.coverage_grp = 'amino'
        self.pal_dict = {'T':'#ea4335','A':'#34a853','C':'#4285f4', 'G':'#fbbc05'}

        # self.cutoff = 0 
        # self.cov_mask = False
        # self.divergence_plot = False

    def generate_plots(self):
        for i in sorted(self.adata.obs[self.coverage_grp].values.unique()):
            seq_list = [i for i in self.adata.obs[self.adata.obs[self.coverage_grp] == i].refseq.values]
            df_seqinfo = self.build_seqdata(seq_list)
            self.consensus_plot(df_seqinfo, seq_list, i)

            # Create a matrix of base counts at each position
            # self.logo_counts(seq_list, 'deseq', self.coverage_grp, i)

            # self.logo_counts(seq_list, 'robust', self.coverage_grp, i)
            # self.logo_counts(seq_list, 'normalizer', self.coverage_grp, i)
            # self.logo_counts(seq_list, 'minmax', self.coverage_grp, i)
            # self.logo_counts(seq_list, 'standard', self.coverage_grp, i)

            # self.logo_counts(seq_list, 'raw', self.coverage_grp, i)

    def build_seqdata(self, seq_list):
        # Create a matrix of 0s for the sequence reads
        df_seqreads = pd.DataFrame(np.full((len(seq_list[0]),4), 0), columns=['A','C','G','T'])
        # Iterate through the sequence list and add the counts to the matrix to generate a count matrix
        for seq in seq_list:
            for i in range(len(seq)):
                if seq[i] != '-':
                    df_seqreads.at[i,seq[i]] += 1
        # Remove the last 3 positions from the matrix to remove the CCA tail
        df_seqreads = df_seqreads.iloc[:-3,:]
        # Convert the counts to information
        df_seqinfo = logomaker.transform_matrix(df_seqreads, from_type='counts', to_type='information')

        return df_seqinfo
    
    def logo_counts(self, seq_list, norm_type, category, unit):
        # Create a matrix of 0s and add the pseudocount
        df_background = pd.DataFrame(np.full((len(seq_list[0]),4), self.psuedocount), columns=['A','C','G','T'])
        df_sumcounts = pd.DataFrame(np.full((len(seq_list[0]),4), self.psuedocount), columns=['A','C','G','T'])
        df_meancounts = pd.DataFrame(np.full((len(seq_list[0]),4), self.psuedocount), columns=['A','C','G','T'])
        # df_mediancounts = pd.DataFrame(np.full((len(seq_list[0]),4), self.psuedocount), columns=['A','C','G','T'])
        df_actualreads = pd.DataFrame(np.full((len(seq_list[0]),4), self.psuedocount), columns=['A','C','G','T'])
        df_seqreads = pd.DataFrame(np.full((len(seq_list[0]),4), 0), columns=['A','C','G','T'])
        # Subset the adata matrix to match the sequence list
        ad = self.adata[:,np.isin(self.adata.var.gap.values, [False])]
        mx = ad[np.isin(ad.obs[category], [unit]),:]

        # Apply the scikit-learn robust scaler to normalize the counts of mx
        mx.layers['deseq'] = mx.X
        # elif norm_type == 'robust':
        #     mx.layers['robust'] = RobustScaler().fit_transform(mx.X)
        # elif norm_type == 'minmax':
        #     mx.layers['minmax'] = MinMaxScaler().fit_transform(mx.X)
        # elif norm_type == 'standard':
        #     mx.layers['standard'] = StandardScaler().fit_transform(mx.X)
        # elif norm_type == 'normalizer':
        #     mx.layers['normalizer'] = Normalizer().fit_transform(mx.X)

        print(unit, norm_type, mx.layers[norm_type].min(), mx.layers[norm_type].max(), mx.layers[norm_type].mean(), '\n')
        
        # Generate a matrix of base counts at each position as a factor of the total counts
        base_list = ['adenines','cytosines','guanines','thymines']
        for i in range(len(base_list)):
            adx = ad[:,np.isin(ad.var.coverage.values, [base_list[i]])]
            mxx = mx[:,np.isin(mx.var.coverage.values, [base_list[i]])]



            # df_background.iloc[:,i] += adx.layers[norm_type].sum(axis=0)
            df_sumcounts.iloc[:,i] += mxx.layers[norm_type].sum(axis=0)
            df_meancounts.iloc[:,i] += mxx.layers[norm_type].mean(axis=0)
            # df_mediancounts.iloc[:,i] += np.median(mxx.layers[norm_type], axis=0)
            # replace all values over 0 with 1
            mxx.layers[norm_type][mxx.layers[norm_type] > 0] = 1
            df_actualreads.iloc[:,i] += mxx.layers[norm_type].sum(axis=0)

        # Iterate through the sequence list and add the counts to the matrix to generate a count matrix
        for seq in seq_list:
            for i in range(len(seq)):
                if seq[i] != '-':
                    df_seqreads.at[i,seq[i]] += 1

        # Remove the last 3 positions from the matrix to remove the CCA tail
        # df_background = df_background.iloc[:-3,:]
        df_sumcounts = df_sumcounts.iloc[:-3,:]
        df_meancounts = df_meancounts.iloc[:-3,:]
        # df_mediancounts = df_mediancounts.iloc[:-3,:]
        df_actualreads = df_actualreads.iloc[:-3,:]
        df_seqreads = df_seqreads.iloc[:-3,:]

        # print(df_background)
        # print(df_readcounts)
        # print(df_seqreads)
        # print(df_seqreads_alt)
        # exit()

        # Create a mask for the coverage
        # if self.cov_mask:
        #     cov_df = pd.DataFrame(mx.X).mean(axis=0)
        #     self.cov_mask = [True if i >= self.cov_mask else False for i in cov_df/cov_df.max()]

        title_base = f'{self.coverage_grp}_{unit}'

        # if self.cov_mask:
        #     for x in range(len(df_readcounts)):
        #         if self.cov_mask[x] == False: 
        #             df_readcounts.iloc[x] = [0,0,0,0]
        #     title = title_base + f'_frequency_cutoff_{self.cutoff}_cov_testcovmask'

        # Convert the counts to information
        df_info_actualreads = logomaker.transform_matrix(df_actualreads, from_type='counts', to_type='information')
        df_info_seqreads = logomaker.transform_matrix(df_seqreads, from_type='counts', to_type='information')
        df_info_sumcounts = logomaker.transform_matrix(df_sumcounts, from_type='counts', to_type='information')
        df_info_meancounts = logomaker.transform_matrix(df_meancounts, from_type='counts', to_type='information')
        # df_info_mediancounts = logomaker.transform_matrix(df_mediancounts, from_type='counts', to_type='information')
        # df_counts = self.freq_matrix(df)

        self.logo_plot(df_info_actualreads, 'information', title=title_base+f'_information_actualreads_{norm_type}_pseudo{self.psuedocount}')
        self.logo_plot(df_info_seqreads, 'information', title=title_base+'_information_seqreads')
        self.logo_plot(df_info_sumcounts, 'information', title=title_base+f'_information_sumcounts_{norm_type}_pseudo{self.psuedocount}')
        self.logo_plot(df_info_meancounts, 'information', title=title_base+f'_information_meancounts_{norm_type}_pseudo{self.psuedocount}')
        # self.logo_plot(df_info_mediancounts, 'information', title=title_base+f'_information_mediancounts_pseudo{self.psuedocount}')
        self.logo_plot(df_actualreads, 'counts', title=title_base+f'_counts_actualreads_{norm_type}_pseudo{self.psuedocount}')
        self.logo_plot(df_seqreads, 'counts', title=title_base+f'_counts_seqreads')
        self.logo_plot(df_sumcounts, 'counts', title=title_base+f'_counts_sumcounts_{norm_type}_pseudo{self.psuedocount}')
        self.logo_plot(df_meancounts, 'counts', title=title_base+f'_counts_meancounts_{norm_type}_pseudo{self.psuedocount}')
        # self.logo_plot(df_mediancounts, 'counts', title=title_base+f'_counts_mediancounts_pseudo{self.psuedocount}')

    def logo_plot(self, df, seq_type, title):
        fig, ax = plt.subplots(figsize=(9,3))
        logomaker.Logo(df, color_scheme=self.pal_dict, ax=ax)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xlim([-1, 73])
        ax.set_xticks([])
        ax.set_title(title, fontsize='x-large')
        
        # if self.divergence_plot == True:
        #     ax.spines['left'].set_bounds(-2, 2)
        #     y = -2.25
        #     yys = 0.935
        #     ax.set_ylim([-2.5, 2])
        #     ax.set_yticks([-2, -1, 0, 1, 2])
        #     ax.set_yticklabels(['-2','-1','0','1','2'])
        #     ax.set_ylabel("       Bits", fontsize=12, verticalalignment='center')
        # else:

        ax.spines['left'].set_bounds(0, 2)
        y = -.16
        yys = 0.6
        # ax.set_ylim([-0.25, 2])
        if seq_type == 'information':
            ax.set_ylim([-0.25, 2])
            ax.set_yticks([0, 1, 2])
            ax.set_yticklabels(['0','1','2'])
            ax.set_ylabel("       Bits", fontsize=12, verticalalignment='center')
        elif seq_type == 'counts':
            pass
            # ax.set_yticks([0,2])
            # ax.set_yticklabels(['0%','100%'])
            # ax.set_ylabel("       Frequency", fontsize=12, verticalalignment='center')    
            # if self.cutoff > 0:
            #     ax.plot([-1,73],[self.cutoff*2,self.cutoff*2],linewidth=1,ls='--',color='black')
        
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

        plt.savefig(f'{self.output}/{title}.pdf', bbox_inches='tight')

    def consensus_plot(self, df, seqs, unit):
        # Create a figure with 4 subplots
        fig = plt.figure(figsize=(10, 2))
        # Create 4 subplots with the first 3 stacked on top of each other and the last one spanning all 3 on the right
        gs = gridspec.GridSpec(3, 2, height_ratios=[1, 0.25, 0.25], width_ratios=[1, 0.025])
        # Create a seqlogo in the top subplot
        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0])
        ax3 = plt.subplot(gs[2, 0])
        ax4 = plt.subplot(gs[:, 1])
        # Create a seqlogo in the top subplot
        logomaker.Logo(df, color_scheme=self.pal_dict, ax=ax1)
        ax1.set_ylim([0, 2])
        ax1.set_yticks([0, 1, 2])
        ax1.set_yticklabels(['0','1','2'])
        ax1.set_ylabel("Bits", rotation=0, verticalalignment='center')
        ax1.get_yaxis().set_label_coords(-0.1,0.5)
        ax1.set_xticks([])
        # Create a consensus sequence from the df where the base with the highest information is selected
        consensus = df.apply(lambda x: x.idxmax(), axis=1)
        # Create a matrix from the sequence list where each row is a sequence and each column is a position
        # The matrix is a factor of 1 if the base at that position matches the consensus and 0 if it does not
        df_seq = pd.DataFrame(np.full((len(seqs),len(seqs[0])), 0), columns=[i for i in range(len(seqs[0]))])
        for seq in range(len(seqs)):
            for pos in range(len(seqs[seq])-3): # Remove the last 3 positions from the matrix to remove the CCA tail
                if seqs[seq][pos] == consensus[pos]:
                    df_seq.at[seq,pos] = 1
        # Collapse the matrix into a single row by summing the columns then divide by the number of sequences
        df_seq = df_seq.sum(axis=0)/len(seqs)

        print(consensus)

        # For each position in df_seq equal to 0, set the consensus base to a dash
        for i in range(len(df_seq)):
            if df_seq[i] == 0:
                consensus = consensus[:i] + '-' + consensus[i+1:]

                
        print(consensus)
        # Create a heatmap in the lower subplot
        ax2.imshow([df_seq], cmap='Blues', aspect='auto', interpolation='nearest')
        # Print the consensus horizontally aligned with the heatmap below it
        for i in range(len(consensus)):
            if df_seq[i] > 0.5:
                ax2.text(i, 0, consensus[i], fontsize=10, verticalalignment='center', horizontalalignment='center', color='white')
            else:
                ax2.text(i, 0, consensus[i], fontsize=10, verticalalignment='center', horizontalalignment='center', color='black')
        ax2.set_yticks([])
        ax2.set_ylabel("Consensus", rotation=0, verticalalignment='center')
        ax2.get_yaxis().set_label_coords(-0.1,0.5)
        ax2.set_xticks([])
        # Create a legend in the right subplot for the heatmap
        ax4.imshow(np.linspace(0, 100, 101).reshape(-1, 1)[::-1], cmap='Blues', aspect='auto', interpolation='nearest')
        ax4.set_xticks([])
        ax4.set_yticks([0, 100])
        ax4.set_yticklabels(['100%', '0%'])
        ax4.set_ylabel("Frequency", rotation=90, verticalalignment='center')
        # Plot dashed lines
        for i in [ax1, ax2]:
            i.plot([6.5, 6.5], [i.get_ylim()[0], i.get_ylim()[1]], color='k', linewidth=1, linestyle='--', zorder=0)
            i.plot([8.5, 8.5], [i.get_ylim()[0], i.get_ylim()[1]], color='k', linewidth=1, linestyle='--', zorder=0)
            i.plot([24.5, 24.5], [i.get_ylim()[0], i.get_ylim()[1]], color='k', linewidth=1, linestyle='--', zorder=0)
            i.plot([25.5, 25.5], [i.get_ylim()[0], i.get_ylim()[1]], color='k', linewidth=1, linestyle='--', zorder=0)
            i.plot([42.5, 42.5], [i.get_ylim()[0], i.get_ylim()[1]], color='k', linewidth=1, linestyle='--', zorder=0)
            i.plot([47.5, 47.5], [i.get_ylim()[0], i.get_ylim()[1]], color='k', linewidth=1, linestyle='--', zorder=0)
            i.plot([64.5, 64.5], [i.get_ylim()[0], i.get_ylim()[1]], color='k', linewidth=1, linestyle='--', zorder=0)
        # Create the tRNA scheme in the bottom subplot
        ax3.set_ylim([-1, 1])
        ax3.plot([-0.5, 72.5], [0, 0], color='k', linewidth=1, linestyle='-', zorder=-3)
        ax3.axvspan(8.5, 24.5, color='#ff9999', zorder=-2) # D-Arm/Loop Outter
        ax3.axvspan(12.5, 20.5, color='white', zorder=-1) # D-Arm/Loop Inner
        ax3.axvspan(25.5, 42.5, color='#99ff99', zorder=-2) # A-Arm/Loop Outter
        ax3.axvspan(30.5, 37.5, color='white', zorder=-1) # A-Arm/Loop Inner
        ax3.axvspan(47.5, 64.5, color='#99ffff', zorder=-2) # T-Arm/Loop Outter
        ax3.axvspan(52.5, 59.5, color='white', zorder=-1) # T-Arm/Loop Inner
        ax3.axvspan(-0.5, 6.5, color='#98c0c0', zorder=-2) # Stem 1
        ax3.axvspan(64.5, 71.5, color='#98c0c0', zorder=-2) # Stem 2
        ax3.text(16.5, 0, 'D-Arm/Loop', fontsize=10, verticalalignment='center', horizontalalignment='center', color='k')
        ax3.text(34, 0, 'A-Arm/Loop', fontsize=10, verticalalignment='center', horizontalalignment='center', color='k')
        ax3.text(56, 0, 'T-Arm/Loop', fontsize=10, verticalalignment='center', horizontalalignment='center', color='k')
        # Set the x-axis limits so each plot is aligned
        ax1.set_xlim([-1, 73])
        ax2.set_xlim([-1, 73])
        ax3.set_xlim([-1, 73])
        # Remove the top and right spines for ax1 and ax2
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        # Remove the axis spines for ax3
        ax3.axis('off')
        # Add a title to the figure
        fig.suptitle(f'Consensus Logo {self.coverage_grp}_{unit}', fontsize='x-large')
        # save the figure
        plt.savefig(f'{self.output}/consensus_{self.coverage_grp}_{unit}.pdf', bbox_inches='tight')

if __name__ == '__main__':
    pass
    # parser = argparse.ArgumentParser(
    #     prog='logo_tools.py',
    #     description='Generate seqlogos from anndata.',
    # )

    # parser.add_argument('-i', '--anndata', help='Specify AnnData input', required=True)
    # parser.add_argument('-o', '--output', help='Specify output directory', default='seqlogo', required=False)

    # args = parser.parse_args()

    # # Create output directory if it doesn't exist
    # toolsTG.builder(args.output)

    # adata = ad.read_h5ad(args.anndata)

    # visualizer(adata, args.output)
