#!/usr/bin/env python3

import numpy as np
import pandas as pd
# import anndata as ad
# import argparse

import logomaker

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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

        # self.cutoff = 0 
        # self.cov_mask = False
        # self.divergence_plot = False

    def generate_plots(self):
        for i in sorted(self.adata.obs[self.coverage_grp].values.unique()):
            seq_list = [i for i in self.adata.obs[self.adata.obs[self.coverage_grp] == i].refseq.values]
            seq_adata = self.adata[self.adata.obs[self.coverage_grp] == i]
            # Remove gaps from the adata
            seq_adata = seq_adata[:,seq_adata.var['gap'] == False]
            seq_size = seq_adata[:,seq_adata.var['coverage'] == 'adenines'].shape[1]
            # Create an information matrix from the sequence reads
            df_seqinfo, df_consensus = self.build_seqdata(seq_list, seq_adata, seq_size)
            self.consensus_plot(df_seqinfo, df_consensus, seq_list, i)


    def build_seqdata(self, seq_list, seq_adata, seq_size):
        # Create a matrix of 0s for the sequence reads
        df_seqreads = pd.DataFrame(np.full((seq_size,4),0), columns=['A','C','G','T'])
        # Iterate through the reads in A,G,C,T order and add the counts to the matrix from the adata
        for k,base in {'A':'adenines', 'C':'cytosines', 'G':'guanines', 'T':'thymines'}.items():
            base_list = seq_adata[:,seq_adata.var['coverage'] == base]
            # multiple the base_list by the deseq2_sizefactor to unnormalize the counts
            # base_list.X = base_list.X * base_list.obs['deseq2_sizefactor'].values.reshape(-1,1)
            df_seqreads[k] = base_list.X.mean(axis=0)
            # df_seqreads[k] = base_list.X.sum(axis=0)


        # Add a psuedocount to the matrix
        df_seqreads += self.psuedocount

        # Create a consensus sequence from the df where the base with the highest information is selected
        consensus = df_seqreads.apply(lambda x: x.idxmax(), axis=1)
        # Create a list for the consensus scores
        consensus_score = []
        for row in df_seqreads.iterrows():
            # For each row in the df, calculate the consensus score by dividing the highest value by the sum of all values
            consensus_score.append(row[1][row[1].idxmax()]/row[1].sum())
        # Combine the consensus scores with the consensus sequence into a single df
        df_consensus = pd.concat([pd.DataFrame(consensus_score), consensus], axis=1)
        df_consensus.columns = ['score', 'base']
        df_consensus['position'] = df_seqreads.index + 1
        # For each position in df_combine equal to 0, set the consensus base to a dash
        for i in range(len(df_consensus['score'])):
            if df_consensus['score'][i] == 0:
                df_consensus.at[i,'base'] = '-'

        # # Iterate through the sequence list and add the counts to the matrix to generate a count matrix
        # for seq in seq_list:
        #     for i in range(len(seq)):
        #         if seq[i] != '-':
        #             df_seqreads.at[i,seq[i]] += 1
        # # Remove the last 3 positions from the matrix to remove the CCA tail
        # # df_seqreads = df_seqreads.iloc[:-3,:]

        # # Create a consensus sequence from the df where the base with the highest information is selected
        # consensus = df_seqinfo.apply(lambda x: x.idxmax(), axis=1)
        # # Create a matrix from the sequence list where each row is a sequence and each column is a position
        # # The matrix is a factor of 1 if the base at that position matches the consensus and 0 if it does not
        # df_seq = pd.DataFrame(np.full((len(seqs),len(seqs[0])), 0), columns=[i for i in range(len(seqs[0]))])
        # for seq in range(len(seqs)):
        #     for pos in range(len(seqs[seq])-3): # Remove the last 3 positions from the matrix to remove the CCA tail
        #         if seqs[seq][pos] == consensus[pos]:
        #             df_seq.at[seq,pos] = 1
        # # Collapse the matrix into a single row by summing the columns then divide by the number of sequences
        # df_seq = df_seq.sum(axis=0)/len(seqs)
        # # combine the df_seq and consensus into a single df
        # df_combine = pd.concat([df_seq, consensus], axis=1)
        # df_combine.columns = ['score', 'base']
        # # For each position in df_combine equal to 0, set the consensus base to a dash
        # for i in range(len(df_combine['score'])):
        #     if df_combine['score'][i] == 0:
        #         df_combine.at[i,'base'] = '-'
    
        # Convert the counts to information
        df_seqinfo = logomaker.transform_matrix(df_seqreads, from_type='counts', to_type='information')

        return df_seqinfo, df_consensus

    def consensus_plot(self, df_seqinfo, df_consensus, seqs, unit):
        # Create a figure with 4 subplots
        fig = plt.figure(figsize=(10, 2))
        # Create 4 subplots with the first 3 stacked on top of each other and the last one spanning all 3 on the right
        gs = gridspec.GridSpec(3, 2, height_ratios=[1, 0.25, 0.25], width_ratios=[1, 0.025])
        # Create a seqlogo in the top subplot
        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0])
        ax3 = plt.subplot(gs[2, 0])
        ax4 = plt.subplot(gs[:, 1])
        # Create a seqlogo in the top subplot - Account for CCA tail
        logomaker.Logo(df_seqinfo.iloc[:-3,:], color_scheme=self.pal_dict, ax=ax1)
        ax1.set_ylim([0, 2])
        ax1.set_yticks([0, 1, 2])
        ax1.set_yticklabels(['0','1','2'])
        ax1.set_ylabel("Bits", rotation=0, verticalalignment='center')
        ax1.get_yaxis().set_label_coords(-0.1,0.5)
        ax1.set_xticks([])
        # Subset the df_consensus to remove the last 3 positions from the matrix to remove the CCA tail
        df_consensus = df_consensus.iloc[:-3,:]
        # Create a heatmap in the lower subplot
        ax2.imshow([df_consensus['score']], cmap='Blues', aspect='auto', interpolation='nearest')
        # Print the consensus horizontally aligned with the heatmap below it
        for i in df_consensus.iterrows():
            if i[1]['score'] > 0.7:
                ax2.text(i[0], 0, i[1]['base'], fontsize=10, verticalalignment='center', horizontalalignment='center', color='white')
            else:
                ax2.text(i[0], 0, i[1]['base'], fontsize=10, verticalalignment='center', horizontalalignment='center', color='black')
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
        plt.close()

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
