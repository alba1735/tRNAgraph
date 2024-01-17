#!/usr/bin/env python3

import numpy as np
import pandas as pd

import logomaker

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

import time

class visualizer():
    def __init__(self, adata, grp, manual_grp, manual_name, pseudocount, logosize, ccatail, pseudogenes, output):
        self.adata = adata
        self.adata = self.adata[:,self.adata.var['positions'] != '-1']
        self.output = output
        self.coverage_grp = grp
        self.manual_grp = manual_grp
        self.manual_name = manual_name
        self.pseudocount = pseudocount
        self.logosize = logosize
        self.ccatail = ccatail
        self.pseudogenes = pseudogenes
        self.pal_dict = {'T':'#ea4335','A':'#34a853','C':'#4285f4', 'G':'#fbbc05', '-':'#ffffff'}
        # Drop CCA tail if specified
        if self.ccatail:
            for i in ['74','75','76']:
                self.adata = self.adata[:,self.adata.var['positions'] != i]
        # Drop pseudogenes if specified, by dropping the rows where obs.pseudogene == 'tRX'
        if self.pseudogenes:
            self.adata = self.adata[self.adata.obs['pseudogene'] != 'tRX']
        # Create a list of the positions from the adata
        adata_pos = self.adata.var[self.adata.var['coverage'] == 'coverage']
        if self.logosize == 'full' or self.logosize == 'noloop':
            self.positions_list = adata_pos['positions'].tolist()
        else:
            self.positions_list = adata_pos[adata_pos['gap'] == False]['positions'].tolist()

    def generate_plots(self):
        if self.manual_grp:
            # Subset the adata to the manual tRNA list
            seq_adata = self.adata[self.adata.obs['trna'].isin(self.manual_grp)]
            # Create a df of the tRNA sequences for the map plot
            df_seqinfo, df_consensus, seq_df = self.build_seqdata(seq_adata)
            # Create the plots
            outputgrp = 'manual_{}'.format(time.strftime("%Y%m%d-%H%M%S"))
            if self.manual_name:
                outputgrp = 'manual_{}'.format(self.manual_name)
            self.consensus_plot(df_seqinfo, df_consensus, outputgrp)
            self.map_plot(seq_df, outputgrp)
        else:
            for i in sorted(self.adata.obs[self.adata.obs[self.coverage_grp].notna()][self.coverage_grp].unique()):
                # Subset the adata to the current coverage group
                seq_adata = self.adata[self.adata.obs[self.coverage_grp] == i]
                # Create a df of the tRNA sequences for the map plot
                df_seqinfo, df_consensus, seq_df = self.build_seqdata(seq_adata)
                # Create the plots
                self.consensus_plot(df_seqinfo, df_consensus, i)
                self.map_plot(seq_df, i)

    def build_seqdata(self, seq_adata):
        trna_list = seq_adata.obs.trna.unique().tolist()
        seq_list = seq_adata.obs.refseq.values.tolist()
        # Create a dictionary of the tRNA sequences for the map plot
        seq_dict = {}
        for j in trna_list:
            if self.logosize == 'full' or self.logosize == 'noloop':
                if self.ccatail:
                    seq_dict[j] = self.adata.obs[self.adata.obs['trna'] == j].refseq_full.unique()[0][1:-3]
                else:
                    seq_dict[j] = self.adata.obs[self.adata.obs['trna'] == j].refseq_full.unique()[0][1:]
            else:
                if self.ccatail:
                    seq_dict[j] = self.adata.obs[self.adata.obs['trna'] == j].refseq.unique()[0][:-3]
                else:
                    seq_dict[j] = self.adata.obs[self.adata.obs['trna'] == j].refseq.unique()[0]
        # Create a df of the tRNA sequences for the map plot from the dictionary
        seq_df = pd.DataFrame(columns=['position']+trna_list)
        for j in range(len(seq_dict[trna_list[0]])):
            seq_df.loc[j] = [self.positions_list[j]]+[seq_dict[k][j] for k in trna_list]
        seq_df['match'] = [True]*len(seq_df)
        # Iterate through the df and check if the base at each position matches the base in the first tRNA for that position
        for j, row in seq_df.iterrows():
            for k in trna_list:
                if row[k] != row[trna_list[0]]:
                    seq_df.at[j,'match'] = False
                    break
        # Remove the extension loop from the df if specified
        if self.logosize == 'noloop':
            # Subset the adata to the extension loop
            seq_adata = seq_adata[:,seq_adata.var['location'] != 'extensionloop']
            # Subset the seq_df to the extension loop
            seq_df['extensionloop'] = self.adata.var[self.adata.var['coverage'] == 'coverage']['location'].isin(['extensionloop']).tolist()
            seq_df = seq_df[seq_df['extensionloop'] == False]
            seq_df = seq_df.drop(['extensionloop'], axis=1)
            seq_df = seq_df.reset_index(drop=True)
        if self.logosize == 'sprinzl':
            seq_adata = seq_adata[:,seq_adata.var['gap'] == False]
        # Create an information matrix from the sequence reads
        # df_seqinfo, df_consensus = self.build_seqdata(seq_adata)
        # Create a matrix of 0s for the actual base reads
        df_seqreads = pd.DataFrame(np.full((seq_adata[:,seq_adata.var['coverage'] == 'coverage'].shape[1],5),0), columns=['A','C','G','T','M']) 
        # Iterate through the reads in A,G,C,T order and add the counts to the matrix from the adata
        for k,base in {'A':'adenines', 'C':'cytosines', 'G':'guanines', 'T':'thymines', 'M':'mismatchedbases'}.items():
            base_list = seq_adata[:,seq_adata.var['coverage'] == base]
            df_seqreads[k] = base_list.X.mean(axis=0)
            # df_seqreads[k] = base_list.X.sum(axis=0)
        # Drop the M column into its own df
        mismatchedbases = df_seqreads['M'].tolist()
        df_seqreads = df_seqreads.drop(['M'], axis=1)
        # Add a pseudocount to the matrix
        df_seqreads += self.pseudocount
        # create df from seq_df position 
        df_consensus = pd.DataFrame(columns=['position'])
        df_consensus['position'] = seq_adata[:,seq_adata.var['coverage'] == 'coverage'].var['positions'].tolist()
        # Create a consensus sequence from the df where the most common base is selected from seqdf column 1:-1 and if under 50% set to dash
        df_consensus['ref_score'] = seq_df.iloc[:,1:-1].apply(lambda x: x.value_counts().max()/x.value_counts().sum(), axis=1)
        df_consensus['ref_base'] = seq_df.iloc[:,1:-1].apply(lambda x: x.value_counts().idxmax(), axis=1)
        df_consensus['ref_base'] = df_consensus.apply(lambda x: ' ' if x['ref_score'] < 0.35 else x['ref_base'], axis=1)
        df_consensus['actual_score'] = df_seqreads.apply(lambda x: x.max()/x.sum(), axis=1)
        # Replace NaN values with 0 to account for 0 reads at a position giving a zero division error
        df_seqreads = df_seqreads.fillna(0)
        df_consensus['actual_base'] = df_seqreads.apply(lambda x: x.idxmax(), axis=1)
        # If ref base is a dash and actual base is not a dash, set the actual base to a dash
        df_consensus['actual_base'] = df_consensus.apply(lambda x: '-' if x['ref_base'] == '-' and x['actual_base'] != '-' else x['actual_base'], axis=1)
        df_consensus['actual_base'] = df_consensus.apply(lambda x: ' ' if x['actual_score'] < 0.35 else x['actual_base'], axis=1)
        # Add the mismatched bases to the df scaled by the max value
        df_consensus['mismatchedbases'] = mismatchedbases
        df_consensus['mismatchedbases'] = df_consensus.apply(lambda x: x['mismatchedbases']/df_consensus['mismatchedbases'].max(), axis=1)
        # Add the index to the df as plotposition
        df_consensus['plotposition'] = df_consensus.index
        # Convert the counts to information
        df_seqinfo = logomaker.transform_matrix(df_seqreads, from_type='counts', to_type='information')
    
        return df_seqinfo, df_consensus, seq_df

    def consensus_plot(self, df_seqinfo, df_consensus, unit):
        # Create a figure with 4 subplots
        fig = plt.figure(figsize=(20, 3))
        # Create 4 subplots with the first 3 stacked on top of each other and the last one spanning all 3 on the right
        gs = gridspec.GridSpec(3, 2, height_ratios=[1, 0.75, 0.25], width_ratios=[1, 0.025])
        # Create axes
        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[1, 0])
        ax3 = plt.subplot(gs[2, 0])
        ax4 = plt.subplot(gs[:, 1])
        # Create a seqlogo in the top subplot
        logomaker.Logo(df_seqinfo, color_scheme=self.pal_dict, ax=ax1)
        ax1.set_ylim([0, 2])
        ax1.set_yticks([0, 1, 2])
        ax1.set_yticklabels(['0','1','2'])
        ax1.set_ylabel("Bits", rotation=90, verticalalignment='center')
        ax1.get_yaxis().set_label_coords(-0.1,0.5)
        ax1.set_xticks([])
        # Create a heatmap in the lower subplot
        # scores_df = pd.DataFrame([df_consensus['actual_score'], df_consensus['ref_score'], df_consensus['mismatchedbases']])
        scores_df = pd.DataFrame([df_consensus['actual_score'], df_consensus['ref_score']])
        sns.heatmap(scores_df, cmap='Blues', cbar=False, square=True, linewidths=0, linecolor='black', ax=ax2)
        # Print the consensus horizontally aligned with the heatmap below it
        for i in df_consensus.iterrows():
            ascore_color, rscore_color, fweight = 'black', 'black', 'normal'
            if i[1]['actual_score'] > 0.7:
                ascore_color = 'white'
            if i[1]['ref_score'] > 0.7:
                rscore_color = 'white'
            if i[1]['actual_base'] != i[1]['ref_base']:
                if i[1]['actual_base'] != ' ' or i[1]['ref_base'] != ' ':
                    ascore_color = 'red'
                    rscore_color = 'red'
                    fweight='bold'
            ax2.text(i[0]+0.1, 0.55, i[1]['actual_base'], fontsize=12, verticalalignment='center', color=ascore_color, weight=fweight)
            ax2.text(i[0]+0.1, 1.55, i[1]['ref_base'], fontsize=12, verticalalignment='center', color=rscore_color, weight=fweight)
        ax2.set_yticks([0.5, 1.5])
        ax2.set_yticklabels(['Reads', 'Reference'])
        ax2.tick_params(axis='y', which='both', left=False, right=False, labelleft=True)
        ax2.set_ylabel("Consensus", rotation=90, verticalalignment='center')
        ax2.get_yaxis().set_label_coords(-0.1,0.5)
        ax2.set_xticks([])
        # Create a legend in the right subplot for the heatmap
        ax4.imshow(np.linspace(0, 100, 101).reshape(-1, 1)[::-1], cmap='Blues', aspect='auto', interpolation='nearest')
        ax4.set_xticks([])
        ax4.set_yticks([0, 100])
        ax4.set_yticklabels(['100%', '0%'])
        ax4.set_ylabel("Frequency", rotation=90, verticalalignment='center')  
        # Create the tRNA scheme in the bottom subplot
        ax3.set_ylim([-1, 1])
        ax3.set_xlim(ax2.get_xlim())
        ax3.axhspan(-0.125, 0.125, color='black', zorder=-3)
        # Plot the tRNA scheme
        plotposition = (self.position_swap(df_consensus, '10'), self.position_swap(df_consensus, '26'))
        ax3.axvspan(plotposition[0], plotposition[1], color='#ff9999', zorder=-2) # D-Arm/Loop Outer
        plotposition = (self.position_swap(df_consensus, '14'), self.position_swap(df_consensus, '22'))
        ax3.axvspan(plotposition[0], plotposition[1], color='white', alpha=0.8, zorder=-1) # D-Arm/Loop Inner
        plotposition = (self.position_swap(df_consensus, '27'), self.position_swap(df_consensus, '44'))
        ax3.axvspan(plotposition[0], plotposition[1], color='#99ff99', zorder=-2) # A-Arm/Loop Outter
        plotposition = (self.position_swap(df_consensus, '32'), self.position_swap(df_consensus, '39'))
        ax3.axvspan(plotposition[0], plotposition[1], color='white', alpha=0.8, zorder=-1) # A-Arm/Loop Inner
        plotposition = (self.position_swap(df_consensus, '49'), self.position_swap(df_consensus, '66'))
        ax3.axvspan(plotposition[0], plotposition[1], color='#99ffff', zorder=-2) # T-Arm/Loop Outter
        plotposition = (self.position_swap(df_consensus, '54'), self.position_swap(df_consensus, '61'))
        ax3.axvspan(plotposition[0], plotposition[1], color='white', alpha=0.8, zorder=-1) # T-Arm/Loop Inner
        plotposition = (ax3.get_xlim()[0], self.position_swap(df_consensus, '8'))
        ax3.axvspan(plotposition[0], plotposition[1], color='#98c0c0', zorder=-2) # Stem 1
        if not self.ccatail:
            plotposition = (self.position_swap(df_consensus, '66'), ax3.get_xlim()[1]-3)
        else:
            plotposition = (self.position_swap(df_consensus, '66'), ax3.get_xlim()[1])
        ax3.axvspan(plotposition[0], plotposition[1], color='#98c0c0', zorder=-2) # Stem 2
        # Plot the tRNA scheme text
        ax3.text(self.position_swap(df_consensus, '18'), 0, 'D-Arm/Loop', fontsize=10, verticalalignment='center', horizontalalignment='center', color='k')
        ax3.text(self.position_swap(df_consensus, '35'), 0, 'A-Arm/Loop', fontsize=10, verticalalignment='center', horizontalalignment='center', color='k')
        ax3.text(self.position_swap(df_consensus, '57'), 0, 'T-Arm/Loop', fontsize=10, verticalalignment='center', horizontalalignment='center', color='k')
        # Plot dashed lines
        for i in [ax1, ax2, ax3]:
            for j in ['8','10','26','27','44','49','66']:
                plotposition = self.position_swap(df_consensus, j)
                if i == ax1:
                    plotposition -= 0.5
                i.plot([plotposition, plotposition], [i.get_ylim()[0], i.get_ylim()[1]], color='k', linewidth=1, linestyle='--', zorder=1)
        # Remove the top and right spines for ax1 and ax2
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        # Remove the axis spines for ax3
        ax3.axis('off')
        # Add a title to the figure and save the figure
        if self.manual_grp:
            fig.suptitle(f'Consensus Logo {unit}', fontsize='x-large')
            plt.savefig(f'{self.output}/consensus_{unit}.pdf', bbox_inches='tight')
        else:
            fig.suptitle(f'Consensus Logo {self.coverage_grp}_{unit}', fontsize='x-large')
            plt.savefig(f'{self.output}/consensus_{self.coverage_grp}_{unit}.pdf', bbox_inches='tight')
        plt.close()

    def map_plot(self, seq_df, unit):
        trna_list = seq_df.columns[1:-1].tolist()
        # Make a all white heatmap in the shape of the sequence
        fig = plt.figure(figsize=(18,len(trna_list)))

        seq_hm = np.zeros((len(trna_list), len(seq_df)))
        sns.heatmap(seq_hm, cmap='Greys', cbar=False, square=True, linewidths=0, linecolor='black')
        # Populate the heatmap boxes with the sequence matching the color from the self.pal_dict or gray if match is True
        tick_pos, tick_names = [],[]

        for i, row in seq_df.iterrows():
            for j in trna_list:
                if row['match']:
                    plt.text(i+0.5, trna_list.index(j)+0.5, row[j], ha='center', va='center', fontsize=12, color='lightgray')
                else:
                    plt.text(i+0.5, trna_list.index(j)+0.5, row[j], ha='center', va='center', fontsize=16, color=self.pal_dict[row[j]], weight='bold')
                    tick_pos.append(i+0.5)
                    tick_names.append(row['position'])

        # Add the xticks
        plt.xticks(tick_pos, tick_names, rotation=90)
        # Add the sequence list to the yticks
        plt.yticks(np.arange(len(trna_list))+0.5, trna_list, rotation=0)
        plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=True)
        plt.xlabel('Position')
        plt.title('Sequence Logo')
        # Add a title to the figure and save the figure
        if self.manual_grp:
            # fig.suptitle(f'Map {unit}', fontsize='x-large')
            plt.savefig(f'{self.output}/map_{unit}.pdf', bbox_inches='tight')
        else:
            # fig.suptitle(f'Map {self.coverage_grp}_{unit}', fontsize='x-large')
            plt.savefig(f'{self.output}/map_{self.coverage_grp}_{unit}.pdf', bbox_inches='tight')
        plt.close()

    def position_swap(self, df, pos):
        # Return the plotposition for a given position
        return df[df['position'] == pos]['plotposition'].values[0]
        
if __name__ == '__main__':
    pass

    # def build_seqdata(self, seq_adata):
    #     # Create a matrix of 0s for the sequence reads
    #     df_seqreads = pd.DataFrame(np.full((seq_adata[:,seq_adata.var['coverage'] == 'coverage'].shape[1],4),0), columns=['A','C','G','T'])
    #     # Iterate through the reads in A,G,C,T order and add the counts to the matrix from the adata
    #     for k,base in {'A':'adenines', 'C':'cytosines', 'G':'guanines', 'T':'thymines'}.items():
    #         base_list = seq_adata[:,seq_adata.var['coverage'] == base]
    #         # multiple the base_list by the deseq2_sizefactor to unnormalize the counts
    #         # base_list.X = base_list.X * base_list.obs['deseq2_sizefactor'].values.reshape(-1,1)
    #         df_seqreads[k] = base_list.X.mean(axis=0)
    #         # df_seqreads[k] = base_list.X.sum(axis=0)
        
    #     # Create a consensus sequence from the df where the base with the highest information is selected
    #     consensus = df_seqreads.apply(lambda x: x.idxmax(), axis=1)
    #     # Create a list for the consensus scores
    #     consensus_score = []
    #     for row in df_seqreads.iterrows():
    #         # For each row in the df, calculate the consensus score by dividing the highest value by the sum of all values
    #         consensus_score.append(row[1][row[1].idxmax()]/row[1].sum())
    #     # Combine the consensus scores with the consensus sequence into a single df
    #     df_consensus = pd.concat([pd.DataFrame(consensus_score), consensus], axis=1)
    #     df_consensus.columns = ['score', 'base']
    #     df_consensus = df_consensus.fillna(0)
    #     # Add the position column to the df
    #     df_consensus['position'] = seq_adata[:,seq_adata.var['coverage'] == 'coverage'].var['positions'].tolist()
    #     # For each position in df_combine equal to 0, set the consensus base to a dash
    #     for i in range(len(df_consensus['score'])):
    #         if df_consensus['score'][i] <= 0.5:
    #             df_consensus.at[i,'base'] = '-'

    #     for i in df_consensus.iterrows():
    #         print(i[1]['base'], i[1]['score'], i[1]['position'])

    #     # Add a pseudocount to the matrix
    #     df_seqreads += self.pseudocount

    #     # # Iterate through the sequence list and add the counts to the matrix to generate a count matrix
    #     # for seq in seq_list:
    #     #     for i in range(len(seq)):
    #     #         if seq[i] != '-':
    #     #             df_seqreads.at[i,seq[i]] += 1
    #     # # Remove the last 3 positions from the matrix to remove the CCA tail
    #     # # df_seqreads = df_seqreads.iloc[:-3,:]

    #     # # Create a consensus sequence from the df where the base with the highest information is selected
    #     # consensus = df_seqinfo.apply(lambda x: x.idxmax(), axis=1)
    #     # # Create a matrix from the sequence list where each row is a sequence and each column is a position
    #     # # The matrix is a factor of 1 if the base at that position matches the consensus and 0 if it does not
    #     # df_seq = pd.DataFrame(np.full((len(seqs),len(seqs[0])), 0), columns=[i for i in range(len(seqs[0]))])
    #     # for seq in range(len(seqs)):
    #     #     for pos in range(len(seqs[seq])-3): # Remove the last 3 positions from the matrix to remove the CCA tail
    #     #         if seqs[seq][pos] == consensus[pos]:
    #     #             df_seq.at[seq,pos] = 1
    #     # # Collapse the matrix into a single row by summing the columns then divide by the number of sequences
    #     # df_seq = df_seq.sum(axis=0)/len(seqs)
    #     # # combine the df_seq and consensus into a single df
    #     # df_combine = pd.concat([df_seq, consensus], axis=1)
    #     # df_combine.columns = ['score', 'base']
    #     # # For each position in df_combine equal to 0, set the consensus base to a dash
    #     # for i in range(len(df_combine['score'])):
    #     #     if df_combine['score'][i] == 0:
    #     #         df_combine.at[i,'base'] = '-'
    
    #     # Convert the counts to information
    #     df_seqinfo = logomaker.transform_matrix(df_seqreads, from_type='counts', to_type='information')

    #     return df_seqinfo, df_consensus