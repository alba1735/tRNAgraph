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

def visualizer(adata, output):
    '''
    Generate radar plots for each sample in an AnnData object.
    '''
    amino_dict = {'Ala':['AGC','CGC','TGC'], 'Arg':['ACG','CCG','CCT','TCG','TCT'], 'Gly':['CCC','GCC','TCC'], 'Ile':['AAT','GAT','TAT'], 'Leu':['AAG','CAG','TAG','CAA','TAA'], 'Pro':['AGG','CGG','TGG'], 'Ser':['AGA','CGA','GCT','TGA'], 'Thr':['AGT','CGT','TGT'], 'Val':['AAC','CAC','TAC']}
    df = adata.uns['anticodon_counts'].copy()

    # create dict from sample names to group names
    # sample_dict = dict(zip(adata.obs['sample'], adata.obs['group']))
    # combine columns in tdf where the sample names belong to the same group
    # df = df.groupby(sample_dict, axis=1).sum()

    # for aa, codons, in amino_dict.items():
    #     # Subset dataframe to only include codons for the current amino acid
    #     tdf = df[df.index.isin(codons)]

    #     # print('raw counts from -anticodoncounts.txt:')
    #     # print(tdf)
    #     # print('normalized per column between 0-100 summed reads from -anticodoncounts.txt:')
    #     # print(100*tdf/tdf.sum(axis=0))
    #     # break

    #     # Normalize counts to 100
    #     tdf = 100*tdf/tdf.sum(axis=0)
    #     tdf = tdf.sort_index()
    #     # Replace negative values with 0
    #     tdf[tdf<0] = 0
    #     # Set angles for each category
    #     angles = np.linspace(0, 2*np.pi, num=len(list(tdf.T))+1)[:-1]
    #     angles = angles.tolist()
    #     # Create figure
    #     fig, ax = plt.subplots(figsize=(6,6), subplot_kw=dict(polar=True))
    #     # Set theta zero location so that anticodons are in the correct position
    #     ax.set_theta_zero_location('N')
    #     ax.set_theta_direction(-1)
    #     # Set plot parameters
    #     plt.xticks(angles, list(tdf.T), color='black', size=12)
    #     plt.yticks(color='dimgrey', size=8)
    #     ax.set_rlabel_position(0)
    #     ax.xaxis.grid(False)
    #     # Set color palette
    #     pal = dict(zip(tdf.columns, sns.color_palette('husl', len(tdf.columns))))
    #     for i in tdf.columns.values:
    #         v = tdf[i].values.tolist()
    #         v = np.concatenate((v, v[:1]))
    #         a = angles + angles[:1]
    #         # Plot data for each anticodon group
    #         ax.plot(a, v, linewidth=1.5, linestyle='solid', label=i, color=pal[i])
    #         ax.fill(a, v, alpha=0.5/len(tdf.columns.values), color=pal[i])
    #     # Capatilize the legend and move the legend outside the plot and remove the border around it
    #     handles, labels = ax.get_legend_handles_labels()
    #     ax.legend(handles=handles, labels=[x.capitalize() for x in labels])
    #     ax.legend(loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
    #     ax.legend_.set_title('Groups')
    #     # Give plot a title
    #     plt.title(f'{aa} Codon Pool Percent Expression Change')
    #     # Save plot
    #     plt.savefig(f'{output}/{aa}_radar_anticodoncounts.pdf', bbox_inches='tight')
    #     print(f'Plot saved to {output}/{aa}_radar_anticodoncounts.pdf')
    #     plt.close()

    df = pd.DataFrame(adata.obs, columns=['iso', 'group', 'nreads_total_unique_raw'])
    df = pd.pivot_table(df, values='nreads_total_unique_raw', index=['iso'], columns=['group'], aggfunc=np.sum)

    for aa, codons, in amino_dict.items():
        # Subset dataframe to only include codons for the current amino acid
        tdf = df[df.index.isin(codons)]

        # print('raw counts from -trnauniquecounts.txt summed along (fiveprime, threeprime, whole, and other):')
        # print(tdf)
        # print('normalized per column between 0-100 summed reads from -trnauniquecounts.txt:')
        # print(100*tdf/tdf.sum(axis=0))
        # break

        # Normalize counts to 100
        tdf = 100*tdf/tdf.sum(axis=0)
        tdf = tdf.sort_index()
        # Replace negative values with 0
        tdf[tdf<0] = 0
        # Set angles for each category
        angles = np.linspace(0, 2*np.pi, num=len(list(tdf.T))+1)[:-1]
        angles = angles.tolist()
        # Create figure
        fig, ax = plt.subplots(figsize=(6,6), subplot_kw=dict(polar=True))
        # Set theta zero location so that anticodons are in the correct position
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        # Set plot parameters
        plt.xticks(angles, list(tdf.T), color='black', size=12)
        plt.yticks(color='dimgrey', size=8)
        ax.set_rlabel_position(0)
        ax.xaxis.grid(False)
        # Set color palette
        pal = dict(zip(tdf.columns, sns.color_palette('husl', len(tdf.columns))))
        for i in tdf.columns.values:
            v = tdf[i].values.tolist()
            v = np.concatenate((v, v[:1]))
            a = angles + angles[:1]
            # Plot data for each anticodon group
            ax.plot(a, v, linewidth=1.5, linestyle='solid', label=i, color=pal[i])
            ax.fill(a, v, alpha=0.5/len(tdf.columns.values), color=pal[i])
        # Capatilize the legend and move the legend outside the plot and remove the border around it
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=handles, labels=[x.capitalize() for x in labels])
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
        ax.legend_.set_title('Groups')
        # Give plot a title
        plt.title(f'{aa} Codon Pool Percent Expression Change')
        # Save plot
        plt.savefig(f'{output}/{aa}_radar.pdf', bbox_inches='tight')
        print(f'Plot saved to {output}/{aa}_radar.pdf')
        plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='radar_tools.py',
        description='Generate radial plots for each sample in an AnnData object.',
    )

    parser.add_argument('-i', '--anndata', help='Specify AnnData input', required=True)
    parser.add_argument('-o', '--output', help='Specify output directory', default='radial', required=False)

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    directory_tools.builder(args.output)

    adata = ad.read_h5ad(args.anndata)

    visualizer(adata, args.output)