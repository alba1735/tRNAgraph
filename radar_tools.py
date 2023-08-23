#!/usr/bin/env python3

import numpy as np
import pandas as pd
import anndata as ad
import argparse

import directory_tools

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns

def visualizer(adata, radargrp, colormap, output):
    '''
    Generate radar plots for each sample in an AnnData object.
    '''
    # Amino acid dictionary for testing in hg38
    #amino_dict = {'Ala':['AGC','CGC','TGC'], 'Arg':['ACG','CCG','CCT','TCG','TCT'], 'Gly':['CCC','GCC','TCC'], 'Ile':['AAT','GAT','TAT'], 'Leu':['AAG','CAG','TAG','CAA','TAA'], 'Pro':['AGG','CGG','TGG'], 'Ser':['AGA','CGA','GCT','TGA'], 'Thr':['AGT','CGT','TGT'], 'Val':['AAC','CAC','TAC']}
    df = adata.uns['anticodon_counts'].copy()
    df = pd.DataFrame(adata.obs, columns=['iso', 'amino', radargrp, 'nreads_total_unique_raw'])
    # Drop rows where iso is NNN
    df = df[df['iso']!='NNN']
    # Create dict from all iso and amino where each amino is a key and each value is a list of iso that have that amino
    amino_dict = {}
    for amino in df['amino'].unique():
        amino_dict[amino] = df[df['amino']==amino]['iso'].unique().tolist()
    df = pd.pivot_table(df, values='nreads_total_unique_raw', index=['iso'], columns=[radargrp], aggfunc=np.sum)
    # Only create radar plots for amino acids with at least 3 isoforms
    for aa, codons, in amino_dict.items():
        if len(codons) >= 3:
            # Subset dataframe to only include codons for the current amino acid
            tdf = df[df.index.isin(codons)]
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
            # Set colormap if specified
            if colormap != None:
                colormap = {k:v if v[0]!='#' else mplcolors.to_rgb(v) for k,v in colormap.items()}
                for v in tdf.columns:
                    if v not in colormap:
                        print(f'Color {v} not found in colormap. Using default colors instead.')
                        colormap = None
                        break
            # Plot data
            for i in tdf.columns.values:
                v = tdf[i].values.tolist()
                v = np.concatenate((v, v[:1]))
                a = angles + angles[:1]
                # Plot data for each anticodon group
                if colormap != None:
                    ax.plot(a, v, linewidth=1.5, linestyle='solid', label=i, color=colormap[i])
                    ax.fill(a, v, alpha=0.5/len(tdf.columns.values), color=colormap[i])
                else:
                    pal = dict(zip(tdf.columns, sns.color_palette('husl', len(tdf.columns))))
                    ax.plot(a, v, linewidth=1.5, linestyle='solid', label=i, color=pal[i])
                    ax.fill(a, v, alpha=0.5/len(tdf.columns.values), color=pal[i])
            # Capatilize the legend and move the legend outside the plot and remove the border around it
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles=handles, labels=[x.capitalize() for x in labels])
            ax.legend(loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
            ax.legend_.set_title(radargrp.capitalize())
            # Give plot a title
            plt.title(f'{aa} Codon Pool Percent Expression Change')
            # Save plot
            plt.savefig(f'{output}/{aa}_radar_by_{radargrp}.pdf', bbox_inches='tight')
            print(f'Plot saved to {output}/{aa}_radar_by_{radargrp}.pdf')
            plt.close()

    # Create radar plot for all amino acids
    # Normalize counts to 100
    tdf = 100*df/df.sum(axis=0)
    # Drop rows where the mean of the row comprises less than 1% of the total reads
    tdf = tdf[tdf.min(axis=1)>1]
    tdf = tdf.sort_index()
    # Take the log2 value of the dataframe
    # tdf = np.log2(tdf)
    # Replace negative values with 0
    tdf[tdf<0] = 0
    # Set angles for each category
    angles = np.linspace(0, 2*np.pi, num=len(list(tdf.T))+1)[:-1]
    angles = angles.tolist()
    # Create figure
    fig, ax = plt.subplots(figsize=(12,12), subplot_kw=dict(polar=True))
    # Set theta zero location so that anticodons are in the correct position
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    # Set plot parameters
    plt.xticks(angles, list(tdf.T), color='black', size=12)
    plt.yticks(color='dimgrey', size=8)
    ax.set_rlabel_position(0)
    ax.xaxis.grid(False)
    # Set colormap if specified
    if colormap != None:
        colormap = {k:v if v[0]!='#' else mplcolors.to_rgb(v) for k,v in colormap.items()}
        for v in tdf.columns:
            if v not in colormap:
                print(f'Color {v} not found in colormap. Using default colors instead.')
                colormap = None
                break
    # Plot data
    for i in tdf.columns.values:
        v = tdf[i].values.tolist()
        v = np.concatenate((v, v[:1]))
        a = angles + angles[:1]
        # Plot data for each anticodon group
        if colormap != None:
            ax.plot(a, v, linewidth=1.5, linestyle='solid', label=i, color=colormap[i])
            ax.fill(a, v, alpha=0.5/len(tdf.columns.values), color=colormap[i])
        else:
            pal = dict(zip(tdf.columns, sns.color_palette('husl', len(tdf.columns))))
            ax.plot(a, v, linewidth=1.5, linestyle='solid', label=i, color=pal[i])
            ax.fill(a, v, alpha=0.5/len(tdf.columns.values), color=pal[i])
    # Capatilize the legend and move the legend outside the plot and remove the border around it
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=[x.capitalize() for x in labels])
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
    ax.legend_.set_title(radargrp.capitalize())
    # Give plot a title
    plt.title(f'Codon Pool Percent Expression Change')
    # Save plot
    plt.savefig(f'{output}/all_gt1percent_radar_by_{radargrp}.pdf', bbox_inches='tight')
    print(f'Plot saved to {output}/all_radar_by_{radargrp}.pdf')
    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='radar_tools.py',
        description='Generate radial plots for each sample in an AnnData object.',
    )

    parser.add_argument('-i', '--anndata', help='Specify AnnData input', required=True)
    parser.add_argument('--radargrp', help='Specify AnnData column to group by (default: group) (optional)', default='group', required=False)
    parser.add_argument('-o', '--output', help='Specify output directory', default='radial', required=False)
    parser.add_argument('--colormap', help='Specify a colormap for coverage plots (optional)', default=None)

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    directory_tools.builder(args.output)

    adata = ad.read_h5ad(args.anndata)

    visualizer(adata, args.radargrp, args.colormap, args.output)