#!/usr/bin/env python3

import seaborn as sns
import numpy as np
import pandas as pd

import toolsTG

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


def visualizer(adata, grp, readtypes, cutoff, heatbound, heatsubplots, output, threaded=False, config_name='default', overwrite=False):
    '''
    Generate heatmap visualizations for each group in an AnnData object.
    '''
    if grp not in adata.obs.columns:
        raise ValueError('Specified group not found in AnnData object.')

    # Create an empty df to store the log2FC values for each group so a combined heatmap can be generated as well
    df_combine = pd.DataFrame()
    # Create a heatmap for each group
    for readtype in readtypes:
        readtype = f'nreads_{readtype}_norm'
        # Create a color palette for the heatmap
        cmap = sns.diverging_palette(255, 85, s=255, l=70, sep=20, as_cmap=True)
        # Create a correlation matrix from reads stored in adata observations
        df, log2fc_dict = toolsTG.adataLog2FC(adata, grp, readtype, readcount_cutoff=cutoff, config_name=config_name, overwrite=overwrite).main()
        df['readtype'] = readtype
        # combine df with df_combine by stacking them vertically if readtype is not total_unique or total
        if readtype != 'nreads_total_unique_norm' and readtype != 'nreads_total_norm':
            df_combine = pd.concat([df_combine, df], axis=0)
        # save df to csv
        df.to_csv(f'{output}{grp}_{readtype}_{cutoff}_{heatbound}_heatmap.csv')
        # Create a pdf with a heatmap for sorted by each group on each page
        with PdfPages(f'{output}{grp}_{readtype}_{cutoff}_{heatbound}_heatmap.pdf') as pdf:
            for col in [i for i in df.columns.tolist() if 'log2' in i]:
                plt = heatmap_plot(df, col, cmap, heatbound)
                # Save figure
                if threaded:
                    threaded += f'Saving heatmap for {readtype} {col}...\n'
                else:
                    print(f'Saving heatmap for {readtype} {col}...')
                if heatsubplots:
                    plt.savefig(f'{output}{grp}_{readtype}_{cutoff}_{heatbound}_{col}_heatmap.pdf', bbox_inches='tight')
                pdf.savefig(bbox_inches='tight')
                plt.close()
    # Create a heatmap for the combined groups
    if not df_combine.empty:
        with PdfPages(f'{output}{grp}_combine_{cutoff}_{heatbound}_heatmap.pdf') as pdf:
            for col in [i for i in df_combine.columns.tolist() if 'log2' in i]:
                plt = heatmap_plot(df_combine, col, cmap, heatbound)
                # Save figure
                if threaded:
                    threaded += f'Saving heatmap for combine {col}...\n'
                else:
                    print(f'Saving heatmap for combine {col}...')
                if heatsubplots:
                    plt.savefig(f'{output}{grp}_combine_{cutoff}_{heatbound}_{col}_heatmap.pdf', bbox_inches='tight')
                pdf.savefig(bbox_inches='tight')
                plt.close()
    if threaded:
        return threaded

def heatmap_plot(df, col, cmap, heatbound):
    # Sort df by the sum of the log2FC values for each row
    tdf = df.sort_values(by=col, ascending=False)
    # subset df to only include the top and bottom heatbound values only if the heatmap is larger than heatbound
    if len(tdf) > int(heatbound):
        tdf = pd.concat([tdf.iloc[:int(heatbound), :], tdf.iloc[-int(heatbound):, :]])  
    # Create a heatmap
    if len(tdf.columns)>=20:
        fig, axs = plt.subplots(1, 2, figsize=(16,12))
    else:
        fig, axs = plt.subplots(1, 2, figsize=(6,12))

    log_tdf = tdf[[i for i in tdf.columns if 'log2' in i]]
    pval_tdf = tdf[[i for i in tdf.columns if 'pval' in i]]
    sns.heatmap(log_tdf, ax=axs[0], cmap=cmap, center=0, vmax=4, vmin=-4, cbar=True, square=True, cbar_kws={'fraction':0.05, 'pad':0.05})
    sns.heatmap(-np.log10(pval_tdf), ax=axs[1], cmap='Greens', vmin=0, vmax=3, cbar=True, yticklabels=False, square=True, cbar_kws={'fraction':0.05, 'pad':0.05})
    axs[0].tick_params(axis='x', labelrotation=90)
    axs[1].tick_params(axis='x', labelrotation=90)
    axs[1].set_ylabel('')

    # # Set a column of symbols to the left of the heatmap based on readtype
    # for i, row in enumerate(tdf.index.tolist()):
    #     if tdf.loc[row, 'readtype'] == 'nreads_fiveprime_unique_norm' or tdf.loc[row, 'readtype'] == 'nreads_fiveprime_norm':
    #         axs[0].text(-0.5, i+0.5, '\u25CF', fontsize=8, horizontalalignment='center', verticalalignment='center')
    #     elif tdf.loc[row, 'readtype'] == 'nreads_threeprime_unique_norm' or tdf.loc[row, 'readtype'] == 'nreads_threeprime_norm':
    #         axs[0].text(-0.5, i+0.5, '\u25A0', fontsize=8, horizontalalignment='center', verticalalignment='center')
    #     elif tdf.loc[row, 'readtype'] == 'nreads_other_unique_norm' or tdf.loc[row, 'readtype'] == 'nreads_other_norm':
    #         axs[0].text(-0.5, i+0.5, '\u25B2', fontsize=8, horizontalalignment='center', verticalalignment='center')
    #     elif tdf.loc[row, 'readtype'] == 'nreads_whole_unique_norm' or tdf.loc[row, 'readtype'] == 'nreads_whole_norm':
    #         axs[0].text(-0.5, i+0.5, '\u25C6', fontsize=8, horizontalalignment='center', verticalalignment='center')
    #     elif tdf.loc[row, 'readtype'] == 'nreads_precounts_norm':
    #         axs[0].text(-0.5, i+0.5, '\u2715', fontsize=8, horizontalalignment='center', verticalalignment='center')
    # # Adjust xticklabels to the left of the heatmap to make room for the symbols
    # axs[0].set_xticklabels(['']+tdf.columns.tolist(), fontsize=8)

    cbar = axs[0].collections[0].colorbar
    # cbar.set_ticks([-2, -1, 0, 1, 2])
    # cbar.set_ticklabels(['<=-2', '-1', '0', '1', '>=2'])
    cbar.set_ticks([-4,-3,-2,-1,0,1,2,3,4])
    cbar.set_ticklabels(['<=-4','-3','-2','-1','0','1','2','3','>=4'])
    cbar.ax.tick_params(labelsize=8) 
    cbar = axs[1].collections[0].colorbar
    cbar.set_ticks([0, 1.3, 3])
    cbar.set_ticklabels(['1', '0.05', '<=0.001'])
    cbar.ax.tick_params(labelsize=8)

    # Add a title to the figure
    plt.suptitle(f'Heatmap of log2FC of tRNA read counts between groups \nsorted by {col}', fontsize=12, y=1.15)
    axs[0].set_title('log2FC', fontsize=10)
    axs[1].set_title('pval', fontsize=10)

    return plt

if __name__ == '__main__':
    pass
