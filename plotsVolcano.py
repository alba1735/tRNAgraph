#!/usr/bin/env python3

import seaborn as sns
import numpy as np

import toolsTG

import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


def visualizer(adata, grp, readtype, cutoff, output, threaded=True):
    '''
    Generate volcano visualizations for each group in an AnnData object.
    '''
    # Create a correlation matrix from reads stored in adata observations
    df = toolsTG.adataLog2FC(adata, grp, readtype, overwrite=False, readcount_cutoff=cutoff).main()

    pairs = ['_'.join(i.split('_')[1:]) for i in df.columns]
    for pair in pairs:
        plt.figure(figsize=(8, 8))
        x = df[f'log2_{pair}']
        y = -np.log10(df[f'pval_{pair}'])
        ax = sns.scatterplot(x=x, y=y)
        # Add a line at log2FC = 1.5 and -1.5 and pval = 0.05 and 0.001
        ax.axvline(x=-1.5, color='black', linestyle='--')
        ax.axvline(x=1.5, color='black', linestyle='--')
        ax.axhline(y=1.3, color='black', linestyle='--')
        ax.axhline(y=3, color='black', linestyle='--')
        # Set axis limits so that the plot is square
        if ax.get_xlim()[1] < 3 and ax.get_xlim()[0] > -3:
            ax.set_xlim(-3, 3)
        else:
            ax.set_xlim(-1.1*max(abs(x)), 1.1*max(abs(x)))
        if ax.get_ylim()[1] < 10:
            ax.set_ylim(0, 10)
        # Add labels
        ax.set_xlabel('log2(fold change)')
        ax.set_ylabel('-log10(p-value)')
        ax.set_title(f'Volcano plot of {pair} tRNA read counts')
        # Save figure
        if threaded:
            threaded += f'Saving figure: {output}{pair}_{readtype}_{cutoff}_volcano.pdf\n'
        else:   
            print(f'Saving figure: {output}{pair}_{readtype}_{cutoff}_volcano.pdf')
        plt.savefig(f'{output}{pair}_{readtype}_{cutoff}_volcano.pdf', bbox_inches='tight')
        plt.close()

    if threaded:
        return threaded

if __name__ == '__main__':
    pass