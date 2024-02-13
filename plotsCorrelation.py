#!/usr/bin/env python3

import pandas as pd
import anndata as ad

import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns

def visualizer(adata, corr_method, corr_group, output, threaded=True):
    '''
    Generate correlation graphs for each sample in an AnnData object.
    '''

    # print(adata.obs.columns) # add sample group to adata.obs.columns from samples file so this can be used for a secondary set of plots
    # Add shorter sample names flag to allow for cleaner plots

    # Create a correlation matrix from reads stored in adata observations
    df = pd.DataFrame(adata.obs, columns=['trna', corr_group] + [i for i in adata.obs.columns if '_norm' in i])

    for i in df.columns[2:]:
        df_corr = df.pivot_table(index='trna', columns=corr_group, values=i, observed=True)
        # Only plot correlation matrices with more than 20 samples will be generated
        if df_corr.max().max() < 20:
            if threaded:
                threaded += f'Not enough samples to generate correlation matrix for {i}\n'
            else:
                print(f'Not enough samples to generate correlation matrix for {i}')
        else:
            if threaded:
                threaded += f'Generating correlation matrix for {i}\n'
            else:
                print(f'Generating correlation matrix for {i}')
            df_corr = df_corr.corr(method=corr_method)
            # Plot the correlation matrix
            plt.figure(figsize=(6, 6))
            ax = sns.heatmap(df_corr**2, square=True, cmap='Blues', cbar_kws={'label': f'{corr_method} R^2'})
            # Remove the axis labels and set the title
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_title(f'{corr_method} {corr_group} {i.split("_")[1]} Correlation Matrix'.title())
            # Set the box aspect ratio to 1 so the plot is square
            plt.gca().set_box_aspect(1)
            # Save the plot
            plt.savefig(f'{output}{corr_method}_{corr_group}_{i.split("_")[1]}_correlation_matrix.pdf', bbox_inches='tight')
            if threaded:
                threaded += f'Correlation matrix for {i} saved to {output}{corr_method}_{corr_group}_{i.split("_")[1]}_correlation_matrix.pdf\n'
            else:
                print(f'Correlation matrix for {i} saved to {output}{corr_method}_{corr_group}_{i.split("_")[1]}_correlation_matrix.pdf')
            plt.close()

    if threaded:
        return threaded


if __name__ == '__main__':
    pass