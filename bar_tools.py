#!/usr/bin/env python3

import seaborn as sns
import anndata as ad
import argparse

import directory_tools

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


def visualizer(adata, colormap, output):
    '''
    Generate barplots for readtypes and isoforms.
    '''
    # Dictionary for readtypes to combine misc into other
    # type_dict = {'tRNA':'tRNA', 'pretRNA':'pretRNA', 'snoRNA':'snoRNA', 'snRNA':'snRNA', 'scaRNA':'scaRNA', 'sRNA':'other', 'miRNA':'miRNA', 
    #              'misc_RNA':'other', 'rRNA_pseudogene':'rRNA_pseudogene', 'ribozyme':'other', 'lncRNA':'other', 'protein_coding':'other', 
    #              'Mt_rRNA':'other', 'Mt_tRNA':'Mt_tRNA', 'rRNA':'rRNA', 'other':'other', 'Y_RNA':'other', 'snoRNA_pseudogene':'snoRNA_pseudogene',
    #              'vault_RNA':'vault_RNA', 'snRNA_pseudogene':'snRNA_pseudogene', 'misc_RNA_pseudogene':'misc_RNA_pseudogene', 'pretRNA_amtosemse':'pretRNA_amtosemse',
    #              'tRNA_antisense':'tRNA_antisense'}
    # # Add column for readtypes and convert to nice names and categories from dictionary
    # df['type'] = [type_dict.get(i) for i in df.index]

    # Load data from AnnData
    for count_type in ['type_counts', 'amino_counts']:
        df = adata.uns[count_type]

        if colormap != None:
            colormap = {k:v if v[0]!='#' else mplcolors.to_rgb(v) for k,v in colormap.items()}
            use_colormap = True
            for v in df.columns.unique():
                if v not in colormap:
                    print(f'Color {v} not found in colormap. Using default colors instead.')
                    use_colormap = False
                    break

        if use_colormap:
            print('cmap')
            split_barplots(df.copy(), count_type, output, colormap=colormap)
            split_barplots(df*100/df.sum(), count_type, output, colormap=colormap, percent=True)
        else:
            split_barplots(df.copy(), count_type, output)
            split_barplots(df*100/df.sum(), count_type, output, percent=True)
        
        stacked_barplots(df.copy(), count_type, output)
        stacked_barplots(df*100/df.sum(), count_type, output, percent=True)

def split_barplots(df, count_type, output, colormap=None, percent=False):
    '''
    Create split barplots for readtypes.
    '''
    df['type'] = df.index
    # Create split barplots
    fig, ax = plt.subplots(figsize=(10, 10))
    # Melt the dataframe so it can be used for plotting
    df = df.melt(id_vars='type', var_name='group', value_name='count')
    # Create barplot for readtypes
    if colormap != None:
        sns.barplot(x='type', y='count', hue='group', errorbar=None, palette=colormap, data=df, ax=ax)
    else:
        sns.barplot(x='type', y='count', hue='group', errorbar=None, palette=sns.husl_palette(len(df['group'].unique())), data=df, ax=ax)
    # Legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels, loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
    ax.legend_.set_title('Groups')
    # Labels
    if percent:
        ax.set_ylabel('Percentage of Total Reads')
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%'])
        ax.set_title('Percentage of Total Reads')
    else:
        ax.set_ylabel('Readcounts')
        ax.set_title('Total Readcounts')
    ax.set_xlabel(f'{count_type} Group')
    # rotate x labels
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    # Set the box aspect ratio to 1 so the plot is square
    plt.gca().set_box_aspect(1)
    # Save
    if percent:
        plt.savefig(f'{output}/percent_{count_type}_split.pdf', bbox_inches='tight')
        print(f'Plot saved to {output}/percent_{count_type}_split.pdf')
    else:
        plt.savefig(f'{output}/total_{count_type}_split.pdf', bbox_inches='tight')
        print(f'Plot saved to {output}/total_{count_type}_split.pdf')
    plt.close()

def stacked_barplots(df, count_type, output, percent=False):
    '''
    Create stacked barplots for readtypes.
    '''
    # Create stacked barplots
    fig, ax = plt.subplots(figsize=(10, 10))
    # Create parameters for barplot
    bar_bottom = len(df.columns) * [0]
    # Create palette for barplot
    pal = sns.husl_palette(len(df.index))
    pal_position=0
    # Create barplot for readtypes
    for bar in df.index.values[::-1]:
        ax.bar(df.columns, df.loc[bar], 0.9, bottom=bar_bottom, color=pal[pal_position], label=bar, linewidth=1, edgecolor='black', clip_on=False)
        bar_bottom += df.loc[bar]
        pal_position += 1
    # Spines
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    # Labels
    if percent:
        ax.set_ylabel('Percentage of Total Reads')
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%'])
        ax.set_title('Percentage of Total Reads')
    else:
        ax.set_ylabel('Readcounts')
        ax.set_title(f'Total {count_type}')
    ax.set_xlabel('Group')
    # rotate x labels
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    # Legend - Reverse order to match barplot
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(reversed(handles), reversed(labels), loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
    ax.legend_.set_title(f'{count_type} Group')
    # Set the box aspect ratio to 1 so the plot is square
    plt.gca().set_box_aspect(1)
    # Save
    if percent:
        plt.savefig(f'{output}/percent_{count_type}_stacked.pdf', bbox_inches='tight')
        print(f'Plot saved to {output}/percent_{count_type}_stacked.pdf')
    else:
        plt.savefig(f'{output}/total_{count_type}_stacked.pdf', bbox_inches='tight')
        print(f'Plot saved to {output}/total_{count_type}_stacked.pdf')
    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='bar_tools.py',
        description='Generate barplots for readtypes and isoforms.',
    )

    parser.add_argument('-i', '--anndata', help='Specify AnnData input', required=True)
    parser.add_argument('-o', '--output', help='Specify output directory', default='barplot', required=False)
    parser.add_argument('--colormap', help='Specify a colormap for coverage plots (optional)', default=None)

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    directory_tools.builder(args.output)

    adata = ad.read_h5ad(args.anndata)

    visualizer(adata, args.colormap, args.output)