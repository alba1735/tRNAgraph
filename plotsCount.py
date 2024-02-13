#!/usr/bin/env python3

import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

def visualizer(adata, colormap_tc, colormap_bg, output, threaded=True):
    '''
    Generate barplots for readtypes and isoforms.
    '''
    grp = 'sample' # group by sample by default this is a count plot and thus static
    
    # Get the unique grp as a dictionary with sample names as keys
    # grp_dict = {k:v for k,v in zip(adata.obs['sample'], adata.obs[grp])}

    # Load data from AnnData
    for count_type in ['type_counts', 'amino_counts']:
        use_colormap = False
        df = adata.uns[count_type].copy()
        # create a combine df if count_type is amino_counts
        if count_type == 'amino_counts':
            df_combine = df.copy()
            # merge columns using the timepoint_dict
            # df_combine_mean = df_combine.groupby(grp_dict, axis=1).mean()
            # df_combine_sum = df_combine.groupby(grp_dict, axis=1).sum()
            # The above is deprecated because it doesn't work with the new version of pandas use DataFrame.groupby with axis=1 is deprecated. Do `frame.T.groupby(...)` without axis instead.
            df_combine_mean = df_combine.T.groupby(level=0, observed=False).mean().T
            df_combine_sum = df_combine.T.groupby(level=0, observed=False).sum().T

        if count_type == 'type_counts':
            colormap = colormap_tc
        else:
            colormap = colormap_bg

        if colormap != None:
            colormap = {k:v if v[0]!='#' else mplcolors.to_rgb(v) for k,v in colormap.items()}
            use_colormap = True
            for v in df.columns.unique():
                if v not in colormap:
                    if threaded:
                        threaded += f'Color {v} not found in colormap. Using default colors instead.\n'
                    else:
                        print(f'Color {v} not found in colormap. Using default colors instead.')
                    use_colormap = False
                    break

        if use_colormap:
            split_barplots(df.copy(), count_type, output, colormap=colormap, threaded=threaded)
            split_barplots(df*100/df.sum(), count_type, output, colormap=colormap, percent=True, threaded=threaded)
            if count_type == 'amino_counts' and grp != 'sample':
                split_barplots(df_combine_mean.copy(), count_type, output, title='mean', colormap=colormap, threaded=threaded)
                split_barplots(df_combine_mean*100/df_combine_mean.sum(), count_type, output, title='mean', colormap=colormap, percent=True, threaded=threaded)
                split_barplots(df_combine_sum.copy(), count_type, output, title='sum', colormap=colormap, threaded=threaded)
                split_barplots(df_combine_sum*100/df_combine_sum.sum(), count_type, output, title='sum', colormap=colormap, percent=True, threaded=threaded)
        else:
            split_barplots(df.copy(), count_type, output, threaded=threaded)
            split_barplots(df*100/df.sum(), count_type, output, percent=True, threaded=threaded)
            if count_type == 'amino_counts' and grp != 'sample':
                split_barplots(df_combine_mean.copy(), count_type, output, title='mean', threaded=threaded)
                split_barplots(df_combine_mean*100/df_combine_mean.sum(), count_type, output, title='mean', percent=True, threaded=threaded)
                split_barplots(df_combine_sum.copy(), count_type, output, title='sum', threaded=threaded)
                split_barplots(df_combine_sum*100/df_combine_sum.sum(), count_type, output, title='sum', percent=True, threaded=threaded)
        
        stacked_barplots(df.copy(), count_type, output, threaded=threaded)
        stacked_barplots(df*100/df.sum(), count_type, output, percent=True, threaded=threaded)
        if count_type == 'amino_counts' and grp != 'sample':
            stacked_barplots(df_combine_mean.copy(), count_type, output, title='mean', threaded=threaded)
            stacked_barplots(df_combine_mean*100/df_combine_mean.sum(), count_type, output, title='mean', percent=True, threaded=threaded)
            stacked_barplots(df_combine_sum.copy(), count_type, output, title='sum', threaded=threaded)
            stacked_barplots(df_combine_sum*100/df_combine_sum.sum(), count_type, output, title='sum', percent=True, threaded=threaded)

    if threaded:
        return threaded

def split_barplots(df, count_type, output, title=None, colormap=None, percent=False, threaded=True):
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
        if title:
            plt.savefig(f'{output}percent_{count_type}_split_{title}.pdf', bbox_inches='tight')
            print(f'Plot saved to {output}percent_{count_type}_split_{title}.pdf')
        else:
            plt.savefig(f'{output}percent_{count_type}_split.pdf', bbox_inches='tight')
            print(f'Plot saved to {output}percent_{count_type}_split.pdf')
    else:
        if title:
            plt.savefig(f'{output}total_{count_type}_split_{title}.pdf', bbox_inches='tight')
            print(f'Plot saved to {output}total_{count_type}_split_{title}.pdf')
        else:
            plt.savefig(f'{output}total_{count_type}_split.pdf', bbox_inches='tight')
            print(f'Plot saved to {output}total_{count_type}_split.pdf')
    plt.close()

def stacked_barplots(df, count_type, output, title=None, percent=False, threaded=True):
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
        if title:
            plt.savefig(f'{output}percent_{count_type}_stacked_{title}.pdf', bbox_inches='tight')
            print(f'Plot saved to {output}percent_{count_type}_stacked_{title}.pdf')
        else:
            plt.savefig(f'{output}percent_{count_type}_stacked.pdf', bbox_inches='tight')
            print(f'Plot saved to {output}percent_{count_type}_stacked.pdf')
    else:
        if title:
            plt.savefig(f'{output}total_{count_type}_stacked_{title}.pdf', bbox_inches='tight')
            print(f'Plot saved to {output}total_{count_type}_stacked_{title}.pdf')
        else:
            plt.savefig(f'{output}total_{count_type}_stacked.pdf', bbox_inches='tight')
            print(f'Plot saved to {output}total_{count_type}_stacked.pdf')
    plt.close()

if __name__ == '__main__':
    pass