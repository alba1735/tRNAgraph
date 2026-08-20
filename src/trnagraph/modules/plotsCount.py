#!/usr/bin/env python3

import logging

import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
import pandas as pd
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

logger = logging.getLogger(__name__)

def visualizer(adata, colormap, output, threaded=True):
    '''
    Generate barplots for readtypes and isoforms.
    '''
    # Load data from AnnData
    for count_type in ['type_counts', 'amino_counts']:
        df = adata.uns[count_type].copy()
        
        # We need mapping from sample to group
        grp_dict = {k:v for k,v in zip(adata.obs['sample'], adata.obs['group'])}
        
        is_sample_level = any(col in grp_dict for col in df.columns)
        
        c_type_amino = colormap.get(count_type.split('_')[0], {}) if colormap else {}
        c_group = colormap.get('group', {}) if colormap else {}
        
        # Build c_sample if colormap has 'group' and is not explicitly provided
        c_sample = colormap.get('sample', {}) if colormap else {}
        if not c_sample and c_group:
            for sample, group in grp_dict.items():
                if group in c_group:
                    c_sample[sample] = c_group[group]
        
        # Convert hex correctly
        c_type_amino = {k:v if v[0]!='#' else mplcolors.to_rgb(v) for k,v in c_type_amino.items()}
        c_group = {k:v if v[0]!='#' else mplcolors.to_rgb(v) for k,v in c_group.items()}
        c_sample = {k:v if v[0]!='#' else mplcolors.to_rgb(v) for k,v in c_sample.items()}
        
        if is_sample_level:
            # Group mean DataFrames:
            # Columns are groups instead of samples
            df_group_mean = df.T.groupby(grp_dict, observed=False).mean().T
            df_group_percent = df_group_mean*100/df_group_mean.sum()
            
            # sample-level
            threaded = stacked_barplots(df.copy(), count_type, 'sample', output, colormap=c_type_amino, threaded=threaded)
            threaded = stacked_barplots(df*100/df.sum(), count_type, 'sample', output, percent=True, colormap=c_type_amino, threaded=threaded)
            
            threaded = split_barplots(df.copy(), count_type, 'sample', output, colormap=c_sample, threaded=threaded)
            threaded = split_barplots(df*100/df.sum(), count_type, 'sample', output, percent=True, colormap=c_sample, threaded=threaded)
        else:
            df_group_mean = df.copy()
            df_group_percent = df_group_mean*100/df_group_mean.sum()
        
        # group-level (using mean counts for total)
        threaded = stacked_barplots(df_group_mean.copy(), count_type, 'group', output, colormap=c_type_amino, threaded=threaded)
        threaded = stacked_barplots(df_group_percent.copy(), count_type, 'group', output, percent=True, colormap=c_type_amino, threaded=threaded)

        threaded = split_barplots(df_group_mean.copy(), count_type, 'group', output, colormap=c_group, threaded=threaded)
        threaded = split_barplots(df_group_percent.copy(), count_type, 'group', output, percent=True, colormap=c_group, threaded=threaded)

    if threaded:
        return threaded

def split_barplots(df, count_type, level, output, colormap=None, percent=False, threaded=True):
    '''
    Create split barplots for readtypes.
    X-axis is type (amino/type), Hue is group (which is columns from df, i.e., sample or group).
    '''
    df['type'] = df.index
    fig, ax = plt.subplots(figsize=(12, 8))
    df = df.melt(id_vars='type', var_name='group', value_name='count')
    
    if colormap:
        sns.barplot(x='type', y='count', hue='group', errorbar=None, palette=colormap, data=df, ax=ax)
    else:
        sns.barplot(x='type', y='count', hue='group', errorbar=None, palette=sns.husl_palette(len(df['group'].unique())), data=df, ax=ax)
        
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels, loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
    ax.legend_.set_title(level.capitalize())
    
    # Apply plotsBar styling
    sns.despine(left=True, bottom=True)
    
    if percent:
        ax.set_ylabel('Percentage of Total Reads')
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%'])
        ax.set_title(f'Percentage of Total Reads by {level.capitalize()}')
    else:
        ax.set_ylabel('Readcounts')
        ax.set_title(f'Total Readcounts by {level.capitalize()}')
        
    ax.set_xlabel(f"{count_type.split('_')[0].capitalize()}")
    plt.xticks(rotation=90)
    
    prefix = 'percent' if percent else 'totalcounts'
    c_type = count_type.split('_')[0]
    filename = f"{output}{level}_{prefix}_{c_type}_split.pdf"
    plt.savefig(filename, bbox_inches='tight')
    plt.close()
    
    ps = f'Plot saved to {filename}'
    if threaded:
        threaded += f'{ps}\n'
        return threaded
    else:
        logger.info(ps)
        return None

def stacked_barplots(df, count_type, level, output, colormap=None, percent=False, threaded=True):
    '''
    Create stacked barplots for readtypes.
    X-axis is sample/group (columns), Hue is amino/type (index).
    '''
    fig, ax = plt.subplots(figsize=(12, 8))
    bar_bottom = len(df.columns) * [0]
    
    if colormap:
        pal = [colormap.get(bar, sns.husl_palette(len(df.index))[i]) for i, bar in enumerate(df.index)]
    else:
        pal = sns.husl_palette(len(df.index))
        
    pal_position=0
    for bar in df.index.values[::-1]:
        # apply plotsBar stacking styling: linewidth=0.5
        ax.bar(df.columns, df.loc[bar], 0.9, bottom=bar_bottom, color=pal[pal_position], label=bar, linewidth=0.5, edgecolor='black', clip_on=False)
        bar_bottom += df.loc[bar]
        pal_position += 1
        
    sns.despine(left=True, bottom=True)
    
    if percent:
        ax.set_ylabel('Percentage of Total Reads')
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%'])
        ax.set_title(f'Percentage of Total Reads of {level.capitalize()}')
    else:
        ax.set_ylabel('Readcounts')
        ax.set_title(f'Total {count_type.replace("_", " ").capitalize()}')
        
    ax.set_xlabel('Group' if level == 'group' else 'Sample')
    plt.xticks(rotation=90)
    
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(list(reversed(handles)), list(reversed(labels)), loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
    ax.legend_.set_title(f"{count_type.split('_')[0].capitalize()} Group")
    
    prefix = 'percent' if percent else 'totalcounts'
    c_type = count_type.split('_')[0]
    filename = f"{output}{level}_{prefix}_{c_type}_stacked.pdf"
    
    plt.savefig(filename, bbox_inches='tight')
    plt.close()
    
    ps = f'Plot saved to {filename}'
    if threaded:
        threaded += f'{ps}\n'
        return threaded
    else:
        logger.info(ps)
        return None

if __name__ == '__main__':
    pass