#!/usr/bin/env python3

import logging

import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
import pandas as pd

from . import toolsTG
from . import plotsPalette
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

logger = logging.getLogger(__name__)

# The collapsed band's label and color. Deliberately "non-tRNA" rather than "small RNA":
# these are whatever the GTF carried (rRNA, snoRNA, protein_coding, Unknown, ...), which is
# not a small-RNA set. Gray so it reads as background against the tRNA categories, which are
# the subject of the plot.
NONTRNA_COLLAPSED_LABEL = 'non-tRNA'
NONTRNA_COLLAPSED_COLOR = plotsPalette.NONTRNA_COLLAPSED


def collapse_nontrna_types(df):
    '''
    Fold every non-tRNA row of a type-count frame into one summed `non-tRNA` row.

    The boundary is adataBuild.TRNA_TYPE_LABELS -- the same definition the split-variant
    filters use, which keeps Mt_tRNA on the tRNA side. tRNA rows keep their original order
    and the collapsed band is appended, so the tRNA part of the stack looks unchanged.

    This is a regrouping, not a filter: the column totals are identical afterwards.
    '''
    from .adataBuild import TRNA_TYPE_LABELS

    if df.empty:
        return df
    is_trna = [label in TRNA_TYPE_LABELS for label in df.index]
    trna_part = df[is_trna]
    nontrna_part = df[[not flag for flag in is_trna]]
    if nontrna_part.empty:
        return df

    collapsed = nontrna_part.sum().rename(NONTRNA_COLLAPSED_LABEL).to_frame().T
    return pd.concat([trna_part, collapsed]) if not trna_part.empty else collapsed

def visualizer(adata, colormap, output, threaded=True, settings=None):
    '''
    Generate barplots for readtypes and isoforms.
    '''
    # Load data from AnnData
    for count_type in ['type_counts', 'amino_counts']:
        df = adata.uns[count_type].copy()
        # Only the gene-type breakdown has non-tRNA rows to fold together; amino_counts is
        # tRNA-only by construction, so a collapsed variant of it would be the same plot.
        collapsible = count_type == 'type_counts' and not collapse_nontrna_types(df).equals(df)
        
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
            threaded = stacked_barplots(df.copy(), count_type, 'sample', output, colormap=c_type_amino, threaded=threaded, settings=settings)
            threaded = stacked_barplots(df*100/df.sum(), count_type, 'sample', output, percent=True, colormap=c_type_amino, threaded=threaded, settings=settings)
            if collapsible:
                cdf = collapse_nontrna_types(df)
                threaded = stacked_barplots(cdf.copy(), count_type, 'sample', output, colormap=c_type_amino, threaded=threaded, settings=settings, collapsed=True)
                threaded = stacked_barplots(cdf*100/cdf.sum(), count_type, 'sample', output, percent=True, colormap=c_type_amino, threaded=threaded, settings=settings, collapsed=True)
            
            threaded = split_barplots(df.copy(), count_type, 'sample', output, colormap=c_sample, threaded=threaded, settings=settings)
            threaded = split_barplots(df*100/df.sum(), count_type, 'sample', output, percent=True, colormap=c_sample, threaded=threaded, settings=settings)
        else:
            df_group_mean = df.copy()
            df_group_percent = df_group_mean*100/df_group_mean.sum()
        
        # group-level (using mean counts for total)
        threaded = stacked_barplots(df_group_mean.copy(), count_type, 'group', output, colormap=c_type_amino, threaded=threaded, settings=settings)
        threaded = stacked_barplots(df_group_percent.copy(), count_type, 'group', output, percent=True, colormap=c_type_amino, threaded=threaded, settings=settings)
        if collapsible:
            cgm = collapse_nontrna_types(df_group_mean)
            threaded = stacked_barplots(cgm.copy(), count_type, 'group', output, colormap=c_type_amino, threaded=threaded, settings=settings, collapsed=True)
            threaded = stacked_barplots(cgm*100/cgm.sum(), count_type, 'group', output, percent=True, colormap=c_type_amino, threaded=threaded, settings=settings, collapsed=True)

        threaded = split_barplots(df_group_mean.copy(), count_type, 'group', output, colormap=c_group, threaded=threaded, settings=settings)
        threaded = split_barplots(df_group_percent.copy(), count_type, 'group', output, percent=True, colormap=c_group, threaded=threaded, settings=settings)

    if threaded:
        return threaded

def split_barplots(df, count_type, level, output, colormap=None, percent=False, threaded=True, settings=None):
    '''
    Create split barplots for readtypes.
    X-axis is type (amino/type), Hue is group (which is columns from df, i.e., sample or group).
    '''
    df['type'] = df.index
    fig, ax = plt.subplots(figsize=toolsTG.figsize_for(settings, (12, 8)))
    df = df.melt(id_vars='type', var_name='group', value_name='count')
    
    if colormap:
        sns.barplot(x='type', y='count', hue='group', errorbar=None, palette=colormap, data=df, ax=ax)
    else:
        sns.barplot(x='type', y='count', hue='group', errorbar=None, palette=plotsPalette.categorical(settings, len(df['group'].unique())), data=df, ax=ax)
        
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
    toolsTG.save_current(filename, settings)
    plt.close()
    
    ps = f'Plot saved to {filename}'
    if threaded:
        threaded += f'{ps}\n'
        return threaded
    else:
        logger.info(ps)
        return None

def stacked_barplots(df, count_type, level, output, colormap=None, percent=False, threaded=True, settings=None, collapsed=False):
    '''
    Create stacked barplots for readtypes.
    X-axis is sample/group (columns), Hue is amino/type (index).
    '''
    fig, ax = plt.subplots(figsize=toolsTG.figsize_for(settings, (12, 8)))
    bar_bottom = len(df.columns) * [0]
    
    # Keyed by label, not by position. The bars are drawn in reverse index order while the
    # palette used to be built in forward order, so a colormap entry landed on whichever bar
    # sat at the mirrored position instead of on the bar it named.
    fallback = plotsPalette.categorical(settings, len(df.index))
    bar_colors = {bar: (colormap or {}).get(bar, fallback[i]) for i, bar in enumerate(df.index)}
    if collapsed:
        # The collapsed band is background, not a category competing with the tRNA types, so
        # it takes a fixed gray unless the user's own colormap names it explicitly.
        bar_colors[NONTRNA_COLLAPSED_LABEL] = (colormap or {}).get(
            NONTRNA_COLLAPSED_LABEL, NONTRNA_COLLAPSED_COLOR)
        
    for bar in df.index.values[::-1]:
        # apply plotsBar stacking styling: linewidth=0.5
        ax.bar(df.columns, df.loc[bar], 0.9, bottom=bar_bottom, color=bar_colors[bar], label=bar,
               linewidth=toolsTG.linewidth_for(settings, 0.5), edgecolor='black', clip_on=False)
        bar_bottom += df.loc[bar]
        
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
    # The path names what was plotted, not which flag was typed.
    suffix = '_collapsed' if collapsed else ''
    filename = f"{output}{level}_{prefix}_{c_type}{suffix}_stacked.pdf"
    
    toolsTG.save_current(filename, settings)
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