#!/usr/bin/env python3

# import pandas as pd
import seaborn as sns
# import anndata as ad
# import argparse

# import toolsTG

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


def visualizer(adata, colormap_tc, colormap_bg, output):
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
                    print(f'Color {v} not found in colormap. Using default colors instead.')
                    use_colormap = False
                    break

        if use_colormap:
            split_barplots(df.copy(), count_type, output, colormap=colormap)
            split_barplots(df*100/df.sum(), count_type, output, colormap=colormap, percent=True)
            if count_type == 'amino_counts' and grp != 'sample':
                split_barplots(df_combine_mean.copy(), count_type, output, title='mean', colormap=colormap)
                split_barplots(df_combine_mean*100/df_combine_mean.sum(), count_type, output, title='mean', colormap=colormap, percent=True)
                split_barplots(df_combine_sum.copy(), count_type, output, title='sum', colormap=colormap)
                split_barplots(df_combine_sum*100/df_combine_sum.sum(), count_type, output, title='sum', colormap=colormap, percent=True)
        else:
            split_barplots(df.copy(), count_type, output)
            split_barplots(df*100/df.sum(), count_type, output, percent=True)
            if count_type == 'amino_counts' and grp != 'sample':
                split_barplots(df_combine_mean.copy(), count_type, output, title='mean')
                split_barplots(df_combine_mean*100/df_combine_mean.sum(), count_type, output, title='mean', percent=True)
                split_barplots(df_combine_sum.copy(), count_type, output, title='sum')
                split_barplots(df_combine_sum*100/df_combine_sum.sum(), count_type, output, title='sum', percent=True)
        
        stacked_barplots(df.copy(), count_type, output)
        stacked_barplots(df*100/df.sum(), count_type, output, percent=True)
        if count_type == 'amino_counts' and grp != 'sample':
            stacked_barplots(df_combine_mean.copy(), count_type, output, title='mean')
            stacked_barplots(df_combine_mean*100/df_combine_mean.sum(), count_type, output, title='mean', percent=True)
            stacked_barplots(df_combine_sum.copy(), count_type, output, title='sum')
            stacked_barplots(df_combine_sum*100/df_combine_sum.sum(), count_type, output, title='sum', percent=True)

def split_barplots(df, count_type, output, title=None, colormap=None, percent=False):
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

def stacked_barplots(df, count_type, output, title=None, percent=False):
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
    # parser = argparse.ArgumentParser(
    #     prog='bar_tools.py',
    #     description='Generate barplots for readtypes and isoforms.',
    # )

    # parser.add_argument('-i', '--anndata', help='Specify AnnData input', required=True)
    # parser.add_argument('-o', '--output', help='Specify output directory', default='barplot', required=False)
    # parser.add_argument('--bargrp', help='Specify AnnData column to group by (default: group) (optional)', default='group', required=False)
    # parser.add_argument('--colormap_tc', help='Specify a colormap for coverage plots for type counts (optional)', default=None)
    # parser.add_argument('--colormap_bg', help='Specify a colormap for coverage plots for bargrp (optional)', default=None)

    # args = parser.parse_args()

    # # Create output directory if it doesn't exist
    # toolsTG.builder(args.output)

    # adata = ad.read_h5ad(args.anndata)

    # visualizer(adata, args.bargrp, args.colormap_tc, args.colormap_bg, args.output)

# Incoporate the following

# from statannot import add_stat_annotation

# df = pd.read_csv('/home/jovyan/data/rnaseq/trax/c305/c305_all_sep/c305_all_sep-typecounts.txt', delimiter='\t')
# df = df.drop(columns=['HEK293_s13_pos_dm','HEK293_s14_pos_dm','HEK293_s15_pos_dm','HEK293_s16_pos_dm','HEK293_s17_pos_dm','HEK293_s18_pos_dm',
#                       'HEK293_s1_neg_dm','HEK293_s2_neg_dm','HEK293_s3_neg_dm','HEK293_s4_neg_dm','HEK293_s5_neg_dm','HEK293_s6_neg_dm',
#                       'HEK293_s1_neg_arm','HEK293_s3_neg_arm','HEK293_s5_pos_arm','HEK293_s7_pos_arm'])

# df = df.T

# df_meta = np.array([i.split('_') for i in df.index.values])
# df['time'] = df_meta[:,0]
# df['alkb'] = df_meta[:,2]
# df['seq'] = df_meta[:,3]
# df['alkb_seq'] = df['alkb'] + '_' + df['seq']
# df['unique_samples'] = df.index.values
# df['sample_grp'] = ['_'.join(df_meta[i,0:1]) + '_' + '_'.join(df_meta[i,2:]) for i in range(len(df_meta))]

# df_arm = df[df['seq'] == 'arm']

# fig, ax = plt.subplots(figsize=(6,3))

# odr = ['d0_neg_arm','d14_neg_arm','d35_neg_arm','d70_neg_arm','d0_pos_arm','d14_pos_arm','d35_pos_arm','d70_pos_arm']
# bps = [('d0_neg_arm','d14_neg_arm'),('d0_neg_arm','d35_neg_arm'),('d0_neg_arm','d70_neg_arm'),('d0_pos_arm','d14_pos_arm'),('d0_pos_arm','d35_pos_arm'),('d0_pos_arm','d70_pos_arm')]

# ax = sns.barplot(x="sample_grp", y="tRNA", hue="time", data=df_arm, order=odr)
# add_stat_annotation(ax, data=df_arm, x="sample_grp", y="tRNA", order=odr, box_pairs=bps, test='t-test_ind',
#                    text_format='star', loc='inside', verbose=2, linewidth=1)

# leg = plt.legend(loc='upper left', bbox_to_anchor=(1.03, 1))
# leg.get_frame().set_linewidth(0.0)

# # ax.set_yscale('log')
# # plt.ylim([10**3.5,10**7])

# ax.set_axisbelow(True)
# ax.set_title('tRNA Reads')
# ax.set_xlabel('AlkB Treatment')
# ax.set_ylabel('Normalized Readcounts')
# ax.set_xticklabels(['0 Neg','14 Neg','35 Neg','70 Neg','0 Pos','14 Pos','35 Pos','70 Pos'], fontsize=10, rotation=90)

# plt.savefig('trna_reads_boxplot_time.pdf', bbox_inches='tight')


## Also this 

# readtypes = ['tRNA', 'pretRNA','snoRNA', 'snRNA', 'scaRNA', 'sRNA', 'miRNA', 'misc_RNA', 'rRNA_pseudogene', 'ribozyme', 'lncRNA', 'Mt_rRNA', 'Mt_tRNA', 'rRNA', 'other']

# df_arm = df[df['seq'] == 'arm']
# df_melt = pd.melt(df_arm, id_vars=['unique_samples'], value_vars=readtypes)

# df_meta = np.array([i.split('_') for i in df_melt['unique_samples']])
# df_melt['time'] = df_melt['variable'] + '_' + df_meta[:,0]
# df_melt['alkb'] = df_melt['variable'] + '_' + df_meta[:,2]
# df_melt['seq'] = df_melt['variable'] + '_' + df_meta[:,3]
# df_melt['alkb_seq'] = df_melt['alkb'] + '_' + df_meta[:,3]
# df_melt['sample_grp'] = ['_'.join(df_meta[i,0:1]) + '_' + '_'.join(df_meta[i,2:]) for i in range(len(df_meta))]

# def bpplot(df,readtype):
#     df = df[df_melt['variable'] == readtype]

#     fig, ax = plt.subplots(figsize=(6,3))

#     bps = [('{}_neg'.format(readtype),'{}_pos'.format(readtype))]
    
#     pal = {'{}_d0'.format(readtype):'#007FFF','{}_d14'.format(readtype):'#00FF7F','{}_d35'.format(readtype):'#FF007F','{}_d70'.format(readtype):'#FFD700'}

#     ax = sns.barplot(x="alkb", y="value", hue="time", data=df, palette=pal)
#     add_stat_annotation(ax, data=df, x="alkb", y="value", box_pairs=bps, test='t-test_ind', text_format='star', loc='inside', verbose=2, linewidth=1)

#     leg = plt.legend(loc='upper left', bbox_to_anchor=(1.03, 1))
#     leg.get_frame().set_linewidth(0.0)

#     # ax.set_yscale('log')
#     # plt.ylim([10**3.5,10**7])

#     ax.set_axisbelow(True)
#     ax.set_title('{} Reads'.format(readtype))
#     ax.set_xlabel('Samples')
#     ax.set_ylabel('Normalized Readcounts')
#     #ax.set_xticklabels(['0 Neg','14 Neg','35 Neg','70 Neg','0 Pos','14 Pos','35 Pos','70 Pos'], fontsize=10, rotation=90)

#     plt.savefig('boxplot_melt_{}.pdf'.format(readtype), bbox_inches='tight')
    
# def bpplot_alt(df,readtype):
#     df = df[df_melt['variable'] == readtype]

#     fig, ax = plt.subplots(figsize=(6,3))

#     odr = ['d0_neg_arm','d14_neg_arm','d35_neg_arm','d70_neg_arm','d0_pos_arm','d14_pos_arm','d35_pos_arm','d70_pos_arm']
#     bps = [('d0_neg_arm','d14_neg_arm'),('d0_neg_arm','d35_neg_arm'),('d0_neg_arm','d70_neg_arm'),('d0_pos_arm','d14_pos_arm'),('d0_pos_arm','d35_pos_arm'),('d0_pos_arm','d70_pos_arm')]
    
#     pal = {'{}_d0'.format(readtype):'#007FFF','{}_d14'.format(readtype):'#00FF7F','{}_d35'.format(readtype):'#FF007F','{}_d70'.format(readtype):'#FFD700'}
    
#     ax = sns.barplot(x="sample_grp", y="value", hue="time", data=df, palette=pal)
#     add_stat_annotation(ax, data=df, x="sample_grp", y="value", box_pairs=bps, test='t-test_ind', text_format='star', loc='inside', verbose=2, linewidth=1)

#     leg = plt.legend(loc='upper left', bbox_to_anchor=(1.03, 1))
#     leg.get_frame().set_linewidth(0.0)

#     # ax.set_yscale('log')
#     # plt.ylim([10**3.5,10**7])

#     ax.set_axisbelow(True)
#     ax.set_title('{} Reads'.format(readtype))
#     ax.set_xlabel('Samples')
#     ax.set_ylabel('Normalized Readcounts')
#     ax.set_xticklabels(['0 Neg','14 Neg','35 Neg','70 Neg','0 Pos','14 Pos','35 Pos','70 Pos'], fontsize=10, rotation=90)

#     plt.savefig('boxplot_melt_alt_{}.pdf'.format(readtype), bbox_inches='tight')
#     plt.close()
    
# for rt in readtypes:
#     bpplot(df_melt,rt)
#     bpplot_alt(df_melt,rt)