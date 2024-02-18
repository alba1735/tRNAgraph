#!/usr/bin/env python3

import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats
import itertools
import os

def builder(directory):
    '''
    Function to create output directory if it does not already exist
    '''
    if '/' not in directory: # fix for dumb pathing
        directory = './' + directory
    if not os.path.exists(directory):
        output = f'Creating output directory: {os.path.dirname(directory)}'
        os.makedirs(os.path.dirname(directory), exist_ok=True)
    else:
        output = f'Output directory already exists: {os.path.dirname(directory)}'
    return output

class adataLog2FC():
    def __init__(self, adata, compare, readtype, readcount_cutoff=80, config_name='default'):
        self.adata = adata
        self.compare = compare
        self.readtype = readtype
        self.readcount_cutoff = str(readcount_cutoff)
        self.config_name = config_name

    def main(self):
        # Search nested adata.uns for log2FC for compare, readtype, the readcutoff
        self.log2fc_dict = self.adata.uns.get('log2FC', {})
        self.log2fc_dict[self.config_name] = self.log2fc_dict.get(self.config_name, {})
        self.log2fc_dict[self.config_name][self.compare] = self.log2fc_dict[self.config_name].get(self.compare, {})
        self.log2fc_dict[self.config_name][self.compare][self.readtype] = self.log2fc_dict[self.config_name][self.compare].get(self.readtype, {})
        self.log2fc_dict[self.config_name][self.compare][self.readtype][self.readcount_cutoff] = self.log2fc_dict[self.config_name][self.compare][self.readtype].get(self.readcount_cutoff, {})
        df = self.log2fc_dict[self.config_name][self.compare][self.readtype][self.readcount_cutoff].get('df', pd.DataFrame())
        if df.empty:
            df, pairs = self.log2fc_df()
            self.log2fc_dict[self.config_name][self.compare][self.readtype][self.readcount_cutoff]['df'] = df
            self.log2fc_dict[self.config_name][self.compare]['pairs'] = pairs
            self.adata.uns['log2FC'] = self.log2fc_dict
        else:
            df = self.log2fc_dict[self.config_name][self.compare][self.readtype][self.readcount_cutoff]['df']
            pairs = self.log2fc_dict[self.config_name][self.compare]['pairs']

        return df, self.log2fc_dict

    def log2fc_df(self):
        df = pd.DataFrame(self.adata.obs, columns=['trna', self.compare, self.readtype])
        # Create correlation matrixs from reads stored in adata observations as mean and standard deviation
        sdf = df.pivot_table(index='trna', columns=self.compare, values=self.readtype, aggfunc='std', observed=True)
        mdf = df.pivot_table(index='trna', columns=self.compare, values=self.readtype, aggfunc='mean', observed=True)
        cdf = df.pivot_table(index='trna', columns=self.compare, values=self.readtype, aggfunc='count', observed=True)
        # For rows in df if a value is less than readcount_cutoff, drop the row from df
        mean_drop_list = [True if i >= int(self.readcount_cutoff) else False for i in mdf.mean(axis=1)]
        sdf = sdf[mean_drop_list]
        mdf = mdf[mean_drop_list]
        cdf = cdf[mean_drop_list]
        # Drop rows with NaN values
        sdf = sdf.dropna()
        mdf = mdf.dropna()
        cdf = cdf.dropna()
        # Replace 0 with 0.00000000000000000001 to avoid log2(0) = -inf
        mdf = mdf.replace(0, 0.00000000000000000001)
        # Create permutations of pairings of groups for heatmap
        pairs = list(itertools.combinations(mdf.columns, 2))
        # Create df of log2FC values for each pair from adata.obs nreads_total_raw
        df_pairs = pd.DataFrame()
        for pair in pairs:
            df_pairs[f'log2_{pair[0]}-{pair[1]}'] = np.log2(mdf[pair[1]]) - np.log2(mdf[pair[0]])
            df_pairs[f'pval_{pair[0]}-{pair[1]}'] = stats.ttest_ind_from_stats(mdf[pair[0]], sdf[pair[0]], cdf[pair[0]], mdf[pair[1]], sdf[pair[1]], cdf[pair[1]])[1]
        # sort the columns alphabetically so log2FC are followed by pvals
        df_pairs = df_pairs.reindex(sorted(df_pairs.columns), axis=1)
        return df_pairs, pairs

def log2fc_compare_df(adata, countgrp, comparison_groups, readtype, readcount_cutoff):
    '''
    Calculate log2FC of tRNA read counts between multiple groups.
    '''
    df = pd.DataFrame(adata.obs, columns=[countgrp, *comparison_groups, readtype])
    # Create correlation matrixs from reads stored in adata observations as mean and standard deviation
    sdf = df.pivot_table(index=countgrp, columns=comparison_groups, values=readtype, aggfunc='std', observed=True)
    mdf = df.pivot_table(index=countgrp, columns=comparison_groups, values=readtype, aggfunc='mean', observed=True)
    cdf = df.pivot_table(index=countgrp, columns=comparison_groups, values=readtype, aggfunc='count', observed=True)
    # For rows in df if a value is less than readcount_cutoff, drop the row from df
    mean_drop_list = [True if i >= readcount_cutoff else False for i in mdf.mean(axis=1)]
    sdf = sdf[mean_drop_list]
    mdf = mdf[mean_drop_list]
    cdf = cdf[mean_drop_list]
    # Drop rows with NaN values
    sdf = sdf.dropna()
    mdf = mdf.dropna()
    cdf = cdf.dropna()
    # Make sure all indexs of sdf, mdf, and cdf are the same dropping any that are not
    sdf = sdf[sdf.index.isin(mdf.index)]
    sdf = sdf[sdf.index.isin(cdf.index)]
    mdf = mdf[mdf.index.isin(sdf.index)]
    mdf = mdf[mdf.index.isin(cdf.index)]
    cdf = cdf[cdf.index.isin(sdf.index)]
    cdf = cdf[cdf.index.isin(mdf.index)]
    # Replace 0 with 0.00000000000000000001 to avoid log2(0) = -inf
    mdf = mdf.replace(0, 0.00000000000000000001)
    # Create a dict where the keys are the first level of the multiindex and the values are the second level
    pairs = {i: list(mdf.columns.get_level_values(1)[mdf.columns.get_level_values(0) == i]) for i in mdf.columns.get_level_values(0).unique()}
    # Subset the value lists in each key-value pair to only include values found in all key-pair values
    ppairs = list(itertools.combinations(set.intersection(*map(set,pairs.values())), 2))
    pairs = list(pairs.keys())
    # Sort the tuples in ppairs so the first element is always alphabetical
    ppairs = [tuple(sorted(i)) for i in ppairs]
    # Print combinations of comparison groups as a tuple for multiindex columns
    pppairs = list(itertools.product(pairs,ppairs))
    pppairs = [tuple(['log2',i[0],i[1][0]+'-'+i[1][1]]) for i in pppairs] + [tuple(['pval',i[0],i[1][0]+'-'+i[1][1]]) for i in pppairs]
    # Create empty df of log2FC values for each pair with a multiindex column
    df_pairs = pd.DataFrame(0, index=mdf.index, columns=pd.MultiIndex.from_tuples(pppairs, names=['stats', 'cgrp1', 'cgrp2']))
    for cgrp1 in pairs:
        for cgrp2 in ppairs:
            df_pairs.loc[:, ('log2', cgrp1, cgrp2[0]+'-'+cgrp2[1])] = np.log2(mdf[cgrp1][cgrp2[1]]) - np.log2(mdf[cgrp1][cgrp2[0]])
            df_pairs.loc[:, ('pval', cgrp1, cgrp2[0]+'-'+cgrp2[1])] = stats.ttest_ind_from_stats(mdf[cgrp1][cgrp2[0]], sdf[cgrp1][cgrp2[0]], cdf[cgrp1][cgrp2[0]], 
                                                                                                 mdf[cgrp1][cgrp2[1]], sdf[cgrp1][cgrp2[1]], cdf[cgrp1][cgrp2[1]])[1]

    return df_pairs