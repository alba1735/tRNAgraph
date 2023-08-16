#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy import stats
import itertools

def log2fc_df(adata, comparison_groups, readtype, readcount_cutoff):
    '''
    Calculate log2FC of tRNA read counts between groups.
    '''
    df = pd.DataFrame(adata.obs, columns=['trna', comparison_groups, readtype])
    # Create correlation matrixs from reads stored in adata observations as mean and standard deviation
    sdf = df.pivot_table(index='trna', columns=comparison_groups, values=readtype, aggfunc='std')
    mdf = df.pivot_table(index='trna', columns=comparison_groups, values=readtype, aggfunc='mean')
    cdf = df.pivot_table(index='trna', columns=comparison_groups, values=readtype, aggfunc='count')
    # For rows in df if a value is less than readcount_cutoff, drop the row from df
    mean_drop_list = [True if i >= readcount_cutoff else False for i in mdf.mean(axis=1)]
    sdf = sdf[mean_drop_list]
    mdf = mdf[mean_drop_list]
    cdf = cdf[mean_drop_list]
    # Drop rows with NaN values
    sdf = sdf.dropna()
    mdf = mdf.dropna()
    cdf = cdf.dropna()
    # Create permutations of pairings of groups for heatmap
    pairs = list(itertools.combinations(mdf.columns, 2))
    # Create df of log2FC values for each pair from adata.obs nreads_total_raw
    df_pairs = pd.DataFrame()
    for pair in pairs:
        df_pairs[f'log2_{pair[0]}-{pair[1]}'] = np.log2(mdf[pair[1]]) - np.log2(mdf[pair[0]])
        df_pairs[f'pval_{pair[0]}-{pair[1]}'] = stats.ttest_ind_from_stats(mdf[pair[0]], sdf[pair[0]], cdf[pair[0]], mdf[pair[1]], sdf[pair[1]], cdf[pair[1]])[1]

    # sort the columns alphabetically so log2FC are followed by pvals
    df_pairs = df_pairs.reindex(sorted(df_pairs.columns), axis=1)

    return df_pairs