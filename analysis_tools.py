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
    sdf = sdf[mdf >= readcount_cutoff]
    mdf = mdf[mdf >= readcount_cutoff]
    cdf = cdf[mdf >= readcount_cutoff]
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
    # Sort df_pairs by the sum of the log2FC values for each row
    df_pairs = df_pairs.sort_values(by=df_pairs.columns.tolist(), ascending=False)
    # sort the columns alphabetically
    df_pairs = df_pairs.reindex(sorted(df_pairs.columns), axis=1)

    return df_pairs