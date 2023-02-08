#!/usr/bin/env python3

import pandas as pd
import numpy as np
import itertools

def log2fc_df(adata, comparison_groups, readtype, readcount_cutoff):
    '''
    Calculate log2FC of tRNA read counts between groups.
    '''
    # Create a correlation matrix from reads stored in adata observations
    df = pd.DataFrame(adata.obs, columns=['trna', comparison_groups, readtype])
    df = df.pivot_table(index='trna', columns=comparison_groups, values=readtype, aggfunc='mean')
    # For rows in df if a value is less than readcount_cutoff, drop the row from df
    df = df[df > readcount_cutoff]
    df = df.dropna()
    # Create permutations of pairings of groups for heatmap
    pairs = list(itertools.combinations(df.columns, 2))
    # Create df of log2FC values for each pair from adata.obs nreads_total_raw
    df_pairs = pd.DataFrame()
    for pair in pairs:
        df_pairs[f'{pair[0]}-{pair[1]}'] = np.log2(df[pair[1]]) - np.log2(df[pair[0]])
    # Sort df_pairs by the sum of the log2FC values for each row
    df_pairs = df_pairs.sort_values(by=df_pairs.columns.tolist(), ascending=False)

    return df_pairs