#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
import anndata as ad
import argparse

import directory_tools
import analysis_tools

import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


def visualizer(adata, grp, readtype, cutoff, output):
    '''
    Generate volcano visualizations for each group in an AnnData object.
    '''
    # Create a correlation matrix from reads stored in adata observations
    df = analysis_tools.log2fc_df(adata, grp, readtype, cutoff)
    print(df)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='pca_tools.py',
        description='Generate volcano visualizations for each group in an AnnData object.'
    )

    parser.add_argument('-i', '--anndata', help='Specify AnnData input', required=True)
    parser.add_argument('-o', '--output', help='Specify output directory', default='volcano', required=False)
    parser.add_argument('--volgrp', help='Specify group to use for volcano plot', default='group', required=False)
    parser.add_argument('--volrt', help='Specify readtype to use for volcano plot', default='nreads_total_norm', required=False)
    parser.add_argument('--volcutoff', help='Specify readcount cutoff to use for volcano plot', default=80, required=False)

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    directory_tools.builder(args.output)

    adata = ad.read_h5ad(args.anndata)

    visualizer(adata, args.volgrp, args.volrt, args.volcutoff, args.output)
