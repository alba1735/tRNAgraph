#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
import anndata as ad
import argparse

import directory_tools

import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


def visualizer(adata, output):
    '''
    Generate volcano visualizations for each group in an AnnData object.
    '''
    pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='pca_tools.py',
        description='Generate volcano visualizations for each group in an AnnData object.'
    )

    parser.add_argument('-i', '--anndata',
                        help='Specify AnnData input', required=True)
    parser.add_argument('-o', '--output', 
                        help='Specify output directory', default='heatmap', required=False)

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    directory_tools.builder(args.output)

    adata = ad.read_h5ad(args.anndata)

    visualizer(adata, args.output)
