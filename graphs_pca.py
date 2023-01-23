#!/usr/bin/env python3

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns

import pandas as pd
import anndata as ad
import os
import argparse

class grapher():
    def __init__(self, adata, pcagroup, output):
        self.adata = adata
        # self.pcagroup = pcagroup
        self.pcagroup = 'timepoint'
        self.output = output

    def create(self):
        # subset adata to unique coverage and drop gap positions and convert to df
        adata_u = self.adata[:,self.adata.var.coverage=='uniquecoverage']
        adata_u = adata_u[:,adata_u.var.gap==False]
        df = adata_u.to_df()
        # grab number of reads per sample and add sample groups
        df = df.max(axis=1).to_frame('n_reads')
        df['sample'] = adata_u.obs['sample'].values


        # df = df.T
        # df = StandardScaler().fit_transform(df)

        pca = PCA(n_components=2)
        pca.fit(df)

        print(pca.explained_variance_ratio_)
        print(pca.singular_values_)

        df_pca = pd.DataFrame(pca.transform(df), columns=['PC1', 'PC2'])
        # df_pca['group'] = str(adata_u.obs[self.pcagroup].values) + ' ' + str(adata_u.obs['sample'].values)

        sns.scatterplot(data=df_pca, x='PC1', y='PC2', palette='tab10')

        plt.savefig('{}/{}_pca.pdf'.format(self.output, self.pcagroup))

        print(adata_u.obs['sample'].values.unique())










def directory_builder(directory):
    if not os.path.exists(directory):
        print('Creating output directory: {}'.format(directory))
        os.makedirs(directory, exist_ok=True)
    else:
        print('Output directory already exists: {}'.format(directory))       

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='graphs_pca',
        description='Generate PCA graphs for each sample in an AnnData object.'
    )

    parser.add_argument('-i', '--anndata', help='Specify AnnData input', required=True)
    parser.add_argument('-x', '--pcagroup', help='Specify PCA grouping', required=True)
    parser.add_argument('-o', '--output', help='Specify output directory', default='pca', required=False)

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    directory_builder(args.output)

    adata = ad.read_h5ad(args.anndata)

    grapher(adata, args.pcagroup, args.output).create()