#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

class visualizer:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file

    def plot(self):
        if not os.path.exists(self.input_file):
            print(f"Input file not found: {self.input_file}")
            return

        try:
            # Input is likely a counts matrix: Genes x Samples
            df = pd.read_csv(self.input_file, sep='\t', index_col=0)
        except Exception as e:
            print(f"Error reading input file: {e}")
            return

        # Transpose so samples are rows
        df_t = df.T
        
        # Filter out zero variance columns
        df_t = df_t.loc[:, (df_t != df_t.iloc[0]).any()]
        
        if df_t.empty:
            print("No data for PCA.")
            return

        # Standardize
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(df_t)
        
        # PCA
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(scaled_data)
        
        pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'], index=df_t.index)
        pca_df['Sample'] = pca_df.index
        
        # Explained variance
        var_exp = pca.explained_variance_ratio_
        
        # Plotting
        plt.figure(figsize=(8, 8))
        sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Sample', s=100)
        
        plt.xlabel(f'PC1 ({var_exp[0]*100:.2f}%)')
        plt.ylabel(f'PC2 ({var_exp[1]*100:.2f}%)')
        plt.title("PCA of Read Counts")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        
        plt.savefig(self.output_file, bbox_inches='tight')
        print(f"Plot saved to {self.output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 plotsLegacyPCA.py <input_file> <output_file>")
        sys.exit(1)
    
    viz = visualizer(sys.argv[1], sys.argv[2])
    viz.plot()
