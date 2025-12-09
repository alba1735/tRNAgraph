#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
import numpy as np

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
            # Expected columns: log2FoldChange, pvalue, padj, etc.
            df = pd.read_csv(self.input_file, sep='\t')
        except Exception as e:
            print(f"Error reading input file: {e}")
            return

        # Check for required columns
        if 'log2FoldChange' not in df.columns:
             print("Input file missing 'log2FoldChange' column.")
             return
             
        if 'pvalue' not in df.columns:
             if 'padj' in df.columns:
                 df.rename(columns={'padj': 'pvalue'}, inplace=True)
             else:
                 print("Input file missing 'pvalue' or 'padj' column.")
                 return

        # Remove NA
        df = df.dropna(subset=['log2FoldChange', 'pvalue'])
        
        # Calculate -log10(pvalue)
        df['neg_log10_pvalue'] = -np.log10(df['pvalue'])
        
        # Color by significance
        # p < 0.05 and abs(log2FC) > 1
        df['Significance'] = 'NS'
        df.loc[(df['pvalue'] < 0.05) & (df['log2FoldChange'] > 1), 'Significance'] = 'Up'
        df.loc[(df['pvalue'] < 0.05) & (df['log2FoldChange'] < -1), 'Significance'] = 'Down'
        
        colors = {'NS': 'grey', 'Up': 'red', 'Down': 'blue'}
        
        plt.figure(figsize=(8, 8))
        sns.scatterplot(data=df, x='log2FoldChange', y='neg_log10_pvalue', hue='Significance', palette=colors, alpha=0.7)
        
        plt.title("Volcano Plot")
        plt.xlabel("Log2 Fold Change")
        plt.ylabel("-Log10 P-value")
        plt.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.8)
        plt.axvline(1, color='black', linestyle='--', linewidth=0.8)
        plt.axvline(-1, color='black', linestyle='--', linewidth=0.8)
        
        plt.tight_layout()
        plt.savefig(self.output_file, bbox_inches='tight')
        print(f"Plot saved to {self.output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 plotsLegacyVolcano.py <input_file> <output_file>")
        sys.exit(1)
    
    viz = visualizer(sys.argv[1], sys.argv[2])
    viz.plot()
