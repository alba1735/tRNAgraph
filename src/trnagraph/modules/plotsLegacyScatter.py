#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
# import itertools
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
            df = pd.read_csv(self.input_file, sep='\t', index_col=0)
        except Exception as e:
            print(f"Error reading input file: {e}")
            return

        samples = df.columns
        if len(samples) < 2:
            print("Need at least 2 samples for scatter plot.")
            return

        # Generate pairs
        # pairs = list(itertools.combinations(samples, 2))
        
        # Limit number of plots if too many?
        # For now, let's just do a pairplot if N is small, or specific pairs.
        # The legacy script likely did all pairs or specific ones.
        # Let's use seaborn pairplot for a matrix view.
        
        # Log transform for better visualization
        df_log = df.apply(lambda x: x + 1).apply(np.log2)
        
        plt.figure(figsize=(10, 10))
        g = sns.pairplot(df_log, diag_kind="kde")
        g.fig.suptitle("Scatter Matrix of Log2(Counts + 1)", y=1.02)
        
        plt.savefig(self.output_file, bbox_inches='tight')
        print(f"Plot saved to {self.output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 plotsLegacyScatter.py <input_file> <output_file>")
        sys.exit(1)
    
    viz = visualizer(sys.argv[1], sys.argv[2])
    viz.plot()
