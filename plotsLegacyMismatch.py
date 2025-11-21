#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

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
            df = pd.read_csv(self.input_file, sep='\t')
        except Exception as e:
            print(f"Error reading input file: {e}")
            return

        # Expected columns: count, type, Sample1, Sample2, ...
        if 'count' not in df.columns or 'type' not in df.columns:
            print("Input file missing required columns (count, type)")
            return

        # Melt the dataframe
        id_vars = ['count', 'type']
        value_vars = [col for col in df.columns if col not in id_vars]
        
        melted_df = df.melt(id_vars=id_vars, value_vars=value_vars, var_name='Sample', value_name='Value')

        # Plotting
        # Facet by Sample, x=count, y=Value, hue=type
        g = sns.FacetGrid(melted_df, col="Sample", col_wrap=4, sharex=True, sharey=False)
        g.map(sns.barplot, "count", "Value", "type", order=sorted(melted_df['count'].unique()), palette="viridis")
        g.add_legend()
        
        g.set_axis_labels("Mismatches", "Count (Normalized)")
        g.fig.suptitle("Read Mismatches", y=1.02)
        
        plt.savefig(self.output_file, bbox_inches='tight')
        print(f"Plot saved to {self.output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 plotsLegacyMismatch.py <input_file> <output_file>")
        sys.exit(1)
    
    viz = visualizer(sys.argv[1], sys.argv[2])
    viz.plot()
