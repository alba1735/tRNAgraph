#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
# import numpy as np

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

        # Expected columns: position, coverage, mismatchedbases
        required_cols = ['position', 'coverage', 'mismatchedbases']
        if not all(col in df.columns for col in required_cols):
            print(f"Input file missing required columns: {required_cols}")
            return

        # Calculate percent mismatch
        # Avoid division by zero
        df = df[df['coverage'] > 0]
        df['percent_mismatch'] = (df['mismatchedbases'] / df['coverage']) * 100

        # Filter positions?
        # The R script had a specific position order.
        # We'll just sort numerically if possible.
        
        def parse_position(pos):
            try:
                return int(pos)
            except ValueError:
                return -1 # e.g. 'e1', 'e2' -> -1 or handle separately

        df['pos_numeric'] = df['position'].apply(parse_position)
        
        # Filter out non-numeric positions for now (or handle 'e' positions if needed)
        # The R script handled 'e' positions (extensions?).
        # Let's keep only numeric positions 1-76 for simplicity, or all.
        
        # Plotting
        # Boxplot of percent_mismatch vs position
        # Facet by Sample
        
        # To make it readable, maybe just plot mean/median line with error band?
        # Boxplot for every position is crowded.
        # But the script is named 'boxplotmismatches'.
        
        g = sns.FacetGrid(df, col="Sample", col_wrap=4, sharex=True, sharey=True)
        
        # Using lineplot as a cleaner alternative to boxplot for many positions
        g.map(sns.lineplot, "pos_numeric", "percent_mismatch")
        
        g.set_axis_labels("Position", "Mismatch %")
        g.fig.suptitle("Mismatch Percentage per Position", y=1.02)
        
        plt.savefig(self.output_file, bbox_inches='tight')
        print(f"Plot saved to {self.output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 plotsLegacyMismatchBoxplot.py <input_file> <output_file>")
        sys.exit(1)
    
    viz = visualizer(sys.argv[1], sys.argv[2])
    viz.plot()
