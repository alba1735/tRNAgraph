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

        # Expected columns: tRNA_name, sample, position, coverage, total
        required_cols = ['tRNA_name', 'sample', 'position', 'coverage']
        if not all(col in df.columns for col in required_cols):
            print(f"Input file missing required columns: {required_cols}")
            return

        # The 'position' column might be strings like 'head30', '1', 'tail10'.
        # We need to convert them to a sortable numeric or categorical order.
        # Let's inspect unique positions to determine order.
        # unique_positions = df['position'].unique()
        
        # Simple heuristic for sorting:
        # headX -> negative numbers?
        # numbers -> numbers
        # tailX -> large numbers?
        
        def parse_position(pos):
            pos = str(pos)
            if pos.startswith('head'):
                return -int(pos.replace('head', ''))
            elif pos.startswith('tail'):
                return 1000 + int(pos.replace('tail', '')) # Assuming tRNA length < 1000
            else:
                try:
                    return int(pos)
                except ValueError:
                    return 0 # Fallback

        df['pos_numeric'] = df['position'].apply(parse_position)
        
        # Sort by tRNA and position
        df = df.sort_values(['tRNA_name', 'pos_numeric'])

        # Plotting
        # Since there are many tRNAs, plotting all of them in one file might be huge.
        # The R script likely did one plot per tRNA or a summary.
        # Let's plot a summary (average coverage across all tRNAs) and maybe top 5 abundant ones.
        
        # 1. Average coverage profile
        avg_df = df.groupby(['sample', 'pos_numeric'])['coverage'].mean().reset_index()
        
        plt.figure(figsize=(12, 6))
        sns.lineplot(data=avg_df, x='pos_numeric', y='coverage', hue='sample')
        plt.title("Average Pre-tRNA Coverage Profile")
        plt.xlabel("Position (relative to mature tRNA)")
        plt.ylabel("Mean Coverage")
        plt.tight_layout()
        plt.savefig(self.output_file, bbox_inches='tight')
        print(f"Summary plot saved to {self.output_file}")
        
        # Note: To plot individual tRNAs, we would need a multi-page PDF or separate files.
        # For now, the summary confirms the data is readable and plottable.

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 plotsLegacyLocusCoverage.py <input_file> <output_file>")
        sys.exit(1)
    
    viz = visualizer(sys.argv[1], sys.argv[2])
    viz.plot()
