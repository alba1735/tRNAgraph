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
            # Input format: Feature, Sample, EndType, Count?
            # Or a matrix?
            # Let's assume it's a matrix or long format.
            df = pd.read_csv(self.input_file, sep='\t')
        except Exception as e:
            print(f"Error reading input file: {e}")
            return

        # Check columns
        # If 'end' column exists (common in tRAX output where first col is feature index and second is 'end')
        if 'end' in df.columns:
            df.rename(columns={'end': 'EndType'}, inplace=True)
        elif 'EndType' not in df.columns and df.index.name != 'EndType':
             # Try to guess
             if df.shape[1] > 1 and df.iloc[:, 0].dtype == object:
                 # Maybe first column is EndType?
                 # But here we have Feature as index likely.
                 pass

        if 'EndType' in df.columns:
            # Melt
            id_vars = ['EndType']
            value_vars = [col for col in df.columns if col not in id_vars and col != 'Feature']
            
            melted_df = df.melt(id_vars=id_vars, value_vars=value_vars, var_name='Sample', value_name='Count')
            
            # Aggregate counts by Sample and EndType
            agg_df = melted_df.groupby(['Sample', 'EndType'])['Count'].sum().reset_index()
            
            # Calculate percentages
            # Total count per sample
            sample_totals = agg_df.groupby('Sample')['Count'].transform('sum')
            agg_df['Percentage'] = (agg_df['Count'] / sample_totals) * 100
            
            # Plot
            plt.figure(figsize=(12, 8))
            sns.barplot(data=agg_df, x="Sample", y="Percentage", hue="EndType")
            plt.xticks(rotation=90)
            plt.title("3' End Composition")
            plt.ylabel("Percentage of Reads")
            plt.tight_layout()
            
            plt.savefig(self.output_file, bbox_inches='tight')
            print(f"Plot saved to {self.output_file}")
        else:
            print("Could not identify 'EndType' column.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 plotsLegacyCCA.py <input_file> <output_file>")
        sys.exit(1)
    
    viz = visualizer(sys.argv[1], sys.argv[2])
    viz.plot()
