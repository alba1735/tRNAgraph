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
            # The file format has N columns in header and N+1 in data.
            # Pandas usually handles this by using the first column as index.
            df = pd.read_csv(self.input_file, sep='\t')
        except Exception as e:
            print(f"Error reading input file: {e}")
            return

        # If the index name is None, give it a name
        if df.index.name is None:
            df.index.name = 'FeatureType'
        
        # Reset index to make FeatureType a column
        df = df.reset_index()
        
        # Melt
        id_vars = ['FeatureType']
        value_vars = [col for col in df.columns if col not in id_vars]
        
        melted_df = df.melt(id_vars=id_vars, value_vars=value_vars, var_name='Sample', value_name='Count')

        # Plotting
        # Stacked bar chart per sample? Or grouped bar chart?
        # The R script 'featuretypesreal.R' did a bar plot.
        
        # Let's do a stacked bar plot where x=Sample, y=Count, fill=FeatureType
        
        plt.figure(figsize=(12, 8))
        sns.barplot(data=melted_df, x="Sample", y="Count", hue="FeatureType")
        plt.xticks(rotation=90)
        plt.title("Feature Types (Real Counts)")
        plt.tight_layout()
        
        plt.savefig(self.output_file, bbox_inches='tight')
        print(f"Plot saved to {self.output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 plotsLegacyFeatureTypesReal.py <input_file> <output_file>")
        sys.exit(1)
    
    viz = visualizer(sys.argv[1], sys.argv[2])
    viz.plot()
