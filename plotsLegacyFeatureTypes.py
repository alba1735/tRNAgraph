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
        except Exception:
            # Handle R-style dataframe
            try:
                df = pd.read_csv(self.input_file, sep='\t', index_col=0)
                df.reset_index(inplace=True)
                df.rename(columns={df.columns[0]: 'FeatureType'}, inplace=True)
            except Exception as e:
                print(f"Error reading input file with index_col=0: {e}")
                return
        except Exception as e:
            print(f"Error reading input file: {e}")
            return

        # Similar structure to GeneFeatures: Type, Sample1, Sample2...
        if df.columns[0].lower() in ['type', 'featuretype', 'unnamed: 0']:
            df.rename(columns={df.columns[0]: 'FeatureType'}, inplace=True)
        else:
             if df.index.name is None:
                df.index.name = 'FeatureType'
             df = df.reset_index()

        # Melt
        id_vars = ['FeatureType']
        value_vars = [col for col in df.columns if col not in id_vars]
        
        melted_df = df.melt(id_vars=id_vars, value_vars=value_vars, var_name='Sample', value_name='Count')

        # Plotting
        plt.figure(figsize=(14, 8))
        sns.barplot(data=melted_df, x="FeatureType", y="Count", hue="Sample")
        plt.xticks(rotation=90)
        plt.title("Read Counts by Feature Type")
        plt.tight_layout()
        
        plt.savefig(self.output_file, bbox_inches='tight')
        print(f"Plot saved to {self.output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 plotsLegacyFeatureTypes.py <input_file> <output_file>")
        sys.exit(1)
    
    viz = visualizer(sys.argv[1], sys.argv[2])
    viz.plot()
