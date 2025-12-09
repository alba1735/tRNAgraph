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

        # Check for 4-column no-header format (Feature, GeneType, Chromosome, Count)
        try:
            df_check = pd.read_csv(self.input_file, sep='\t', header=None, nrows=5)
            if df_check.shape[1] == 4 and pd.api.types.is_numeric_dtype(df_check.iloc[:, 3]):
                df = pd.read_csv(self.input_file, sep='\t', header=None)
                df.columns = ['Feature', 'GeneType', 'Chromosome', 'Count']
                
                # Aggregate by GeneType
                agg_df = df.groupby('GeneType')['Count'].sum().reset_index()
                agg_df['Sample'] = 'Sample' # Dummy sample name
                
                # Plot
                plt.figure(figsize=(12, 8))
                sns.barplot(data=agg_df, x="Sample", y="Count", hue="GeneType")
                plt.xticks(rotation=90)
                plt.title("Read Counts by Gene Type")
                plt.tight_layout()
                plt.savefig(self.output_file, bbox_inches='tight')
                print(f"Plot saved to {self.output_file}")
                return
        except Exception:
            pass

        try:
            # Assume wide format: Type, Sample1, Sample2...
            try:
                df = pd.read_csv(self.input_file, sep='\t')
            except Exception:
                # Handle R-style dataframe (header has one less column)
                try:
                    df = pd.read_csv(self.input_file, sep='\t', index_col=0)
                    df.reset_index(inplace=True)
                    # The index column name might be empty or 'index', rename it to 'Type' or similar if needed
                    if df.columns[0] == 'index':
                        df.rename(columns={'index': 'Type'}, inplace=True)
                    else:
                        # If the index didn't have a name, reset_index names it 'index'.
                        # If it had a name, it keeps it.
                        # We want a generic name for the feature column
                        df.rename(columns={df.columns[0]: 'Type'}, inplace=True)
                except Exception as e:
                     print(f"Error reading input file with index_col=0: {e}")
                     return
        except Exception as e:
            print(f"Error reading input file: {e}")
            return

        # If the first column is unnamed or 'Type', assume it's the feature type
        if df.columns[0].lower() in ['type', 'feature', 'unnamed: 0']:
            df.rename(columns={df.columns[0]: 'GeneType'}, inplace=True)
        else:
            # If index was read as a column
            if df.index.name is None:
                df.index.name = 'GeneType'
            df = df.reset_index()
            
        # Melt
        id_vars = ['GeneType']
        value_vars = [col for col in df.columns if col not in id_vars]
        
        melted_df = df.melt(id_vars=id_vars, value_vars=value_vars, var_name='Sample', value_name='Count')

        # Plotting
        plt.figure(figsize=(12, 8))
        sns.barplot(data=melted_df, x="Sample", y="Count", hue="GeneType")
        plt.xticks(rotation=90)
        plt.title("Read Counts by Gene Type")
        plt.tight_layout()
        
        plt.savefig(self.output_file, bbox_inches='tight')
        print(f"Plot saved to {self.output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 plotsLegacyGeneFeatures.py <input_file> <output_file>")
        sys.exit(1)
    
    viz = visualizer(sys.argv[1], sys.argv[2])
    viz.plot()
