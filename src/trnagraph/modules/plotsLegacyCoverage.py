#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
from matplotlib.backends.backend_pdf import PdfPages

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

        # Expected columns: Feature, Sample, position, coverage, etc.
        # The exact columns depend on the output of toolsGetCoverage.py
        # Based on tRAX/newcoverageplots.R, it expects: Feature, Sample, position, coverage
        
        required_cols = ['Feature', 'Sample', 'position', 'coverage']
        if not all(col in df.columns for col in required_cols):
            print(f"Input file missing required columns: {required_cols}")
            # Try to guess if columns are different
            print(f"Available columns: {df.columns}")
            return

        # Sort positions
        # Positions can be integers or strings (e.g. 'head', 'tail')
        # We need a robust sorter.
        
        def parse_position(pos):
            pos = str(pos)
            if pos.startswith('head'):
                return -int(pos.replace('head', ''))
            elif pos.startswith('tail'):
                return 1000 + int(pos.replace('tail', ''))
            else:
                try:
                    return int(pos)
                except ValueError:
                    return 0

        df['pos_numeric'] = df['position'].apply(parse_position)
        df = df.sort_values(['Feature', 'pos_numeric'])

        # Get unique features (tRNAs)
        features = df['Feature'].unique()
        
        print(f"Generating coverage plots for {len(features)} features...")
        
        with PdfPages(self.output_file) as pdf:
            # 1. Combined Summary Plot (Average coverage across all tRNAs)
            plt.figure(figsize=(12, 6))
            avg_df = df.groupby(['Sample', 'pos_numeric'])['coverage'].mean().reset_index()
            sns.lineplot(data=avg_df, x='pos_numeric', y='coverage', hue='Sample')
            plt.title("Average tRNA Coverage Profile")
            plt.xlabel("Position")
            plt.ylabel("Mean Coverage")
            plt.tight_layout()
            pdf.savefig()
            plt.close()
            
            # 2. Individual Plots per tRNA
            # To avoid making a massive file if there are hundreds of tRNAs, 
            # we might want to limit this or just do it as requested (legacy behavior was all of them).
            # We'll do all of them.
            
            for feature in features:
                feat_df = df[df['Feature'] == feature]
                
                plt.figure(figsize=(10, 5))
                sns.lineplot(data=feat_df, x='pos_numeric', y='coverage', hue='Sample')
                plt.title(f"Coverage: {feature}")
                plt.xlabel("Position")
                plt.ylabel("Coverage")
                plt.tight_layout()
                pdf.savefig()
                plt.close()
                
        print(f"Coverage plots saved to {self.output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 plotsLegacyCoverage.py <input_file> <output_file>")
        sys.exit(1)
    
    viz = visualizer(sys.argv[1], sys.argv[2])
    viz.plot()
