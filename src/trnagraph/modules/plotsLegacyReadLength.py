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

        # Expected columns: Length, Sample, other, trnas, pretrnas
        if 'Length' not in df.columns or 'Sample' not in df.columns:
            print("Input file missing required columns (Length, Sample)")
            return

        # Binning (bin width = 2, similar to R script)
        bin_width = 2
        df['bin'] = (df['Length'] // bin_width) * bin_width

        # Aggregate tRNAs
        # The R script aggregated 'trnas' by bin and Sample
        plot_df = df.groupby(['bin', 'Sample'])['trnas'].sum().reset_index()

        # Plotting
        # samples = plot_df['Sample'].unique()
        
        # Create a plot per sample or a combined plot?
        # The R script seemed to produce one PDF.
        # Let's make a faceted plot or a single plot with hue.
        
        g = sns.FacetGrid(plot_df, col="Sample", col_wrap=4, sharex=True, sharey=False)
        g.map(sns.barplot, "bin", "trnas", order=sorted(plot_df['bin'].unique()))
        
        # Adjust labels
        for ax in g.axes.flat:
            for label in ax.get_xticklabels():
                label.set_rotation(90)
        
        g.set_axis_labels("Read Length (bin)", "Count")
        g.fig.suptitle("tRNA Read Length Distribution", y=1.02)
        
        plt.savefig(self.output_file, bbox_inches='tight')
        print(f"Plot saved to {self.output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 plotsLegacyReadLength.py <input_file> <output_file>")
        sys.exit(1)
    
    viz = visualizer(sys.argv[1], sys.argv[2])
    viz.plot()
