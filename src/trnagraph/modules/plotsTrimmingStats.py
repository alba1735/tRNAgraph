#!/usr/bin/env python3

import logging

import numpy as np
import pandas as pd

from . import plotsPalette
import matplotlib.pyplot as plt
import seaborn as sns
import os

plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

class visualizer:
    def __init__(self, stats_file, output_file, colormap=None, settings=None):
        self.stats_file = stats_file
        self.output_file = output_file
        self.colormap = colormap
        self.settings = settings
        self.logger = logging.getLogger(__name__)

    def plot(self):
        if not os.path.exists(self.stats_file):
            self.logger.error(f"Stats file not found: {self.stats_file}")
            return

        try:
            df = pd.read_csv(self.stats_file)
        except Exception as e:
            self.logger.error(f"Error reading stats file: {e}")
            return

        # Calculate Discarded
        # Ensure columns exist
        required_cols = ['Sample', 'Raw_Reads', 'Clean_Reads', 'Merged_Reads', 'Unmerged_Reads', 'Trimmed_Reads']
        if not all(col in df.columns for col in required_cols):
            self.logger.error(f"Stats file missing required columns: {required_cols}")
            return

        # 'Merged'/'Unmerged' only apply to paired-end samples (fastp's merge step never runs
        # on single-end input); single-end samples' filter-passing reads land in 'Trimmed'
        # instead, so the two library types never get mislabeled under a shared bucket.
        df['Discarded'] = df['Raw_Reads'] - df['Clean_Reads']
        df['Merged'] = df['Merged_Reads']
        df['Unmerged'] = df['Unmerged_Reads']
        df['Trimmed'] = df['Trimmed_Reads']

        # Select columns for plotting
        plot_df = df[['Sample', 'Merged', 'Unmerged', 'Trimmed', 'Discarded']]

        # Melt
        df_melt = plot_df.melt(id_vars='Sample', var_name='ReadType', value_name='Count')

        # Filter out 0 counts
        df_melt = df_melt[df_melt['Count'] > 0]

        # Filter < 1% (replicating R script logic)
        # Calculate totals per sample
        sample_totals = df_melt.groupby('Sample')['Count'].transform('sum')
        df_melt = df_melt[df_melt['Count'] > sample_totals / 100]

        if df_melt.empty:
            self.logger.warning("No data to plot after filtering.")
            return

        # Pivot back for pandas plotting (stacked bar)
        df_pivot = df_melt.pivot(index='Sample', columns='ReadType', values='Count')

        # Reorder columns if needed: Merged, Unmerged, Trimmed, Discarded
        cols = [c for c in ['Merged', 'Unmerged', 'Trimmed', 'Discarded'] if c in df_pivot.columns]
        df_pivot = df_pivot[cols]

        # Calculate width based on number of samples
        # R script: width = 3 + .25*length(unique(countsmelt$variable))
        num_samples = len(df_pivot.index)
        width = 3 + 0.25 * num_samples
        height = 6 # Default height

        # Create plot
        fig, ax = plt.subplots(figsize=(width, height))

        # Plot stacked bars, husl palette (or --colormap 'trimtype' overrides) to match the
        # rest of tRNAgraph's plots*.py modules -- see plotsCount.py's stacked_barplots
        palette = plotsPalette.categorical_palette(len(cols))
        if self.colormap:
            palette = [self.colormap.get(c, palette[i]) for i, c in enumerate(cols)]
        df_pivot.plot(kind='bar', stacked=True, ax=ax, width=0.8, edgecolor=plotsPalette.BAR_EDGE, linewidth=0.5,
                      color=palette)

        # Styling
        ax.set_xlabel("Sample")
        ax.set_ylabel("Total Reads")
        ax.set_title("Read Types per Sample")

        # X-axis labels
        ax.set_xticklabels(df_pivot.index, rotation=90)

        # Legend, matching the plotsCount.py stacked-barplot convention
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=handles, labels=labels, loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
        ax.legend_.set_title('Read Type')

        # Theme, matching the rest of tRNAgraph's plots*.py modules
        sns.despine(left=True, bottom=True)

        # Adjust layout to prevent clipping
        plt.tight_layout()

        self.logger.info(f"Saving plot to {self.output_file}")
        plt.savefig(self.output_file, bbox_inches='tight')
        plt.close()

if __name__ == '__main__':
    pass
