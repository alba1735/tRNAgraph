#!/usr/bin/env python3

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns

class visualizer:
    def __init__(self, adata, radargrp, radarmethod, radarscaled, colormap, output, threaded=False):
        self.adata = adata
        self.radargrp = radargrp
        self.radarmethod = radarmethod
        self.radarscaled = radarscaled
        self.colormap = colormap
        self.output = output
        self.threaded = threaded

    def isotype_plots(self):
        for readtype in ['nreads_total_unique_norm', 'nreads_total_norm']:
            df = self.adata.uns['anticodon_counts'].copy()
            df = pd.DataFrame(self.adata.obs, columns=['iso', 'amino', self.radargrp, readtype])
            # Drop rows where iso is NNN
            df = df[df['iso']!='NNN']
            # Create dict from all iso and amino where each amino is a key and each value is a list of iso that have that amino
            amino_dict = {}
            for amino in df['amino'].unique():
                amino_dict[amino] = df[df['amino']==amino]['iso'].unique().tolist()
            # Convert the df to a pivot table and use the specified stats method to aggregate the data
            df = pd.pivot_table(df, values=readtype, index=['iso'], columns=[self.radargrp], aggfunc=self.radarmethod, observed=True)
            # Only create radar plots for amino acids with at least 3 isoforms

            # Drop ACT if in the dataframe -- issue with non-standard anticodon
            if 'ACT' in df.index:
                df = df.drop('ACT')

            for aa, codons, in amino_dict.items():
                if len(codons) >= 3:
                    # Subset dataframe to only include codons for the current amino acid
                     self.generate_plot(df[df.index.isin(codons)], readtype, aminoacid=aa)
            # Create radar plot for all isoforms
            self.generate_plot(df, readtype)

    def generate_plot(self, df, readtype, aminoacid=False):
        '''
        Generate radar plots for each sample in an AnnData object.
        '''
        # Normalize counts to 100
        tdf = 100*df/df.sum(axis=0)
        if not aminoacid:
            # Drop rows where the mean of the row comprises less than 1% of the total reads
            tdf = tdf[tdf.min(axis=1)>1]
        # Sort the dataframe by index
        tdf = tdf.sort_index()
        # Replace negative values with 0
        # tdf[tdf<0] = 0
        # Set angles for each category
        angles = np.linspace(0, 2*np.pi, num=len(list(tdf.T))+1)[:-1]
        angles = angles.tolist()
        # Create figure
        fsize = (12,12) if not aminoacid else (6,6)
        fig, ax = plt.subplots(figsize=fsize, subplot_kw=dict(polar=True))
        # Set theta zero location so that anticodons are in the correct position
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        # Set plot parameters
        plt.xticks(angles, list(tdf.T), color='black', size=12)
        plt.yticks(color='dimgrey', size=8)
        ax.set_rlabel_position(0)
        ax.xaxis.grid(False)
        # Set the maximum value for the y-axis to 100
        if self.radarscaled:
            ax.set_ylim(0,100)
        # Set colormap if specified
        if self.colormap != None:
            self.colormap = {k:v if v[0]!='#' else mplcolors.to_rgb(v) for k,v in self.colormap.items()}
            for v in tdf.columns:
                if v not in self.colormap:
                    if self.threaded:
                        self.threaded += f'Color {v} not found in colormap. Using default colors instead.\n'
                    else:
                        print(f'Color {v} not found in colormap. Using default colors instead.')
                    self.colormap = None
                    break
        # Plot data
        for i in tdf.columns.values:
            v = tdf[i].values.tolist()
            v = np.concatenate((v, v[:1]))
            a = angles + angles[:1]
            # Plot data for each anticodon group
            if self.colormap != None:
                ax.plot(a, v, linewidth=1.5, linestyle='solid', label=i, color=self.colormap[i])
                ax.fill(a, v, alpha=0.5/len(tdf.columns.values), color=self.colormap[i])
            else:
                pal = dict(zip(tdf.columns, sns.color_palette('husl', len(tdf.columns))))
                ax.plot(a, v, linewidth=1.5, linestyle='solid', label=i, color=pal[i])
                ax.fill(a, v, alpha=0.5/len(tdf.columns.values), color=pal[i])
        # Capatilize the legend and move the legend outside the plot and remove the border around it
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=handles, labels=[x.capitalize() for x in labels])
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
        ax.legend_.set_title(self.radargrp.capitalize())
        # Give plot a title
        if aminoacid:
            plt.title(f'{aminoacid} Codon Pool Percent Expression Change')
        else:
            plt.title('Codon Pool Percent Expression Change')
        # Save plot
        outname = ''
        if aminoacid:
            outname += f'{self.output}{aminoacid}_radar_by_{self.radargrp}'
        else:
            outname += f'{self.output}all_gt1percent_radar_by_{self.radargrp}'
        if self.radarscaled:
            outname += '_scaled'
        if readtype == 'nreads_total_unique_norm':
            outname += '_unique'
        outname += f'_{self.radarmethod}.pdf'
        if self.threaded:
            self.threaded += f'Plot saved to {outname}\n'
        else:   
            print(f'Plot saved to {outname}')
        plt.savefig(outname, bbox_inches='tight')
        plt.close()
        if self.threaded:
            return self.threaded

if __name__ == '__main__':
    pass