#!/usr/bin/env python3

import logging

import numpy as np
import pandas as pd

from . import toolsTG
from . import plotsPalette

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns

logger = logging.getLogger(__name__)

def _log2fc_by_stratum(adata, countgrp, comparegrp1, comparegrp2, readtype):
    '''
    Fold changes between --comparegrp2 values, computed separately within each --comparegrp1
    value, as the MultiIndex ('stats', 'cgrp1', 'cgrp2') frame this plot draws from.

    The fit is adataLog2FC's PyDESeq2 negative-binomial model -- the same engine the heatmap
    and volcano use -- rather than the t-test over per-group means and SDs this plot used to
    carry, so two plots of one dataset can no longer disagree about the same contrast. It is
    called through log2fc_df() rather than main() deliberately: the uns['log2FC'] cache is
    keyed [config][compare][readtype][cutoff] with no slot for a feature axis or a stratum,
    and adding a level would invalidate the cached entries in every already-built object, for
    a plot that is opt-in, excluded from -g all, and cheap to fit at 22-54 features.

    Stratifying by subsetting, rather than modelling ~grp1+grp2, is also deliberate: a real
    interaction term belongs to the planned multivariate engine, not here.
    '''
    strata = sorted(adata.obs[comparegrp1].dropna().unique())
    per_stratum = {}
    for cgrp1 in strata:
        subset = adata[adata.obs[comparegrp1] == cgrp1]
        # A negative-binomial fit needs residual degrees of freedom to estimate dispersion
        # from, so a stratum with one sample per --comparegrp2 value cannot be modelled at
        # all. Checked here, by name, rather than left to PyDESeq2's "no replicates to
        # estimate the dispersion" ValueError, which names neither the column nor the stratum
        # and would abort the whole graph run rather than this one plot. Note the previous
        # t-test engine did not handle this case either -- its per-cell standard deviation was
        # NaN and dropna() emptied the frame, so it wrote nothing and said nothing.
        samples = subset.obs.drop_duplicates('sample')
        levels = samples[comparegrp2].nunique()
        if len(samples) <= levels:
            logger.warning(
                f'Skipping the {countgrp} compare plot: {comparegrp1} "{cgrp1}" has '
                f'{len(samples)} samples across {levels} {comparegrp2} values, leaving no '
                f'replicates to estimate dispersion from. Fold changes need at least one '
                f'{comparegrp2} value with two or more samples inside each {comparegrp1}.'
            )
            return pd.DataFrame()
        # readcount_cutoff=0: compare has never applied one, and the axis is amino/iso, where
        # a cutoff tuned for per-tRNA counts would mean something quite different.
        frame, pairs = toolsTG.adataLog2FC(
            subset, compare=comparegrp2, readtype=readtype, readcount_cutoff=0,
        ).log2fc_df(index_col=countgrp)
        per_stratum[cgrp1] = (frame, [f'{a}-{b}' for a, b in pairs])

    # Only pairs present in EVERY stratum, and only features present in every stratum: the
    # plot indexes (feature, cgrp1, cgrp2) for all three axes, so a pair or feature missing
    # from one stratum has no cell to draw. This is what the previous implementation's
    # set-intersection over pairs and dropna() over features amounted to.
    shared_pairs = sorted(set.intersection(*(set(p) for _, p in per_stratum.values()))) if per_stratum else []
    shared_index = None
    for frame, _ in per_stratum.values():
        shared_index = frame.index if shared_index is None else shared_index.intersection(frame.index)
    if not shared_pairs or shared_index is None or len(shared_index) == 0:
        return pd.DataFrame()

    columns = pd.MultiIndex.from_tuples(
        [(stat, cgrp1, pair) for stat in ('log2', 'pval') for cgrp1 in strata for pair in shared_pairs],
        names=['stats', 'cgrp1', 'cgrp2'],
    )
    out = pd.DataFrame(0.0, index=shared_index, columns=columns)
    for cgrp1, (frame, _) in per_stratum.items():
        for pair in shared_pairs:
            out.loc[:, ('log2', cgrp1, pair)] = frame[f'log2_{pair}'].reindex(shared_index)
            out.loc[:, ('pval', cgrp1, pair)] = frame[f'pval_{pair}'].reindex(shared_index)
    return out


def visualizer(adata, comparegrp1, comparegrp2, colormap, output, threaded=True, read_basis=toolsTG.READ_BASIS_UNIQUE, settings=None):
    # Fall back to 'sample' (with a warning) if the specified columns aren't in the AnnData object
    comparegrp1 = toolsTG.resolve_grp_column(adata, comparegrp1, 'comparegrp1')
    comparegrp2 = toolsTG.resolve_grp_column(adata, comparegrp2, 'comparegrp2')
    # Create a color palette for the p
    if colormap != None:
        pal = {k:v if v[0]!='#' else mplcolors.to_rgb(v) for k,v in colormap.items()}
    else:
        pal = plotsPalette.categorical(settings, len(adata.obs[comparegrp1].unique()))
        pal = dict(zip(sorted(adata.obs[comparegrp1].unique()), pal))

    # Was hardcoded to all reads while plotsCluster was hardcoded to unique, so two plots of
    # one dataset rested on different denominators with nothing on either saying so.
    readtype = toolsTG.resolve_readtype('total', read_basis, adata)

    for countgrp in ['amino','iso']:
        # Get log2 fold change dataframe from analysis_tools
        df = _log2fc_by_stratum(adata, countgrp, comparegrp1, comparegrp2, readtype)
        # Zero valid comparegrp1/comparegrp2 pairs is legitimate (e.g. comparegrp1 has no
        # comparegrp2 value shared across every one of its own values -- notably always true
        # when comparegrp1 is a per-observation-unique column like 'sample', the Parameter
        # Fallback default). The frame is empty and has no column levels at all in that case,
        # so skip this countgrp instead of a KeyError on the .loc[:, ('log2')] below.
        if df.empty or 'log2' not in df.columns.get_level_values('stats'):
            logger.warning(
                f'WARNING: no valid {comparegrp1}/{comparegrp2} pairs found for {countgrp}; skipping compare plot.'
            )
            continue
        # Sort the df by the mean of the log2 fold change
        df = df.loc[df.loc[:, ('log2')].abs().mean(axis=1).sort_values(ascending=True).index, :]
        cgrp1list = df.columns.get_level_values('cgrp1').unique()
        cgrp2list = df.columns.get_level_values('cgrp2').unique()
        # Create a bar widths table for amount of values in cgrp1
        barwidths = np.linspace(0.1, 0.9, len(cgrp1list)+1)
        bardiff = np.diff(barwidths)[0]
        barwidths = {k:v for k,v in zip(cgrp1list, barwidths)}
        # Bar edges default to the bar spacing itself, which is how this plot was tuned; a
        # --style line_width replaces that with an absolute width.
        edgewidth = toolsTG.linewidth_for(settings, bardiff)
        # Enumerate the index of the dataframe so that the bar heights can be adjusted
        en_dict = dict(enumerate(df.index.values))
        for cgrp2 in cgrp2list:
            # One figure per comparison, created INSIDE this loop. Created once outside it,
            # the axes were reused for every pair, so the second figure carried the first
            # pair's bars as well as its own and the third carried both -- every file after
            # the first held cumulative content that looked like real data.
            fig, ax = plt.subplots(figsize=toolsTG.figsize_for(settings, (8, 12)))
            xminmax = tuple([-np.abs(df.loc[:, ('log2')]).max().max()*1.1, np.abs(df.loc[:, ('log2')]).max().max()*1.1])
            for cgrp1 in cgrp1list:
                for y,posname in en_dict.items():
                    if abs(df.loc[posname, ('log2',cgrp1,cgrp2)]) >= 1:
                        ax.barh(y+barwidths[cgrp1], df.loc[posname, ('log2',cgrp1,cgrp2)], color=pal[cgrp1], align='edge',
                                height=bardiff, linewidth=edgewidth, edgecolor=pal[cgrp1], label=cgrp1)
                        # ax.hlines(y+0.25, xmin=xminmax[0], xmax=xminmax[1], color=plotsPalette.GRID_LINE, linewidth=bardiff, zorder=-1)
                        # ax.hlines(y+0.75, xmin=xminmax[0], xmax=xminmax[1], color=plotsPalette.GRID_LINE, linewidth=bardiff, zorder=-1)
                        ax.barh(y+0.4, xminmax, color=plotsPalette.GRID_LINE, align='edge', height=0.2, linewidth=0, zorder=-2)
                    else:
                        ax.barh(y+barwidths[cgrp1], df.loc[posname, ('log2',cgrp1,cgrp2)], color='white', align='edge',
                                height=bardiff, linewidth=edgewidth, edgecolor=pal[cgrp1], label=cgrp1)
                    ax.hlines(y+0.5, xmin=xminmax[0], xmax=xminmax[1], color=plotsPalette.GRID_LINE, linewidth=0.5, zorder=1)
            # Set the xlim to the xminmax
            ax.set_xlim(xminmax)
            ax.set_xlabel('Log2 Fold-Change')
            # Set the yticks to the index values + 0.5
            ax.set_yticks(np.arange(len(df.index))+0.5)
            # Set the yticklabels to the index values
            ax.set_yticklabels(df.index)
            ax.set_ylabel(countgrp.capitalize())
            # Set ymin and ymax to -0.5 and len(df.index)+0.5
            ax.set_ylim(-0.5, len(df.index)+0.5)
            # Add light gray vertical lines at each integer
            ax.vlines(np.arange(round(xminmax[0]), round(xminmax[1])), ymin=-0.5, ymax=len(df.index)+0.5, color=plotsPalette.GRID_LINE,linewidth=0.5, zorder=-2)
            # Add a legend made manually from the bar colors
            handles = [plt.Rectangle((0,0),1,1, color=pal[i]) for i in cgrp1list] + \
                [plt.Rectangle((0,0),1,1, color='white')] + \
                [plt.Rectangle((0,0),1,1, color='black', linewidth=bardiff)] + \
                [plt.Rectangle((0,0),1,1, facecolor='white', linewidth=bardiff, edgecolor='black')]
            labels = cgrp1list.to_list() + ['','log2FC >= 1','log2FC < 1']
            ax.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.25, 1), frameon=False)
            # Add a title
            ax.set_title(f'{cgrp2} by {comparegrp1} {countgrp.capitalize()} Log2 Fold-Change')
            # Save the figure
            outpath = f'{output}{comparegrp2}_{cgrp2}_by_{comparegrp1}_{countgrp}_log2fc.pdf'
            toolsTG.save_current(outpath, settings)
            plt.close(fig)
            # Accumulated, not returned: returning here ended the whole function after the
            # first saved figure whenever `threaded` was truthy, so a pooled run produced one
            # file where the serial run produced all of them.
            if threaded:
                threaded += f'Saving figure to {outpath}...\n'
            else:
                logger.info(f'Saving figure to {outpath}...')

    return threaded

if __name__ == '__main__':
    pass