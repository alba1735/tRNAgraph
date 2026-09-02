#!/usr/bin/env python3

import logging

import seaborn as sns
import numpy as np
import pandas as pd
import os

from . import toolsTG
from . import plotsThresholds
from . import plotsPalette

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

logger = logging.getLogger(__name__)


#: Ceiling of the significance panel's colour scale: the scale saturates exactly at the
#: strong-evidence threshold, so a fully-saturated cell reads as "at least this significant".
SIGNIFICANCE_VMAX = plotsThresholds.neglog10(plotsThresholds.PVAL_STRONG)


def significance_colorbar_ticks():
    '''
    Positions and labels for the significance panel's colorbar.

    The panel plots -log10(padj), so the ticks are derived positions rather than the p-values
    themselves. Written as literals (`[0, 1.3, 3]`) these restated plotsVolcano's own copy of
    the same pair and inherited its rounded -log10(0.05); deriving them keeps the tick, the
    scale ceiling and the volcano's reference line tied to one definition.
    '''
    ticks = [0.0,
             plotsThresholds.neglog10(plotsThresholds.PVAL_SIG),
             plotsThresholds.neglog10(plotsThresholds.PVAL_STRONG)]
    labels = ['1', f'{plotsThresholds.PVAL_SIG:g}', f'<={plotsThresholds.PVAL_STRONG:g}']
    return ticks, labels


def visualizer(adata, grp, readtypes, cutoff, heatbound, heatsubplots, output, threaded=False, config_name='default', overwrite=False, settings=None, orientation='vertical'):
    '''
    Generate heatmap visualizations for each group in an AnnData object. `readtypes` are
    fully-resolved obs column names (e.g. 'nreads_total_unique_norm'), already carrying the
    read basis chosen once for the whole `graph` command -- see toolsTG.resolve_readtype().

    The combined,
    multi-comparison heatmap PDFs are saved at the top level; when --heatsubplots is set, the
    individual per-comparison heatmaps are also saved to an `individual/` subfolder, mirroring
    plotsCoverage.py's split-vs-combined output layout.
    '''
    grp = toolsTG.resolve_grp_column(adata, grp, 'group')

    # A non-default orientation is a different PICTURE of the same data, so its files carry a
    # token and can sit beside the default's rather than overwriting them. The CSV written
    # alongside is the DATA, so it takes no token and is never transposed -- two runs of one
    # analysis must still produce exports that diff against each other.
    suffix = '' if orientation == 'vertical' else f'_{orientation}'
    individual_output = f'{output}individual/'
    if heatsubplots:
        logger.info(toolsTG.builder(individual_output))

    # Create an empty df to store the log2FC values for each group so a combined heatmap can be generated as well
    df_combine = pd.DataFrame()
    # Create a heatmap for each group
    for readtype in readtypes:
        # `readtypes` arrives as resolved obs column names -- adataGraph.resolved_diffrts()
        # applies the command-wide read basis once, so no module decides its own denominator.
        # Create a color palette for the heatmap
        cmap = plotsPalette.gradient(settings, 'lfc')
        # The significance panel sits beside the log2FC panel in the same figure, so it
        # is resolved here and threaded through rather than read inside heatmap_plot.
        sig_cmap = plotsPalette.gradient(settings, 'significance')
        # Create a correlation matrix from reads stored in adata observations
        df, log2fc_dict = toolsTG.adataLog2FC(adata, grp, readtype, readcount_cutoff=cutoff, config_name=config_name, overwrite=overwrite).main()
        if df.empty:
            if threaded:
                threaded += f'No data for {readtype} heatmap.\n'
            else:
                logger.warning(f'No data for {readtype} heatmap.')
            continue
        df['readtype'] = readtype
        # combine df with df_combine by stacking them vertically if readtype is not total_unique or total
        if readtype not in ('nreads_total_unique_norm', 'nreads_total_norm'):
            df_combine = pd.concat([df_combine, df], axis=0)
        # save df to csv
        csv_output = output.replace('graphs', 'results')
        os.makedirs(csv_output, exist_ok=True)
        df.to_csv(f'{csv_output}{grp}_{readtype}_{cutoff}_{heatbound}_heatmap.csv')
        # Create a pdf with a heatmap for sorted by each group on each page
        with PdfPages(f'{output}{grp}_{readtype}_{cutoff}_{heatbound}_heatmap{suffix}.pdf') as pdf:
            for col in [i for i in df.columns.tolist() if 'log2' in i]:
                plt = heatmap_plot(df, col, cmap, sig_cmap, heatbound, settings, orientation=orientation)
                # Save figure
                if threaded:
                    threaded += f'Saving heatmap for {readtype} {col}...\n'
                else:
                    logger.info(f'Saving heatmap for {readtype} {col}...')
                if heatsubplots:
                    plt.savefig(f'{individual_output}{grp}_{readtype}_{cutoff}_{heatbound}_{col}_heatmap{suffix}.pdf', bbox_inches='tight')
                pdf.savefig(bbox_inches='tight')
                plt.close()
    # Create a heatmap for the combined groups
    if not df_combine.empty:
        with PdfPages(f'{output}{grp}_combine_{cutoff}_{heatbound}_heatmap{suffix}.pdf') as pdf:
            for col in [i for i in df_combine.columns.tolist() if 'log2' in i]:
                plt = heatmap_plot(df_combine, col, cmap, sig_cmap, heatbound, settings, orientation=orientation)
                # Save figure
                if threaded:
                    threaded += f'Saving heatmap for combine {col}...\n'
                else:
                    logger.info(f'Saving heatmap for combine {col}...')
                if heatsubplots:
                    plt.savefig(f'{individual_output}{grp}_combine_{cutoff}_{heatbound}_{col}_heatmap{suffix}.pdf', bbox_inches='tight')
                pdf.savefig(bbox_inches='tight')
                plt.close()
    if threaded:
        return threaded

def _cbar_kws(horizontal):
    '''
    Colorbar placement for one panel. Vertical, at the right, in both layouts.

    A horizontal colorbar under each stacked panel was tried first, on the theory that a
    vertical one would tower over a wide, short panel the way the correlation matrix's used
    to. Rendering it showed the opposite: a right-hand colorbar takes the panel's own height,
    so beside a short panel it is short, while the one under the BOTTOM panel landed on top of
    the rotated feature labels. Only the padding differs.
    '''
    return {'fraction': 0.05, 'pad': 0.02 if horizontal else 0.05}


def _horizontal_figsize(log_tdf):
    '''
    Figure size for the stacked layout, derived from what is being drawn.

    `square=True` ties each panel's height to width * rows/columns, so with many features and
    few comparisons the panels are very short. A fixed height then leaves a band of white
    between them that no layout engine removes -- the same class of problem the suptitle
    anchor below exists for. Height is therefore the two panels plus room for the rotated
    feature labels, estimated from the longest one.

    Underestimating is NOT safe here: the constrained layout keeps decorations inside the
    figure rather than letting them overflow, so a figure too short collapses the panels to
    slivers instead of spilling the labels out for `bbox_inches='tight'` to pick up.
    '''
    n_rows, n_cols = log_tdf.shape
    width = 12.0
    panel_height = (width - 2.0) * n_rows / max(n_cols, 1)
    label_height = 0.085 * max((len(str(c)) for c in log_tdf.columns), default=0)
    return (width, round(2 * panel_height + label_height + 1.0, 2))


def _readtype_labels(tdf, feature_labels, horizontal=False):
    '''
    Feature labels marked with their read type's glyph, or None to leave them alone.

    The glyph sits on the side NEAREST the cells it describes, which is the start of the
    string when the labels run horizontally and the end of it when they are rotated 90
    degrees -- rotation puts the first character at the bottom, furthest from the panel.

    Only the COMBINED heatmap mixes read types; a per-read-type heatmap would get the same
    glyph on every row, which teaches a reader nothing. That is decided from the data rather
    than from a flag: one distinct read type means nothing to disambiguate.

    Read types are taken POSITIONALLY. The combined frame stacks each read type under the same
    tRNA index labels, so `tdf.loc[label, 'readtype']` returns a Series rather than a scalar
    exactly when the marker is needed -- the defect that sank an earlier attempt at this.
    '''
    if 'readtype' not in tdf.columns:
        return None
    readtypes = list(tdf['readtype'])
    if len(set(readtypes)) < 2:
        return None
    glyphs = [plotsPalette.readtype_marker(readtype)['glyph'] for readtype in readtypes]
    if horizontal:
        return [f'{label} {glyph}' for glyph, label in zip(glyphs, feature_labels)]
    return [f'{glyph} {label}' for glyph, label in zip(glyphs, feature_labels)]


def _readtype_legend(fig, tdf, horizontal=False):
    '''
    Key the glyphs, using the same wording the published figure's legend uses.

    Entries cover only the read types actually drawn: naming the whole vocabulary would
    advertise read types this plot does not contain. Built from matplotlib marker codes rather
    than from the unicode glyphs, which is why the two are stored side by side in
    plotsPalette -- a legend drawn from one and labels from the other would drift apart.
    '''
    from matplotlib.lines import Line2D

    # Ordered by the vocabulary, not by which row happened to sort first: the same three read
    # types would otherwise be listed in a different order on the next plot.
    present = {plotsPalette.readtype_marker(rt)['label'] for rt in tdf['readtype']}
    ordered = [m for m in plotsPalette.READTYPE_MARKERS.values() if m['label'] in present]
    ordered += [plotsPalette.readtype_marker(rt) for rt in tdf['readtype']
                if plotsPalette.readtype_marker(rt)['label'] not in
                {m['label'] for m in plotsPalette.READTYPE_MARKERS.values()}]

    seen, handles, labels = set(), [], []
    for marker in ordered:
        if marker['label'] in seen:
            continue
        seen.add(marker['label'])
        handles.append(Line2D([], [], linestyle='none', marker=marker['marker'],
                              color=plotsPalette.AXIS_TEXT,
                              markerfacecolor=None if marker['filled'] else 'none'))
        labels.append(marker['label'])
    # Anchored clear of the axis label, which in the stacked layout sits at the bottom centre
    # and was being written through by the legend.
    offset = -0.12 if horizontal else -0.02
    fig.legend(handles, labels, loc='lower center', ncol=len(labels), frameon=False,
               bbox_to_anchor=(0.5, offset), fontsize=8)


def _readtype_subject(tdf):
    '''
    What this heatmap holds, for its title, or None when the frame does not say.

    Every heatmap used to carry the same title, so the file holding `total` counts and the
    file stacking four read types were indistinguishable once opened -- and their filenames
    are no help, since `total` and `combine` both read as "everything" while meaning opposite
    things (`total` is ONE read type, numerically the sum of the other four; `combine` stacks
    those four as separate rows). Renaming the files belongs to the terminology pass on the
    roadmap; saying what a figure holds does not have to wait for it.

    Derived from the same column the markers are, so the title and the glyphs cannot describe
    different things.
    '''
    if 'readtype' not in tdf.columns:
        return None
    labels, seen = [], set()
    for readtype in tdf['readtype']:
        label = plotsPalette.readtype_marker(readtype)['label']
        if label not in seen:
            seen.add(label)
            labels.append(label)
    if len(labels) == 1:
        return labels[0]
    ordered = [m['label'] for m in plotsPalette.READTYPE_MARKERS.values() if m['label'] in seen]
    ordered += [label for label in labels if label not in ordered]
    return f'Read types combined: {", ".join(ordered)}'


#: Accepted --heatorient values. 'vertical' is the historical layout: features on rows,
#: comparisons on columns, the two panels side by side.
HEAT_ORIENT_CHOICES = ('vertical', 'horizontal')


def heatmap_plot(df, col, cmap, sig_cmap, heatbound, settings=None, orientation='vertical'):
    # Sort df by the sum of the log2FC values for each row
    tdf = df.sort_values(by=col, ascending=False)
    # subset df to only include the top and bottom heatbound values only if the heatmap is larger than heatbound
    if len(tdf) > int(heatbound):
        tdf = pd.concat([tdf.iloc[:int(heatbound), :], tdf.iloc[-int(heatbound):, :]])  
    if orientation not in HEAT_ORIENT_CHOICES:
        raise toolsTG.InvalidParameterError(
            f"--heatorient '{orientation}' is not one of {', '.join(HEAT_ORIENT_CHOICES)}."
        )
    horizontal = orientation == 'horizontal'

    log_tdf = tdf[[i for i in tdf.columns if 'log2' in i]]
    pval_tdf = -np.log10(tdf[[i for i in tdf.columns if 'pval' in i]])
    if horizontal:
        log_tdf, pval_tdf = log_tdf.T, pval_tdf.T

    # Create a heatmap
    # Width follows the comparison count; the height is a default the style file may replace.
    default_size = _horizontal_figsize(log_tdf) if horizontal else (
        (16, 12) if len(tdf.columns) >= 20 else (6, 12))
    # The stacked layout is packed by matplotlib: its panels are short and the space between
    # them is exactly what a fixed layout leaves behind.
    layout = {'layout': 'constrained'} if horizontal else {}
    fig, axs = plt.subplots(*((2, 1) if horizontal else (1, 2)),
                            figsize=toolsTG.figsize_for(settings, default_size), **layout)
    # The shared axis is labelled once, on the panel at the outside edge of the stack: the
    # left panel in the side-by-side layout, the bottom panel when stacked.
    shared_axis_off = {'yticklabels': False} if not horizontal else {'xticklabels': False}
    # The read-type glyph rides IN the tick label rather than being drawn beside it, so it
    # moves with the label under either orientation and cannot fall out of alignment.
    marked = _readtype_labels(tdf, log_tdf.columns if horizontal else log_tdf.index,
                              horizontal=horizontal)
    labelled = {}
    if marked is not None:
        labelled = {'xticklabels': marked} if horizontal else {'yticklabels': marked}
    sns.heatmap(log_tdf, ax=axs[0], cmap=cmap, center=0, vmax=4, vmin=-4, cbar=True, square=True,
                cbar_kws=_cbar_kws(horizontal),
                **(shared_axis_off if horizontal else labelled))
    sns.heatmap(pval_tdf, ax=axs[1], cmap=sig_cmap, vmin=0, vmax=SIGNIFICANCE_VMAX, cbar=True, square=True,
                cbar_kws=_cbar_kws(horizontal),
                **(labelled if horizontal else shared_axis_off))
    # The feature names are rotated on whichever axis carries them; transposing moves them
    # from y to x, and they are long enough to need rotating in either position.
    for ax in axs:
        ax.tick_params(axis='x', labelrotation=90)
    if horizontal:
        axs[0].set_xlabel('')
    else:
        axs[1].set_ylabel('')


    cbar = axs[0].collections[0].colorbar
    # A right-hand colorbar takes its panel's height, and a stacked panel is short -- nine
    # tick labels on it overlap into an unreadable smear, so the scale is labelled at its
    # ends and centre instead. The side-by-side layout has the height for all nine.
    if horizontal:
        cbar.set_ticks([-4, 0, 4])
        cbar.set_ticklabels(['<=-4', '0', '>=4'])
    else:
        cbar.set_ticks([-4,-3,-2,-1,0,1,2,3,4])
        cbar.set_ticklabels(['<=-4','-3','-2','-1','0','1','2','3','>=4'])
    cbar.ax.tick_params(labelsize=8) 
    cbar = axs[1].collections[0].colorbar
    sig_ticks, sig_labels = significance_colorbar_ticks()
    cbar.set_ticks(sig_ticks)
    cbar.set_ticklabels(sig_labels)
    cbar.ax.tick_params(labelsize=8)

    # Panel titles first: the figure title is positioned relative to them below.
    axs[0].set_title('log2FC', fontsize=10)
    axs[1].set_title('pval', fontsize=10)
    # Anchor the figure title to where the panels ACTUALLY end, not to a fixed height.
    # `square=True` sizes the axes from the row count, so in a fixed-height figure a short
    # heatmap leaves a tall band of interior white that bbox_inches='tight' cannot trim,
    # because the title pins the top edge. Measured before this change: the gap ran to 44% of
    # the saved image at 4 rows against 5% at 20 -- exactly the row-count dependence this
    # removes. get_position() already reflects the aspect adjustment; the tight bbox is
    # preferred where a renderer is available because it also covers the panel titles.
    if marked is not None:
        _readtype_legend(fig, tdf, horizontal=horizontal)
    subject = _readtype_subject(tdf)
    title = 'Heatmap of log2FC of tRNA read counts between groups'
    title += f'\n{subject}, sorted by {col}' if subject else f' \nsorted by {col}'
    if horizontal:
        # The stacked layout has no such gap to close -- its height is derived from the panels
        # themselves -- and the constrained layout already reserves room for a suptitle. Naming
        # a y here would fight it, so the title is left where matplotlib puts it.
        plt.suptitle(title, fontsize=12)
        return plt
    try:
        renderer = fig.canvas.get_renderer()
        inv = fig.transFigure.inverted()
        top = max(a.get_tightbbox(renderer).transformed(inv).y1 for a in axs)
    except AttributeError:
        # Backends without get_renderer: fall back to the axes box plus room for the titles.
        top = max(a.get_position().y1 for a in axs) + 0.03
    plt.suptitle(title, fontsize=12, y=min(top + 0.01, 0.999), va='bottom')

    return plt

if __name__ == '__main__':
    pass
