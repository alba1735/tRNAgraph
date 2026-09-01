#!/usr/bin/env python3

import logging
import warnings

import numpy as np
import pandas as pd
import anndata as ad

from . import toolsTG
from . import plotsPalette

from functools import partial
import multiprocessing
from multiprocessing import Pool

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns

logger = logging.getLogger(__name__)

# Figure-fraction gap between the top of the subplot grid and the specificity overview's
# legend -- just enough to clear the top row's axis titles without reopening a dead band.
LEGEND_TITLE_CLEARANCE = 0.018

#: Grid geometry for the specificity pages. Four tRNA rows keeps each panel at the size the
#: arm-band furniture and tick labels were tuned for; eight group columns is as many as fit
#: across a page while staying legible.
#: Below this ceiling a plot is filed under low_coverage/. Shared by the coverage and
#: specificity individuals so one tRNA sorts the same way in both.
LOW_COVERAGE_CEILING = 20

PARTITION_ROWS_PER_PAGE = 4
PARTITION_COLUMNS_PER_PAGE = 8

#: Above this many --covgrp values the grid is refused outright rather than capped. Plotting
#: the first 24 of 436 would be a figure that misrepresents itself as complete.
MAX_PARTITION_COLUMNS = 24


# Pool tasks that reference a BOUND METHOD pickle `self` for every task -- and `self` carries
# the whole AnnData. Measured on the hg38 object that is 549 MB per task, against 436 tasks for
# the per-tRNA split and 28 per combined PDF, so a run spends its time shovelling gigabytes
# through a pipe to do about 1.3 seconds of drawing per page. Under the `fork` start method a
# worker already shares the parent's memory copy-on-write, so the visualizer is handed over
# through this module global instead and nothing large is pickled at all.
#
# It has to be SET BEFORE the Pool is created, since that is when the fork happens.
_WORKER_VISUALIZER = None


def _share_with_workers(visualizer):
    """Publish the visualizer for pool workers to pick up; returns False if they cannot."""
    global _WORKER_VISUALIZER
    _WORKER_VISUALIZER = visualizer
    # Only `fork` gives a child the parent's globals. Under spawn/forkserver the worker would
    # find None, so the caller runs the work serially rather than silently failing.
    return multiprocessing.get_start_method(allow_none=False) == 'fork'


def _combine_page_worker(args):
    ulist, coverage_fill = args
    return _WORKER_VISUALIZER.generate_combine_page(ulist, coverage_fill)


def _split_single_worker(covobs):
    return _WORKER_VISUALIZER.generate_split_single(covobs)


class visualizer():
    '''
    Generate coverage plots for each sample in an AnnData object.
    '''
    def __init__(self, adata, threads, coverage_grp, coverage_obs, coverage_type, coverage_gap, coverage_method, colormap, output, phase_tracker=None, quiet=False, settings=None):
        self.logger = logging.getLogger(__name__)
        self.threads = threads
        self.coverage_obs = coverage_obs
        self.coverage_grp = toolsTG.resolve_grp_column(adata, coverage_grp, 'coveragegrp')
        self.coverage_type = coverage_type
        self.coverage_gap = coverage_gap
        self.coverage_method = coverage_method
        # Shared outer progress tracker (toolsTG.PhaseTracker), passed in from adataGraph.py so
        # generate_split()'s per-plot completions can advance it directly in lockstep -- optional
        # (None) so this class stays usable standalone, e.g. outside adataGraph.py's orchestration.
        self.phase_tracker = phase_tracker
        self.quiet = quiet
        # Clean AnnData object for plotting
        self.adata, self.readstarts, self.readends = self.clean_adata(adata)
        # Held before the coverage_type subset so the stacked specificity overview can read
        # all four partition categories at once; every other plot uses the subset self.adata.
        self.partition_source = adata
        # This class owns its own directory layout: `output` is the coverage base, the
        # per-category plots go under <base>/<covtype-alias>/, and the stacked overview --
        # which shows every category at once and so has no single category to file under --
        # stays at <base>/. See toolsTG.coverage_category_dir().
        self.category_dir = toolsTG.coverage_category_dir(self.coverage_type)
        if colormap != None: #and self.coverage_combine_all == False:
            self.coverage_pal = {k:v if v[0]!='#' else mplcolors.to_rgb(v) for k,v in colormap.items()}
        else:
            coverage_pal = plotsPalette.categorical(settings, len(self.adata.obs[self.coverage_grp].unique()))
            self.coverage_pal = dict(zip(sorted(self.adata.obs[self.coverage_grp].unique()), coverage_pal))
        self.output = output
        self.category_output = f'{output}{self.category_dir}/'
        # Specificity gets its own directory rather than sitting at the coverage base or under
        # a --covtype: it shows every category at once, so it belongs to none of them, and the
        # combined grid now sits beside the individual plots it summarises.
        self.specificity_output = f'{output}specificity/'
        self.settings = settings

    def build_output_dirs(self):
        '''
        Create the directories this class writes into. Called by adataGraph before plotting
        so the per-category tree is made in one place rather than by whichever method
        happens to run first.
        '''
        logger.info(toolsTG.builder(self.category_output))
        logger.info(toolsTG.builder(f'{self.category_output}{self.coverage_obs}/'))
        logger.info(toolsTG.builder(f'{self.category_output}{self.coverage_obs}/low_coverage/'))
        logger.info(toolsTG.builder(self.specificity_output))
        logger.info(toolsTG.builder(f'{self.specificity_output}{self.coverage_obs}/'))
        logger.info(toolsTG.builder(f'{self.specificity_output}{self.coverage_obs}/low_coverage/'))

    def clean_adata(self, adata):
        '''
        Clean AnnData object for plotting.
        '''
        # Subset gaps from AnnData variables
        adata = adata[:,np.isin(adata.var.gap, self.coverage_gap)]
        # Drop nan values from AnnData
        adata = adata[~adata.obs[self.coverage_grp].isna()]
        # Subset just the readstarts and readends from AnnData variables
        readstarts = adata[:,np.isin(adata.var.coverage, ['readstarts'])].copy()
        readends = adata[:,np.isin(adata.var.coverage, ['readends'])].copy()
        # Subset by coverage type from AnnData variables
        adata = adata[:,np.isin(adata.var.coverage, [self.coverage_type])]

        return adata, readstarts, readends
    
    def __coverage_transform__(self, df, singlecol=False):
        '''
        Transform coverage data for plotting.
        '''
        # A grouping column with exactly one sample per label per covobs (e.g. --covgrp sample,
        # the Parameter Fallback default) makes pandas' single-label column selection collapse
        # to a Series rather than a single-column DataFrame -- already exactly the "one value
        # per row" result singlecol is trying to produce, so there's nothing left to aggregate.
        if singlecol and isinstance(df, pd.Series):
            return df
        # If the coverage method is mean, median, max, or min, transform the df by using groupby on column names
        if self.coverage_method == 'mean':
            if singlecol:
                df = df.mean(axis=1)
            else:
                df = df.T.groupby(level=0, observed=False).mean().T
        elif self.coverage_method == 'median':
            if singlecol:
                df = df.median(axis=1)
            else:
                df = df.T.groupby(level=0, observed=False).median().T
        elif self.coverage_method == 'max':
            if singlecol:
                df = df.max(axis=1)
            else:
                df = df.T.groupby(level=0, observed=False).max().T
        elif self.coverage_method == 'min':
            if singlecol:
                df = df.min(axis=1)
            else:
                df = df.T.groupby(level=0, observed=False).min().T
        elif self.coverage_method == 'sum':
            if singlecol:
                df = df.sum(axis=1)
            else:
                df = df.T.groupby(level=0, observed=False).sum().T
        else:
            raise ValueError(f'Invalid coverage method: {self.coverage_method}')
        
        return df

    def _covobs_list(self):
        '''tRNAs to plot, ordered by anticodon then copy number -- shared by the combined
        and specificity-overview pages so both paginate identically.'''
        ulist = sorted(self.adata.obs[self.coverage_obs].unique())
        if self.coverage_obs == 'trna':
            ulist = sorted(ulist, key=lambda x: ('-'.join(x.split('-')[:-1]), int(x.split('-')[-1])))
        return ulist

    def _partition_groups(self):
        '''The --covgrp values to draw, in a stable order.'''
        return sorted(self.partition_source.obs[self.coverage_grp].dropna().unique())

    def _partition_frame(self):
        '''
        Coverage per position for each of tRAX's four read-specificity categories, PER GROUP:
        {covobs: {group: DataFrame(position x category)}}.

        The group axis is the point. This used to aggregate across every --covgrp value at
        once, so the plot ignored --covgrp entirely and a treated and an untreated sample were
        averaged into one trace with nothing on the figure saying so. Within a group the
        samples are still reduced with --covmethod, exactly as the coverage plots do.

        Built once, up front: these are aggregated position vectors, so even a human build
        (436 tRNAs x 4 categories x a handful of groups) is a few hundred thousand floats.
        '''
        views = {}
        for category in toolsTG.COVERAGE_PARTITION:
            sub = self.partition_source[:, np.isin(self.partition_source.var.coverage, [category])]
            sub = sub[:, np.isin(sub.var.gap, self.coverage_gap)]
            views[category] = sub[~sub.obs[self.coverage_grp].isna()]
        groups = self._partition_groups()
        frames = {}
        for covobs in self._covobs_list():
            per_group = {}
            for group in groups:
                columns = {}
                for category, view in views.items():
                    subset = view[(view.obs[self.coverage_obs] == covobs)
                                  & (view.obs[self.coverage_grp] == group)]
                    df = pd.DataFrame(subset.X.T, columns=subset.obs[self.coverage_grp].values)
                    columns[category] = self.__coverage_transform__(df, singlecol=True)
                per_group[group] = pd.DataFrame(columns)
            frames[covobs] = per_group
        return frames

    def generate_partition_split(self):
        '''
        One specificity plot per --covobs value per --covgrp value.

        The grid answers "how do these groups compare"; this answers "what does this one group
        actually look like", at full size. Sorted into low_coverage/ by the same ceiling the
        individual coverage plots use, so a tRNA that is quiet in both views sorts the same way
        in both places.
        '''
        frames = self._partition_frame()
        palette = plotsPalette.discrete_colors(
            plotsPalette.gradient(self.settings, 'ordered'), len(toolsTG.COVERAGE_PARTITION))
        labels = [toolsTG.COVERAGE_CATEGORY_LABELS[c] for c in toolsTG.COVERAGE_PARTITION]
        for covobs, per_group in frames.items():
            for group, frame in per_group.items():
                fig, ax = plt.subplots(figsize=toolsTG.figsize_for(self.settings, (6, 5.5)))
                self.generate_partition_plot(frame, ax, covobs, palette)
                ax.set_title(f'{covobs} read specificity ({group})')
                # Each file needs its own key: unlike a grid page, there is nothing else on
                # the canvas to name the four stacked categories.
                handles = [mpatches.Patch(color=color, label=label)
                           for color, label in zip(palette, labels)]
                ax.legend(handles=handles[::-1], loc='upper left', bbox_to_anchor=(1, 1),
                          borderaxespad=0, frameon=False, title='Read specificity')
                outstart = f'{self.specificity_output}{self.coverage_obs}/'
                if ax.get_ylim()[1] < LOW_COVERAGE_CEILING:
                    outstart += 'low_coverage/'
                toolsTG.save_current(
                    f'{outstart}{covobs}_specificity_by_{self.coverage_grp}_{group}_'
                    f'{self.coverage_method}.pdf', self.settings)
                plt.close(fig)

    def _partition_pages(self, frames):
        '''
        The specificity grid, as one figure per page: tRNAs down, --covgrp values across.

        Paginated in BOTH axes -- four tRNA rows and eight group columns at a time -- so a
        wide experiment continues onto further pages rather than shrinking panels until
        nothing is readable.
        '''
        groups = self._partition_groups()
        if len(groups) > MAX_PARTITION_COLUMNS:
            raise toolsTG.InvalidParameterError(
                f"--covgrp '{self.coverage_grp}' has {len(groups)} values, more than the "
                f"{MAX_PARTITION_COLUMNS} the specificity grid can show. Group by a column "
                f"with fewer values, or filter this one with --config."
            )
        palette = plotsPalette.discrete_colors(
            plotsPalette.gradient(self.settings, 'ordered'), len(toolsTG.COVERAGE_PARTITION))
        labels = [toolsTG.COVERAGE_CATEGORY_LABELS[c] for c in toolsTG.COVERAGE_PARTITION]
        covobs_list = self._covobs_list()
        row_pages = [covobs_list[i:i + PARTITION_ROWS_PER_PAGE]
                     for i in range(0, len(covobs_list), PARTITION_ROWS_PER_PAGE)]
        column_pages = [groups[i:i + PARTITION_COLUMNS_PER_PAGE]
                        for i in range(0, len(groups), PARTITION_COLUMNS_PER_PAGE)]
        return [self._partition_page(rows, columns, frames, groups, palette, labels)
                for rows in row_pages for columns in column_pages]

    def _partition_page(self, rows, columns, frames, all_groups, palette, labels):
        '''One page of the grid. Panels in a row share a y-axis so the groups are comparable.'''
        fig, axes = plt.subplots(len(rows), len(columns), squeeze=False,
                                 figsize=(3.0 * len(columns) + 1.5, 3.0 * len(rows) + 1.5))
        for r, covobs in enumerate(rows):
            # Scaled against EVERY group of this tRNA, not just the ones on this page: two
            # column-pages of one tRNA drawn to different scales would invite exactly the
            # false comparison the grid exists to make possible.
            row_top = max((frames[covobs][group].sum(axis=1).max() for group in all_groups),
                          default=0)
            for c, group in enumerate(columns):
                ax = axes[r][c]
                self.generate_partition_plot(frames[covobs][group], ax, covobs, palette,
                                             xaxis=r == len(rows) - 1, top=row_top)
                # The group names caption the columns once, along the top; the tRNA names
                # caption the rows once, down the left. Repeating either in every panel would
                # crowd out the plots.
                ax.set_title(group if r == 0 else '')
                ax.set_ylabel(covobs if c == 0 else '')
        # The row label replaces each panel's own y-axis label, so the units are stated once
        # for the page rather than lost.
        fig.supylabel('Normalized Readcounts')
        handles = [mpatches.Patch(color=color, label=label)
                   for color, label in zip(palette, labels)]
        fig.legend(handles=handles[::-1], loc='lower center', frameon=False,
                   ncol=len(labels), bbox_to_anchor=(0.5, -0.01))
        fig.tight_layout()
        return fig

    def generate_partition_overview(self):
        '''
        Stacked plot of how specifically each position's reads could be assigned, porting
        tRAX's newcoverageplots.R `makecovplot(allmultmelt, ...)` -- its most informative
        coverage artifact and one tRNAgraph had no equivalent for.

        The four categories are mutually exclusive and sum to total coverage, so stacking
        them shows at a glance what fraction of a tRNA's signal is transcript-specific
        versus merely isodecoder-, isotype- or amino-level assignable. Reading that from the
        per-category directories instead would mean opening four PDFs and comparing heights
        by eye.

        Styled as tRNAgraph's own plots rather than as ggplot2's: a sequential mako ramp
        (light = least specific, dark = most) rather than R's categorical hues, since the
        categories are an ordered specificity scale, not unordered groups.

        Laid out as a grid -- one row per --covobs value, one column per --covgrp value, with
        a row sharing a y-axis -- so the groups can be read against each other. It previously
        aggregated across every group with --covmethod, which meant the figure ignored --covgrp
        entirely and averaged a treated and an untreated sample into one trace. With a single
        group the grid degenerates to exactly that old figure, which is why it replaces it
        rather than joining it.

        Lives under specificity/, not at the coverage base or under a --covtype, because it
        shows every category at once and so belongs to none of them -- and it is identical
        under --allreads, which selects a category rather than changing this partition.
        '''
        present = set(self.partition_source.var['coverage'])
        missing = [c for c in toolsTG.COVERAGE_PARTITION if c not in present]
        if missing:
            logger.warning(
                f'Skipping the coverage specificity overview: {missing} not found in this '
                f'AnnData object. Rebuild with the current `trnagraph analyze build` to '
                f'populate the full read-specificity partition.'
            )
            return
        frames = self._partition_frame()
        outend = (f'combined_{self.coverage_obs}_specificity_by_'
                  f'{self.coverage_grp}_{self.coverage_method}.pdf')
        with PdfPages(f'{self.specificity_output}{outend}') as pdf:
            for page in self._partition_pages(frames):
                pdf.savefig(page, bbox_inches='tight')
                plt.close(page)
        logger.info(f'Coverage specificity overview saved to {self.specificity_output}{outend}')

    def generate_partition_plot(self, df, ax, covobs, palette, xaxis=True, top=None):
        '''
        One stacked specificity plot. Deliberately reuses generate_plot()'s axis furniture --
        the same D/A/T-arm shading, the 37/58 modification guides and the 0-75 position range
        -- so a specificity page reads as the same kind of figure as a coverage page.
        '''
        ax.stackplot(df.index, *[df[c].values for c in toolsTG.COVERAGE_PARTITION],
                     colors=palette, zorder=2, linewidth=0)
        ax.set_title(f'{covobs} read specificity')
        ax.set_ylabel('Normalized Readcounts')
        if top:
            # A caller drawing a row of groups passes one ceiling for all of them, so the
            # panels can be read against each other rather than each being scaled to itself.
            ax.set_ylim(0, top)
        ax.set_xlim(0, 75)
        ax.set_xticks([18.01, 35.01, 37, 57.01, 58])
        ax.set_xticklabels(['\nD-Arm', '\nA-Arm', '37', '\nT-Arm', '58'])
        if xaxis:
            ax.set_xlabel('Positions on tRNA')
        # Same identical-ylim guard as generate_plot(): a tRNA with no coverage at all leaves
        # the top of the range at 0, which matplotlib warns about but which is benign here.
        top = top or ax.get_ylim()[1]
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Attempting to set identical low and high ylims")
            ax.set_ylim(0, top)
        for i in [37, 58]:
            ax.plot([i, i], [0, top], linewidth=1, ls='--', color=plotsPalette.REFERENCE_LINE, zorder=3)
        for i in [[14, 21], [32, 38], [54, 60], [10, 25], [27, 43], [49, 65]]:
            ax.fill_between(i, [top, top], color=plotsPalette.REGION_SHADE, alpha=plotsPalette.REGION_SHADE_ALPHA, zorder=0)
        ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False)

    def generate_combine(self):
        '''
        Generate combined coverage plots for all tRNAs using multiprocessing.
        '''
        ulist = self._covobs_list()
        # Generate list of tRNAs to plot split by 16 for each page
        ulist = [ulist[i*16:(i+1)*16] for i in range((len(ulist)+15)//16)]
        # Use multiprocessing to generate plotsand return them as a list so they can be saved to a pdf in order
        # Generate plots with confidence intervals
        outend = f'combined_{self.coverage_obs}_{self.coverage_type}_by_{self.coverage_grp}_with_ci_{self.coverage_method}.pdf'
        shared = _share_with_workers(self)
        with PdfPages(f'{self.category_output}{outend}') as pdf:
            for page in self._render_pages(ulist, 'ci', shared):
                pdf.savefig(page, bbox_inches='tight')
                plt.close(page)
        # Generate plots with fill
        outend = f'combined_{self.coverage_obs}_{self.coverage_type}_by_{self.coverage_grp}_with_fill_{self.coverage_method}.pdf'
        with PdfPages(f'{self.category_output}{outend}') as pdf:
            for page in self._render_pages(ulist, 'fill', shared):
                pdf.savefig(page, bbox_inches='tight')
                plt.close(page)

    def _render_pages(self, ulist, coverage_fill, shared):
        '''Render the combined pages, in parallel where the workers can share this object.'''
        if not shared or self.threads <= 1:
            return [self.generate_combine_page(chunk, coverage_fill) for chunk in ulist]
        with Pool(self.threads) as p:
            return p.map(_combine_page_worker, [(chunk, coverage_fill) for chunk in ulist])

    def generate_combine_page(self, ulist, coverage_fill):
        '''
        Generate combined coverage plots page for PdfPages via multiprocessing.
        '''
        # Generate figure
        fig_pdf = plt.figure(figsize=(24,22))
        for i, covobs in enumerate(ulist):
            # Turn off xaxis for all but bottom row
            xaxis = True if ulist.index(covobs) > 11 else False
            # Turn off legend for all but right column
            lgnd = True if ulist.index(covobs) % 4 == 3 else False
            # Generate subplot
            ax = fig_pdf.add_subplot(4,4,i+1)
            df = pd.DataFrame(self.adata[self.adata.obs[self.coverage_obs] == covobs].X.T, columns=self.adata[self.adata.obs[self.coverage_obs] == covobs].obs[self.coverage_grp].values)
            self.generate_plot(df, ax, covobs, coverage_fill=coverage_fill, lgnd=lgnd, xaxis=xaxis)
        return fig_pdf

    def generate_split(self):
        ulist = self.adata.obs[self.coverage_obs].unique()
        # imap_unordered (rather than map) so progress_iterator gets a real per-plot completion
        # signal -- order doesn't matter here, each call only has the side effect of writing its
        # own PDF file(s). If a shared PhaseTracker was passed in (the outer graphing bar), tick
        # it once per completed plot so it moves in lockstep with this, the one graph type that
        # can produce hundreds/thousands of plots, rather than only ticking atomically once the
        # whole "coverage" phase finishes.
        shared = _share_with_workers(self)
        if not shared or self.threads <= 1:
            results = (self.generate_split_single(covobs) for covobs in ulist)
            for _ in toolsTG.progress_iterator(
                results, total=len(ulist), desc="Coverage plots", logger=self.logger,
                quiet=self.quiet,
            ):
                if self.phase_tracker is not None:
                    self.phase_tracker.advance(1)
            return
        with Pool(self.threads) as p:
            for _ in toolsTG.progress_iterator(
                p.imap_unordered(_split_single_worker, ulist),
                total=len(ulist), desc="Coverage plots", logger=self.logger, quiet=self.quiet,
            ):
                if self.phase_tracker is not None:
                    self.phase_tracker.advance(1)

    def generate_split_single(self, covobs):
        # Create df for single tRNA
        lgnd = True
        if self.coverage_grp == self.coverage_obs:
            lgnd = False
        df = pd.DataFrame(self.adata[self.adata.obs[self.coverage_obs] == covobs].X.T, columns=self.adata[self.adata.obs[self.coverage_obs] == covobs].obs[self.coverage_grp].values)
        low_coverage, outstart = self._assemble_plot_(df,covobs,lgnd,'ci')
        self._assemble_plot_(df,covobs,lgnd,'fill')
        # Generate plot with readstarts and readends
        rslist = self.adata.obs[self.coverage_grp].unique()
        if self.coverage_grp == self.coverage_obs:
            rslist = [covobs]
        for i in rslist:
            if not low_coverage:
                # Generate plot with readstarts and readends
                fig, ax = plt.subplots(figsize=toolsTG.figsize_for(self.settings, (6, 5.5)))
                # Subset the df to the current group
                df_grp = df[i]
                if df_grp.sum().sum() != 0:
                    # Get the readstarts/ends for the current group and get the method of each position
                    df_readstarts = pd.DataFrame(self.readstarts[self.readstarts.obs[self.coverage_obs] == covobs].X.T, columns=self.readstarts[self.readstarts.obs[self.coverage_obs] == covobs].obs[self.coverage_grp].values)[i]
                    df_readends = pd.DataFrame(self.readends[self.readends.obs[self.coverage_obs] == covobs].X.T, columns=self.readends[self.readends.obs[self.coverage_obs] == covobs].obs[self.coverage_grp].values)[i]
                    df_readstarts = self.__coverage_transform__(df_readstarts, singlecol=True)
                    df_readends = self.__coverage_transform__(df_readends, singlecol=True)
                    # Combine the readstarts and readends into a single df
                    df_readstartends = pd.concat([df_readstarts, df_readends], axis=1)
                    df_readstartends.columns = ['readstarts', 'readends']
                    # Scale the df so that sum of each column is 1
                    df_readstartends = df_readstartends/df_readstartends.sum()
                    # Add position column to df
                    df_readstartends['position'] = df_readstartends.index
                    # Melt the df so that readstarts and readends are in the same column for plotting
                    df_readstartends = df_readstartends.melt(id_vars='position', value_vars=['readstarts', 'readends'])
                    self.generate_plot(df_grp, ax, covobs, coverage_fill='endstarts', rse=df_readstartends, standalone=True)
                    # Save the plot
                    if self.coverage_grp == self.coverage_obs:
                        outend = f'{covobs}_{self.coverage_type}_with_endstarts_{self.coverage_method}.pdf'
                    else:
                        outend = f'{covobs}_{self.coverage_type}_by_{self.coverage_grp}_{i}_with_endstarts_{self.coverage_method}.pdf'
                    toolsTG.save_current(outstart+outend, self.settings)
            plt.close()

    def _assemble_plot_(self, df, covobs, lgnd, cov_fill):
        # Generate plot with confidence intervals
        fig, ax = plt.subplots(figsize=toolsTG.figsize_for(self.settings, (6, 5.5)))
        self.generate_plot(df, ax, covobs, coverage_fill=cov_fill, lgnd=lgnd, standalone=True)
        # Get max y value for plot
        outend = f'{covobs}_{self.coverage_type}_by_{self.coverage_grp}_with_{cov_fill}_{self.coverage_method}.pdf'
        if self.coverage_grp == self.coverage_obs:
            outend = f'{covobs}_{self.coverage_type}_with_{cov_fill}_{self.coverage_method}.pdf'
        if cov_fill == 'ci':
            outend = '_'.join(outend.split('_')[:-1]) + '.pdf'
        outstart = f'{self.category_output}{self.coverage_obs}/'
        low_coverage = False
        if ax.get_ylim()[1] < LOW_COVERAGE_CEILING:
            outstart += 'low_coverage/'
            low_coverage = True
        toolsTG.save_current(outstart+outend, self.settings)
        plt.close()

        return low_coverage, outstart

    def generate_plot(self, df, ax, trna, coverage_fill, lgnd=True, xaxis=True, rse=None, standalone=False):
        '''
        Generate coverage plots for a single tRNA.

        `standalone` says whether this figure is one of the per-tRNA files or a panel of a
        combined page. It governs one thing: whether a --style `line_width` reaches the trace.
        The combined and multi-page pages deliberately keep their tuned widths, exactly as
        they keep their computed figsize -- shrinking a panel grid is not something a style
        file does, so there is nothing there for a stroke width to stay in proportion with.
        '''
        # Only a standalone plot consults the style; a panel keeps the module's own width.
        trace_settings = self.settings if standalone else None
        # Get the cov method of all the columns
        if coverage_fill == 'ci': 
            sns.lineplot(data=df, linewidth=toolsTG.linewidth_for(trace_settings, 2), dashes=False, palette=self.coverage_pal, errorbar=('se',2), zorder=2, ax=ax)
        elif coverage_fill == 'fill':
            df = self.__coverage_transform__(df)
            sns.lineplot(data=df, linewidth=toolsTG.linewidth_for(trace_settings, 2), dashes=False, palette=self.coverage_pal, errorbar=('se',False), zorder=2, ax=ax)
        elif coverage_fill == 'endstarts':
            # Get the cov method of all the columns
            df = self.__coverage_transform__(df, singlecol=True)
            # Scale the df so that sum of all values is 1
            df = df/df.sum()
            sns.lineplot(data=df, linewidth=toolsTG.linewidth_for(trace_settings, 1.5), dashes=False, color=plotsPalette.COVERAGE_TRACE_MUTED, zorder=3, ax=ax)
            sns.histplot(data=rse, x='position', weights='value', hue='variable', palette=[plotsPalette.READSTART_COLOR, plotsPalette.READEND_COLOR], alpha=0.5, discrete=True, zorder=2, \
                         stat='probability', common_norm=False, linewidth=0, legend=True, ax=ax)
        # Fill the area under the curve with the mean by going in order of the column with the highest mean to the lowest if fill/both is specified
        if coverage_fill == 'fill':
            df_mean = df.T.groupby(level=0, observed=False).mean().T # This is the mean of the columns with the same name
            for i in df_mean.mean().sort_values(ascending=False).index:
                ax.fill_between(df_mean.index, df_mean[i], color=self.coverage_pal.get(i), alpha=0.35, zorder=1)
        # Set plot parameters
        plt.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=True, labelleft=True)
        ax.set_title(f'{trna} {self.coverage_type}')
        ax.set_ylabel("Normalized Readcounts")
        if coverage_fill == 'endstarts':
            ylim = [0,1]
            # Set as 0 to 100%
            ax.set_yticks([0,0.25,0.5,0.75,1])
            ax.set_yticklabels(['0%','25%','50%','75%','100%'])
        else:
            ylim = ax.get_ylim()
        # A tRNA with flat/all-zero coverage at this position range can leave ylim[1] at 0,
        # making the floor-and-ceiling re-clamp below a no-op that matplotlib otherwise warns
        # about ("Attempting to set identical low and high ylims") -- benign here, so scope the
        # suppression to just this call rather than silencing it for the whole plotting module.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Attempting to set identical low and high ylims")
            ax.set_ylim(0, ylim[1])
        # Add dashed lines to plot for common tRNA modifications
        if coverage_fill == 'endstarts':
            # Create a df with the same columns as the rse df but empty
            names_df = pd.DataFrame({'position_start':[18.01,35.01,57.01], 'position_end':[18.01,35.01,57.01], 'name':['\nD-Arm','\nA-Arm','\nT-Arm']},\
                                     columns=['position_start','position_end','name'])
            # Add vertical dashed lines to plot for histobars that are above 10%
            for rt in rse['variable'].unique():
                current_start, current_bound = 0, 0
                endstart_switch = False
                for i, row in rse[rse['variable']==rt].iterrows():
                    if row['value'] >= 0.1:
                        if row['position'] - current_start > 1:
                            plt.plot([row['position']-0.5,row['position']-0.5],[0,ylim[1]],linewidth=1,ls='--',color=plotsPalette.REFERENCE_LINE, zorder=3)
                            current_bound = row['position']
                        current_start = row['position']
                        endstart_switch = True
                    else:
                        if endstart_switch:
                            plt.plot([row['position']-0.5,row['position']-0.5],[0,ylim[1]],linewidth=1,ls='--',color=plotsPalette.REFERENCE_LINE, zorder=3)
                            # Set name as the avg of the start and end positions
                            if current_bound+1 == row['position']:
                                names_df = pd.concat([names_df, pd.DataFrame({'position_start':[current_bound-1], 'position_end':[row['position']], 'name':[str(current_bound)]},\
                                                                              columns=['position_start','position_end','name'])])
                            else:
                                names_df = pd.concat([names_df, pd.DataFrame({'position_start':[current_bound-1], 'position_end':row['position'], 'name':[str(current_bound+1)+'-'+str(row['position'])]},\
                                                                            columns=['position_start','position_end','name'])])
                            endstart_switch = False
                if endstart_switch:
                    if current_bound == 75:
                        names_df = pd.concat([names_df, pd.DataFrame({'position_start':[current_bound], 'position_end':[75], 'name':[str(current_bound)]},\
                                                                      columns=['position_start','position_end','name'])])
                    else:
                        names_df = pd.concat([names_df, pd.DataFrame({'position_start':[current_bound], 'position_end':[75], 'name':[str(current_bound+1)+'-'+str(75)]},\
                                                                      columns=['position_start','position_end','name'])])
        else:
            for i in [37, 58]:
                plt.plot([i,i],[0,ylim[1]],linewidth=1,ls='--',color=plotsPalette.REFERENCE_LINE, zorder=3)
        # Set xaxis parameters
        if xaxis:
            ax.set_xlabel("Positions on tRNA")
        if coverage_fill == 'endstarts':
            ax.set_xlim(-0.5, 75.5)
            # Add horizontal dashed line to plot for histo above 10%
            plt.plot([-0.5,75.5],[0.1,0.1],linewidth=1,ls='--',color=plotsPalette.REFERENCE_LINE, zorder=3)
            # Add the xticks and xticklabels
            names_df['avg'] = (names_df['position_end'] + names_df['position_start'])/2
            names_df = names_df.sort_values(by=['avg'])
            ax.set_xticks(names_df['avg'])
            ax.set_xticklabels(names_df['name'])
        else:
            ax.set_xlim(0, 75) # ax.set_xlim(0, 73)
            ax.set_xticks([18.01,35.01,37,57.01,58])
            ax.set_xticklabels(['\nD-Arm','\nA-Arm','37','\nT-Arm','58'])
        # Add coverage regions to background
        for i in [[14,21],[32,38],[54,60],[10,25],[27,43],[49,65]]:
            ax.fill_between(i,[ylim[1],ylim[1]], color=plotsPalette.REGION_SHADE, alpha=plotsPalette.REGION_SHADE_ALPHA, zorder=0)
        # Capatilize the legend and move the legend outside the plot and remove the border around it
        if lgnd:
            if coverage_fill == 'endstarts':
                classes = ['Readstarts', 'Readends']
                # Must come from the same constants the bars are drawn with; these swatches
                # used to be hardcoded separately, so changing the palette recoloured the bars
                # and left the legend claiming the old colors.
                class_colours = [plotsPalette.READSTART_COLOR, plotsPalette.READEND_COLOR]
                recs = []
                for i in range(0, len(class_colours)):
                    recs.append(mpatches.Rectangle((0, 0), 1, 1, fc=class_colours[i]))
                ax.legend(recs, classes, loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
                ax.legend_.set_title('Coverage')
            else:
                handles, labels = ax.get_legend_handles_labels()
                ax.legend(handles=handles, labels=[x.capitalize() for x in labels])
                ax.legend(loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
                ax.legend_.set_title(self.coverage_grp.capitalize())
        else:
            ax.legend_.remove()
        # Set the box aspect ratio to 1 so the plot is square
        plt.gca().set_box_aspect(1)

if __name__ == '__main__':
    pass