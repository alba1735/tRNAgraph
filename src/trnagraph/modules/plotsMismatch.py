#!/usr/bin/env python3

import contextlib
import logging

import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
import matplotlib.ticker as mticker
from matplotlib.backends.backend_pdf import PdfPages

from . import toolsTG
from . import plotsPalette

plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

logger = logging.getLogger(__name__)

# tRAX's boxplotmismatches.R:110 requires more than 50 reads of base composition behind a
# position before it will plot that position's rate, which keeps near-zero-coverage
# positions from producing spurious 100% misincorporation calls.
MIN_BASE_COMPOSITION = 50

# The dashed reference line tRAX draws across the position plots (boxplotmismatches.R:134).
# A visual convention rather than a filter -- nothing is dropped at this value.
REFERENCE_RATE = 0.1

BASE_COVERAGE_TYPES = ('adenines', 'thymines', 'cytosines', 'guanines')

# seaborn's stripplot jitter draws from numpy's global RNG, so an unseeded run produced a
# visibly different figure every time the same data was plotted -- fine for a glance, wrong
# for a figure that goes into a paper and may be regenerated during review. Seeded here, with
# the caller's RNG state restored afterwards so this never perturbs anything else.
JITTER_SEED = 0


@contextlib.contextmanager
def _seeded_jitter(seed=JITTER_SEED):
    state = np.random.get_state()
    np.random.seed(seed)
    try:
        yield
    finally:
        np.random.set_state(state)


def position_mismatch_rates(adata, pseudocount, numerator='mismatchedbases',
                            min_bases=MIN_BASE_COMPOSITION):
    '''
    Per-position misincorporation rate, one row per (feature, position).

    The rate is `numerator / (coverage + pseudocount)` computed from RAW counts and
    maximized across samples, matching tRAX's `FUN=max` aggregation.

    tRAX computes the same quantity from the size-factor-normalized coverage table while
    gating on raw base composition. Because the size factor cancels in `m/c` but not in
    `m/(c + pseudocount)`, that makes the pseudocount worth `pseudocount x sizefactor` in
    raw terms -- roughly 4 to 41 across the hg38 dataset -- and, since the aggregation is a
    maximum, systematically favours the lowest-size-factor sample. Reading raw counts here
    gives the pseudocount one meaning in every sample. The ratio itself is unaffected: the
    size factor divides out of numerator and denominator alike.
    '''
    var = adata.var
    nongap = ~var['gap'].astype(bool)
    raw = adata.layers['raw'] if 'raw' in adata.layers else adata.X

    def matrix(coverage_type):
        mask = (nongap & (var['coverage'] == coverage_type)).values
        values = np.asarray(raw[:, mask], dtype=float)
        return pd.DataFrame(values, index=adata.obs_names,
                            columns=var.loc[mask, 'positions'].tolist())

    coverage = matrix('coverage')
    counts = matrix(numerator)
    base_composition = sum(matrix(base) for base in BASE_COVERAGE_TYPES)

    rates = counts / (coverage + pseudocount)
    rates = rates.where(base_composition > min_bases)

    positions = list(coverage.columns)
    rates = rates.assign(trna=adata.obs['trna'].values)
    tidy = rates.melt(id_vars='trna', var_name='position', value_name='rate').dropna(subset=['rate'])

    # Collapse samples: one row per (feature, position), carrying that feature's worst sample.
    out = tidy.groupby(['trna', 'position'], observed=True)['rate'].max().reset_index()

    labels = adata.obs[['trna', 'amino', 'iso']].drop_duplicates('trna').set_index('trna')
    out['amino'] = out['trna'].map(labels['amino'])
    out['iso'] = out['trna'].map(labels['iso'])

    # var is already in Sprinzl order, so it -- not sorted() -- is the ordering authority.
    out['position'] = pd.Categorical(out['position'], categories=positions, ordered=True)
    return out.sort_values(['trna', 'position']).reset_index(drop=True)


def mismatch_count_shares(adata):
    '''
    Read-level mismatch histogram as per-sample shares, from uns['mismatch_counts'].

    Returns None when the key is absent, which is what any object built before that key
    existed looks like. The caller skips the plot rather than failing the graph run; the
    fix for a user who wants it is a rebuild, not a flag.
    '''
    stored = adata.uns.get('mismatch_counts')
    if stored is None:
        return None

    df = pd.DataFrame(stored).copy()
    samples = [col for col in df.columns if col not in ('count', 'type')]
    long = df.melt(id_vars=['count', 'type'], value_vars=samples,
                   var_name='Sample', value_name='reads')

    totals = long.groupby(['type', 'Sample'], observed=True)['reads'].transform('sum')
    long['share'] = np.where(totals > 0, long['reads'] / totals.where(totals > 0), 0.0)
    return long


def _palette(values, colormap, key, settings=None):
    '''Resolve a palette for `values`, preferring the user's colormap block for `key`.'''
    supplied = (colormap or {}).get(key, {})
    supplied = {k: (v if not str(v).startswith('#') else mplcolors.to_rgb(v))
                for k, v in supplied.items()}
    missing = [v for v in values if v not in supplied]
    if missing:
        generated = plotsPalette.categorical(settings, len(missing))
        supplied.update(dict(zip(missing, generated)))
    return supplied


class visualizer():
    '''
    Generate mismatch/misincorporation plots from an AnnData object.

    Ports tRAX's boxplotmismatches.R (per-position misincorporation and deletion) and
    plotreadmismatches.R (read-level mismatch histogram). The 5'/3' position-start panels
    boxplotmismatches.R also draws are deliberately not reproduced: plotsCoverage already
    renders readstarts/readends from the same object.
    '''
    def __init__(self, adata, colormap, output, pseudocount, threaded=True, settings=None):
        self.adata = adata
        self.colormap = colormap or {}
        self.output = output
        self.pseudocount = pseudocount
        self.threaded = threaded
        self.individual_output = f'{output}individual/'
        self.settings = settings

    def generate_plots(self):
        messages = []

        for numerator, label, stem in (
            ('mismatchedbases', 'Misincorporation', 'positionmismatches'),
            ('deletedbases', 'Deletion', 'positiondeletions'),
        ):
            df = position_mismatch_rates(self.adata, self.pseudocount, numerator=numerator)
            if df.empty:
                messages.append(f'No positions cleared the base-composition threshold for {stem}.\n')
                continue
            plt.close(self._plot_positions(df, hue='amino', palette_key='amino',
                                           title=f'Position {label}',
                                           outfile=f'{self.output}{stem}.pdf'))
            if numerator == 'mismatchedbases':
                self._plot_per_amino(df, label, stem)

        shares = mismatch_count_shares(self.adata)
        if shares is None:
            message = ('No mismatch_counts in this AnnData object -- skipping the read-level '
                       'mismatch histogram. Rebuild with `analyze build` to add it.\n')
            messages.append(message)
            logger.warning(message.strip())
        else:
            self._plot_count_histogram(shares)

        if self.threaded:
            return (self.threaded or '') + ''.join(messages)
        for message in messages:
            logger.info(message.strip())

    def _plot_positions(self, df, hue, palette_key, title, outfile=None, ylabel=None):
        '''
        Draw one position plot. Returns the figure; saves it too when `outfile` is given, so
        a caller can both write the per-amino file and page the same figure into a combined
        PDF without rendering it twice.
        '''
        df = df.copy()
        # seaborn takes its hue levels from a categorical's full category list, not from the
        # values actually present, so a per-amino subset would demand a palette covering every
        # anticodon in the dataset. Drop the unused categories first.
        if isinstance(df[hue].dtype, pd.CategoricalDtype):
            df[hue] = df[hue].cat.remove_unused_categories()
            levels = list(df[hue].cat.categories)
        else:
            levels = sorted(df[hue].dropna().unique())

        palette = _palette(levels, self.colormap, palette_key, self.settings)
        settings = self.settings or {}
        # stripplot's `size` is a diameter in points where scatter's `s` is an area, so the
        # shared marker_size is square-rooted to keep one number meaning the same visual
        # weight across plot types.
        marker_size = (settings['marker_size'] ** 0.5) if settings.get('marker_size') else 4
        rasterize = bool(settings.get('rasterize_over')) and len(df) > settings['rasterize_over']
        fig, ax = plt.subplots(figsize=toolsTG.figsize_for(self.settings, (20, 7)))
        with _seeded_jitter():
            sns.stripplot(data=df, x='position', y='rate', hue=hue, palette=palette,
                          size=marker_size, jitter=0.25, ax=ax, legend='brief',
                          rasterized=rasterize)
        ax.axhline(REFERENCE_RATE, linestyle='--', linewidth=1, color=plotsPalette.REFERENCE_LINE)
        ax.set_ylim(0, 1)
        # Rates are fractions; tRAX labels the same axis with scales::percent_format().
        ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1.0))
        ax.set_xlabel('Sprinzl position')
        ax.set_ylabel(ylabel or f'Maximum percent {title.split()[-1].lower()}')
        ax.set_title(title)
        ax.tick_params(axis='x', rotation=90)
        ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0, frameon=False)
        if outfile:
            outfile = toolsTG.save_current(outfile, self.settings)
            logger.info(f'Mismatch plot saved to {outfile}')
        return fig

    def _plot_per_amino(self, df, label, stem):
        # The multi-page roll-up belongs at the base alongside the overview, with the
        # single-amino files under individual/ -- the same split plotsCoverage uses
        # (combined_*.pdf at the category base, per-obs plots in a subfolder). Naming it
        # for what it holds, rather than `{stem}_combined`, keeps it from reading as a
        # duplicate of the top-level `{stem}.pdf` overview.
        logger.info(toolsTG.builder(self.individual_output))
        with PdfPages(f'{self.output}combined_{stem}_by_amino.pdf') as pdf:
            for amino in sorted(df['amino'].dropna().unique()):
                subset = df[df['amino'] == amino]
                if subset.empty:
                    continue
                fig = self._plot_positions(subset, hue='iso', palette_key='iso',
                                           title=f'Position {label} ({amino})',
                                           outfile=f'{self.individual_output}{amino}_{stem}.pdf')
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)

    def _plot_count_histogram(self, shares):
        outfile = f'{self.output}mismatchcounts.pdf'
        types = [t for t in ('trna', 'nontrna') if t in set(shares['type'])]
        fig, axes = plt.subplots(1, len(types), figsize=(6 * len(types), 6), squeeze=False)

        counts = sorted(shares['count'].unique())
        palette = _palette([str(c) for c in counts], self.colormap, 'mismatch', self.settings)

        for ax, readtype in zip(axes[0], types):
            subset = shares[shares['type'] == readtype]
            wide = subset.pivot(index='Sample', columns='count', values='share').fillna(0)
            bottom = np.zeros(len(wide))
            for count in counts:
                if count not in wide.columns:
                    continue
                ax.bar(wide.index, wide[count].values, bottom=bottom,
                       color=palette[str(count)], edgecolor=plotsPalette.BAR_EDGE, linewidth=0.4,
                       label=str(count))
                bottom = bottom + wide[count].values
            ax.set_title('tRNA' if readtype == 'trna' else 'non-tRNA')
            ax.set_ylabel('Percentage of total reads')
            ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1.0))
            ax.set_xlabel('Sample')
            ax.set_ylim(0, 1)
            ax.tick_params(axis='x', rotation=90)

        axes[0][-1].legend(title='Read\nmismatches', bbox_to_anchor=(1.01, 1),
                           loc='upper left', borderaxespad=0, frameon=False)
        outfile = toolsTG.save_current(outfile, self.settings)
        plt.close(fig)
        logger.info(f'Mismatch count histogram saved to {outfile}')
