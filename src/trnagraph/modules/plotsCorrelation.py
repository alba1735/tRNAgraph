#!/usr/bin/env python3

import logging

import pandas as pd
import anndata as ad

from . import toolsTG
from . import plotsPalette

import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns

logger = logging.getLogger(__name__)

def _plot_corr_matrix(df_wide, corr_method, corr_group, output, filename, title, threaded, settings=None):
    '''
    Compute and plot a correlation matrix (R^2 heatmap) from a wide (feature x sample/group)
    dataframe. Only plots when more than 20 samples/groups' worth of signal is present, matching
    the existing tRNA-only threshold.
    '''
    if df_wide.empty or df_wide.max().max() < 20:
        msg = f'Not enough samples to generate correlation matrix for {title}'
        if threaded:
            threaded += msg + '\n'
        else:
            logger.warning(msg)
        return threaded

    msg = f'Generating correlation matrix for {title}'
    if threaded:
        threaded += msg + '\n'
    else:
        logger.info(msg)

    df_corr = df_wide.corr(method=corr_method)
    plt.figure(figsize=(6, 6))
    ax = sns.heatmap(df_corr**2, square=True, cmap=plotsPalette.gradient(settings, 'correlation'), cbar_kws={'label': f'{corr_method} R^2'})
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title(f'{corr_method} {corr_group} {title} Correlation Matrix'.title())
    plt.gca().set_box_aspect(1)
    outpath = f'{output}{filename}.pdf'
    toolsTG.save_current(outpath, settings)
    msg = f'Correlation matrix for {title} saved to {outpath}'
    if threaded:
        threaded += msg + '\n'
    else:
        logger.info(msg)
    plt.close()
    return threaded

def _collapse_to_corr_group(df_wide, sample_to_group, corr_group):
    '''
    Collapse a wide (feature x sample) dataframe's columns down to corr_group level (mean),
    matching pivot_table's implicit aggregation on the tRNA side when corr_group != 'sample'.
    '''
    if corr_group == 'sample':
        return df_wide
    return df_wide.rename(columns=lambda c: sample_to_group.get(str(c), c)).T.groupby(level=0, observed=True).mean().T

def visualizer(adata, corr_method, corr_group, output, threaded=True, is_full_variant=True, read_basis=toolsTG.READ_BASIS_UNIQUE, settings=None):
    '''
    Generate correlation graphs for each sample in an AnnData object.

    tRNA-driven matrices (per readtype) use tRNAgraph's default tRNA/tRX-controlled DESeq2
    normalization (adata.obs '..._norm' columns). The non-tRNA matrix and the combined
    tRNA + non-tRNA matrix instead use the all-feature-controlled DESeq2 normalization
    (adata.uns['nontRNA_counts'] and adata.uns['deseq2_sizefactors_allfeatures']), since
    tRNA-controlled size factors are not representative of non-tRNA library composition --
    same three-way pattern as plotsPca.py.
    '''

    # Previously every '_norm' obs column got a matrix, so each run emitted both bases with
    # nothing distinguishing them but a filename token. Correlation has no comparative
    # overview page, so it plots one basis, selected once for the command by --allreads.
    is_unique_basis = read_basis == toolsTG.READ_BASIS_UNIQUE
    readtype_columns = [
        i for i in adata.obs.columns
        if '_norm' in i and ('_unique_' in i) == is_unique_basis
    ]
    df = pd.DataFrame(adata.obs, columns=['trna', corr_group] + readtype_columns)

    for i in readtype_columns:
        df_corr = df.pivot_table(index='trna', columns=corr_group, values=i, observed=True)
        title = i.split('_')[1]
        filename = f'{corr_method}_{corr_group}_{title}_correlation_matrix'
        threaded = _plot_corr_matrix(df_corr, corr_method, corr_group, output, filename, title, threaded, settings=settings)

    # Non-tRNA and combined tRNA + non-tRNA correlation matrices. These run automatically
    # alongside the per-readtype matrices above (no separate flag), but only when non-tRNA
    # feature counts are available -- adata.uns['nontRNA_counts'] is only populated with real
    # rows when `trnagraph analyze build` was given a --gtf; otherwise it's present but empty.
    nontrna_df, skip_message = toolsTG.resolve_nontrna_counts(adata, is_full_variant, 'correlation matrices')
    if nontrna_df is None:
        if threaded:
            threaded += skip_message + '\n'
        else:
            logger.warning(skip_message)
    else:
        sample_to_group = dict(zip(adata.obs['sample'].astype(str), adata.obs[corr_group].astype(str)))

        # Non-tRNA-only correlation. adata.uns['nontRNA_counts'] is normalized against the
        # all-feature-controlled DESeq2 size factors (see adataBuild.py), which is the
        # statistically appropriate normalization for non-tRNA feature counts.
        nontrna_grouped = _collapse_to_corr_group(nontrna_df, sample_to_group, corr_group)
        threaded = _plot_corr_matrix(nontrna_grouped, corr_method, corr_group, output,
                                      f'{corr_method}_{corr_group}_nontrna_correlation_matrix',
                                      'Non-tRNA RNAs', threaded, settings=settings)

        # Combined correlation: all tRNA reads (total, not unique-only) + non-tRNA RNAs, both
        # normalized against the all-feature-controlled size factors so the two feature sets
        # are on the same scale. This requires re-deriving tRNA total counts from raw counts,
        # since adata.obs only stores the tRNA/tRX-controlled (default) normalization.
        allfeature_sizefactors = adata.uns.get('deseq2_sizefactors_allfeatures')
        if allfeature_sizefactors is None or 'nreads_total_raw' not in adata.obs.columns:
            msg = ('No all-feature-controlled DESeq2 size factors found in AnnData object '
                   '(uns[\'deseq2_sizefactors_allfeatures\'] missing, or obs[\'nreads_total_raw\'] '
                   'unavailable). Skipping combined tRNA + non-tRNA correlation matrix. Rebuild with '
                   'the current `trnagraph analyze build` to enable this matrix.')
            if threaded:
                threaded += msg + '\n'
            else:
                logger.warning(msg)
        else:
            sf_map = {k: float(v) for k, v in dict(allfeature_sizefactors).items()}
            total_df = pd.DataFrame(adata.obs, columns=['trna', 'sample', 'nreads_total_raw'])
            missing_samples = sorted(set(total_df['sample'].astype(str)) - set(sf_map.keys()))
            if missing_samples:
                msg = f'All-feature size factors not found for samples {missing_samples}; defaulting to 1.0 for the combined correlation matrix.'
                if threaded:
                    threaded += msg + '\n'
                else:
                    logger.warning(msg)
            total_df['nreads_total_norm_allfeatures'] = total_df['nreads_total_raw'] / total_df['sample'].astype(str).map(lambda s: sf_map.get(s, 1.0))
            total_wide = total_df.pivot_table(index='trna', columns='sample', values='nreads_total_norm_allfeatures', observed=True)

            shared_cols = total_wide.columns.intersection(nontrna_df.columns)
            if len(shared_cols) < len(total_wide.columns) or len(shared_cols) < len(nontrna_df.columns):
                missing_from_nontrna = set(total_wide.columns) - set(nontrna_df.columns)
                missing_from_total = set(nontrna_df.columns) - set(total_wide.columns)
                msg = f'Sample columns did not fully match between tRNA and non-tRNA counts for the combined correlation matrix; using intersection of {len(shared_cols)} shared samples.'
                if missing_from_nontrna:
                    msg += f' Missing from nontRNA_counts: {sorted(missing_from_nontrna)}.'
                if missing_from_total:
                    msg += f' Missing from tRNA total counts: {sorted(missing_from_total)}.'
                if threaded:
                    threaded += msg + '\n'
                else:
                    logger.warning(msg)

            if len(shared_cols) == 0:
                msg = 'No overlapping sample columns between tRNA and non-tRNA counts. Skipping combined correlation matrix.'
                if threaded:
                    threaded += msg + '\n'
                else:
                    logger.warning(msg)
            else:
                combined_df = pd.concat([total_wide[shared_cols], nontrna_df[shared_cols]], axis=0)
                combined_grouped = _collapse_to_corr_group(combined_df, sample_to_group, corr_group)
                threaded = _plot_corr_matrix(combined_grouped, corr_method, corr_group, output,
                                              f'{corr_method}_{corr_group}_allrna_correlation_matrix',
                                              'All tRNA + Non-tRNA RNAs', threaded, settings=settings)

    if threaded:
        return threaded


if __name__ == '__main__':
    pass
