#!/usr/bin/env python3

import logging

import seaborn as sns
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

from . import toolsTG

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

logger = logging.getLogger(__name__)


def _generate_pca_plots(df, hue_dict, colormap, pcamarkers, pcacolors, output, basename, title_suffix, threaded):
    '''
    Fit PCA on a wide (feature x sample) dataframe and generate the explained
    variance ratio barplot, PCA scatter plot, and pairplot. `basename` controls
    the output filenames (`{output}{basename}_{evr,pca,pairplot}.pdf`);
    `title_suffix`, when truthy, is appended to each plot title in parentheses.
    '''
    suffix = f' ({title_suffix})' if title_suffix else ''

    # Scale the data
    df = pd.DataFrame(StandardScaler().fit_transform(df), columns=df.columns, index=df.index)
    # Create a PCA object
    pca = PCA(n_components=min(len(df.columns), 5))
    pca.fit_transform(df)
    evr = pca.explained_variance_ratio_
    # Print the explained variance ratio
    if threaded:
        threaded += f'Principal components: {[f"PC{x}" for x in range(1, len(evr)+1)]}\n'
        threaded += f'Explained variance: {[f"{i:.4f}" for i in pca.explained_variance_]}\n'
        threaded += f'Explained variance ratio: {[f"{i*100:.2f}%" for i in evr]}\n'
    else:
        logger.info('Principal components: {}'.format([f'PC{x}' for x in range(1, len(evr)+1)]))
        logger.info('Explained variance: {}'.format([f'{i:.4f}' for i in pca.explained_variance_]))
        logger.info('Explained variance ratio: {}'.format([f'{i*100:.2f}%' for i in evr]))
    # Transform the data and create a new dataframe
    pca_index = ['PC{}'.format(x) for x in range(1, len(evr)+1)]
    df_pca = pd.DataFrame(pca.components_, columns=df.columns, index=pca_index).T

    # Plot the explained variance ratio
    plt.figure(figsize=(6, 6))
    ax = sns.barplot(x=['PC{}'.format(x) for x in range(1, len(evr)+1)], y=evr, palette=sns.husl_palette(len(evr)), hue=['PC{}'.format(x) for x in range(1, len(evr)+1)])
    # Set the x and y labels and title
    ax.set_xlabel('Principal Component')
    ax.set_ylabel('Explained Variance Ratio')
    ax.set_title('Explained Variance Ratio of Principal Components' + suffix)
    # Set the box aspect ratio to 1 so the plot is square
    plt.gca().set_box_aspect(1)
    # Save the plot
    plt.savefig(f'{output}{basename}_evr.pdf', bbox_inches='tight')
    if threaded:
        threaded += f'Explained variance ratio graph saved to {output}{basename}_evr.pdf\n'
    else:
        logger.info(f'Explained variance ratio graph saved to {output}{basename}_evr.pdf')
    plt.close()

    # Plot the data with seaborn
    plt.figure(figsize=(8, 8))
    if colormap:
        ax = sns.scatterplot(data=df_pca, x='PC1', y='PC2', s=100, palette=colormap, hue=hue_dict, legend='full')
    else:
        ax = sns.scatterplot(data=df_pca, x='PC1', y='PC2', s=100, palette=sns.husl_palette(len(set(hue_dict.values()))), hue=hue_dict, legend='full')
    ax.set_xlabel('PC1 ({:.2f}%)'.format(evr[0]*100))
    ax.set_ylabel('PC2 ({:.2f}%)'.format(evr[1]*100))
    # Capatilize the legend and move the legend outside the plot and remove the border around it
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=[x.capitalize() for x in labels])
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False)
    ax.legend_.set_title(pcacolors.capitalize())
    # Give the plot a title
    ax.set_title(f'PCA of {pcamarkers} colored by {pcacolors}' + suffix)
    # Remove the ticks and tick labels
    ax.tick_params(axis='both', which='both', bottom=False, top=False,
                   left=False, right=False, labelbottom=False, labelleft=False)
    # Set the box aspect ratio to 1 so the plot is square
    plt.gca().set_box_aspect(1)
    # Save the plot
    plt.savefig(f'{output}{basename}_pca.pdf', bbox_inches='tight')
    if threaded:
        threaded += f'PCA graph saved to {output}{basename}_pca.pdf\n'
    else:
        logger.info(f'PCA graph saved to {output}{basename}_pca.pdf')
    plt.close()

    # Plot pairplot of the data with seaborn
    plt.figure(figsize=(10, 10))
    # Rename columns to add explained variance ratio to each PC as well as add colors for use in seaborn
    df_pca.columns = ['PC{} ({:.2f}%)'.format(i, evr[i-1]*100) for i in range(1, len(df_pca.columns)+1)]
    df_pca[pcacolors.capitalize()] = [hue_dict[x] for x in df_pca.index]
    # Plot the pairplot
    if colormap:
        ax = sns.pairplot(df_pca, hue=pcacolors.capitalize(), palette=colormap, hue_order=sorted(set(hue_dict.values())))
    else:
        ax = sns.pairplot(df_pca, hue=pcacolors.capitalize(), palette=sns.husl_palette(len(set(hue_dict.values()))), hue_order=sorted(set(hue_dict.values())))
    # Remove the ticks and tick labels
    ax.tick_params(axis='both', which='both', bottom=False, top=False,
                   left=False, right=False, labelbottom=False, labelleft=False)
    # Give the plot a title
    ax.fig.suptitle(f'PCA Pairplot of {pcamarkers} colored by {pcacolors}' + suffix, y=1.02)
    # Set the box aspect ratio to 1 so the plot is square
    plt.gca().set_box_aspect(1)
    plt.savefig(f'{output}{basename}_pairplot.pdf', bbox_inches='tight')
    if threaded:
        threaded += f'Pairplot graph saved to {output}{basename}_pairplot.pdf\n'
    else:
        logger.info(f'Pairplot graph saved to {output}{basename}_pairplot.pdf')
    plt.close()

    return threaded


# PCA's both-bases comparison used to exist only as a side effect of --pcareadtypes
# defaulting to ['total_unique', 'total']. Now that the basis has been lifted out of that
# flag and onto the command-wide --allreads, the comparison is pinned here instead, so it
# survives whatever --pcareadtypes asks for and whichever basis is active. Same rationale
# and same constant name as plotsVolcano.OVERVIEW_TRNA_READTYPES: a labelled side-by-side
# comparison is not the silent cross-plot inconsistency --allreads exists to remove.
OVERVIEW_TRNA_READTYPES = ('nreads_total_unique_norm', 'nreads_total_norm')


def _readtype_label(column):
    '''obs column name -> the token used in filenames/titles, e.g. 'total_unique'.'''
    return column.replace('nreads_', '').replace('_norm', '')


def visualizer(adata, pcamarkers, pcacolors, pcareadtypes, colormap, output, threaded=True, is_full_variant=True, read_basis=toolsTG.READ_BASIS_UNIQUE):
    '''
    Generate PCA visualizations for each sample in an AnnData object.

    tRNA-driven plots (per readtype, e.g. total/total_unique) use tRNAgraph's default
    tRNA/tRX-controlled DESeq2 normalization (adata.obs 'nreads_*_norm', matching adata.X).
    The non-tRNA plot and the combined tRNA + non-tRNA plot instead use the
    all-feature-controlled DESeq2 normalization (adata.uns['nontRNA_counts'] and
    adata.uns['deseq2_sizefactors_allfeatures']), since tRNA-controlled size factors are not
    representative of non-tRNA library composition.
    '''
    # if pcamarkers not in adata.obs.columns:
    #     raise ValueError('Specified pcamarkers not found in AnnData object.')
    # if pcacolors not in adata.obs.columns:
    #     raise ValueError('Specified pcacolor not found in AnnData object.')
    # --pcareadtypes carries bare readtypes; the basis comes from --allreads via read_basis.
    # 'all' means every readtype that exists in both bases (toolsTG.DUAL_BASIS_READTYPES) --
    # the previous literal list named 'whole_unique', which is not a column adataBuild has
    # ever written (it is 'wholecounts_unique'), so that entry silently produced nothing.
    if 'all' in pcareadtypes:
        pcareadtypes = list(toolsTG.DUAL_BASIS_READTYPES)
    readtype_columns = [toolsTG.resolve_readtype(rt, read_basis, adata) for rt in pcareadtypes]
    # Guarantee the both-bases comparison, appending only what was not already requested so a
    # readtype is never plotted twice into the same filename.
    for overview_column in OVERVIEW_TRNA_READTYPES:
        if overview_column not in readtype_columns:
            readtype_columns.append(overview_column)

    # Create dictionary of sample and pcamarkers parameter for use in seaborn, and validate the
    # colormap. This only depends on pcamarkers/pcacolors (not on readtype), so it's computed once
    # up front and reused for the per-readtype plots as well as the non-tRNA / combined plots below.
    if pcamarkers == pcacolors:
        meta_df = pd.DataFrame(adata.obs, columns=['trna', pcamarkers])
        hue_dict = dict(zip(meta_df[pcamarkers], meta_df[pcamarkers]))
    else:
        meta_df = pd.DataFrame(adata.obs, columns=['trna', pcamarkers, pcacolors])
        hue_dict = dict(zip(meta_df[pcamarkers], meta_df[pcacolors]))
    if colormap != None:
        colormap = {k:v if v[0]!='#' else mplcolors.to_rgb(v) for k,v in colormap.items()}
        for v in hue_dict.values():
            if v not in colormap:
                if threaded:
                    threaded += f'Color {v} not found in colormap. Using default colors instead.\n'
                else:
                    logger.warning(f'Color {v} not found in colormap. Using default colors instead.')
                colormap = None
                break

    for rt in readtype_columns:
        readtype = _readtype_label(rt)
        # Create a dataframe with trna, pcamarkers parameter, and nreads from adata
        if pcamarkers == pcacolors:
            df = pd.DataFrame(adata.obs, columns=['trna', pcamarkers, rt])
        else:
            df = pd.DataFrame(adata.obs, columns=['trna', pcamarkers, rt, pcacolors])
        # Pivot the dataframe to have trna as the index, sample as the columns, and nreads as the values for dimensionality reduction
        df = df.pivot_table(index='trna', columns='sample', values=rt, observed=True)
        threaded = _generate_pca_plots(df, hue_dict, colormap, pcamarkers, pcacolors, output, basename=f'tRNA_{pcamarkers}_by_{pcacolors}_{readtype}', title_suffix=readtype.replace('_', ' ').title(), threaded=threaded)

    # Non-tRNA and combined tRNA + non-tRNA PCA plots. These run automatically alongside the
    # per-readtype plots above (no separate flag), but only when non-tRNA feature counts are
    # available -- adata.uns['nontRNA_counts'] is only populated with real rows when
    # `trnagraph analyze build` was given a --gtf; otherwise it's present but empty.
    nontrna_df, skip_message = toolsTG.resolve_nontrna_counts(adata, is_full_variant, 'PCA plots')
    if nontrna_df is None:
        if threaded:
            threaded += skip_message + '\n'
        else:
            logger.warning(skip_message)
    else:
        # Non-tRNA-only PCA. adata.uns['nontRNA_counts'] is normalized against the
        # all-feature-controlled DESeq2 size factors (see adataBuild.py), which is the
        # statistically appropriate normalization for non-tRNA feature counts.
        threaded = _generate_pca_plots(nontrna_df.copy(), hue_dict, colormap, pcamarkers, pcacolors, output, basename=f'nontRNA_{pcamarkers}_by_{pcacolors}', title_suffix='Non-tRNA RNAs', threaded=threaded)

        # Combined PCA: all tRNA reads (total, not unique-only) + non-tRNA RNAs, both normalized
        # against the all-feature-controlled size factors so the two feature sets are on the
        # same scale. This requires re-deriving tRNA total counts from raw counts, since
        # adata.obs only stores the tRNA/tRX-controlled (default) normalization.
        allfeature_sizefactors = adata.uns.get('deseq2_sizefactors_allfeatures')
        if allfeature_sizefactors is None or 'nreads_total_raw' not in adata.obs.columns:
            msg = ('No all-feature-controlled DESeq2 size factors found in AnnData object '
                   '(uns[\'deseq2_sizefactors_allfeatures\'] missing, or obs[\'nreads_total_raw\'] '
                   'unavailable). Skipping combined tRNA + non-tRNA PCA plot. Rebuild with the '
                   'current `trnagraph analyze build` to enable this plot.')
            if threaded:
                threaded += msg + '\n'
            else:
                logger.warning(msg)
        else:
            sf_map = {k: float(v) for k, v in dict(allfeature_sizefactors).items()}
            total_df = pd.DataFrame(adata.obs, columns=['trna', 'sample', 'nreads_total_raw'])
            missing_samples = sorted(set(total_df['sample'].astype(str)) - set(sf_map.keys()))
            if missing_samples:
                msg = f'All-feature size factors not found for samples {missing_samples}; defaulting to 1.0 for the combined PCA plot.'
                if threaded:
                    threaded += msg + '\n'
                else:
                    logger.warning(msg)
            total_df['nreads_total_norm_allfeatures'] = total_df['nreads_total_raw'] / total_df['sample'].astype(str).map(lambda s: sf_map.get(s, 1.0))
            total_df = total_df.pivot_table(index='trna', columns='sample', values='nreads_total_norm_allfeatures', observed=True)

            shared_cols = total_df.columns.intersection(nontrna_df.columns)
            if len(shared_cols) < len(total_df.columns) or len(shared_cols) < len(nontrna_df.columns):
                missing_from_nontrna = set(total_df.columns) - set(nontrna_df.columns)
                missing_from_total = set(nontrna_df.columns) - set(total_df.columns)
                msg = f'Sample columns did not fully match between tRNA and non-tRNA counts for combined PCA; using intersection of {len(shared_cols)} shared samples.'
                if missing_from_nontrna:
                    msg += f' Missing from nontRNA_counts: {sorted(missing_from_nontrna)}.'
                if missing_from_total:
                    msg += f' Missing from tRNA total counts: {sorted(missing_from_total)}.'
                if threaded:
                    threaded += msg + '\n'
                else:
                    logger.warning(msg)

            if len(shared_cols) == 0:
                msg = 'No overlapping sample columns between tRNA and non-tRNA counts. Skipping combined PCA plot.'
                if threaded:
                    threaded += msg + '\n'
                else:
                    logger.warning(msg)
            else:
                combined_df = pd.concat([total_df[shared_cols], nontrna_df[shared_cols]], axis=0)
                threaded = _generate_pca_plots(combined_df, hue_dict, colormap, pcamarkers, pcacolors, output, basename=f'allRNA_{pcamarkers}_by_{pcacolors}', title_suffix='All tRNA + Non-tRNA RNAs', threaded=threaded)

    if threaded:
        return threaded

if __name__ == '__main__':
    pass
