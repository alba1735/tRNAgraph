"""Regression tests for toolsTG.resolve_nontrna_counts, the shared gate used by
plotsVolcano.py/plotsPca.py/plotsCorrelation.py's non-tRNA/combined plots (roadmap.md's
small-RNA/non-tRNA split-variant bug). Non-tRNA reads can't be meaningfully partitioned by the
same length cutoff used to split tRNAs, so non-tRNA/combined plots must only ever be generated
for the complete (non-split) variant -- regardless of whether uns['nontRNA_counts'] happens to be
present for a split variant's view."""
import pandas as pd
import anndata as ad
import numpy as np

from trnagraph.modules.toolsTG import resolve_nontrna_counts


def _make_adata(nontrna_df=None):
    adata = ad.AnnData(X=np.zeros((1, 1)))
    if nontrna_df is not None:
        adata.uns['nontRNA_counts'] = nontrna_df
    return adata


def test_returns_the_counts_when_full_variant_and_data_present():
    df = pd.DataFrame({'count': [1, 2]}, index=['geneA', 'geneB'])
    adata = _make_adata(df)

    result_df, skip_message = resolve_nontrna_counts(adata, is_full_variant=True, feature_label='PCA plots')

    pd.testing.assert_frame_equal(result_df, df)
    assert skip_message is None


def test_skips_with_missing_data_message_when_full_variant_but_no_data():
    adata = _make_adata(None)

    result_df, skip_message = resolve_nontrna_counts(adata, is_full_variant=True, feature_label='PCA plots')

    assert result_df is None
    assert '--gtf' in skip_message


def test_skips_split_variants_even_when_data_is_present():
    """The core fix: a split variant must never get non-tRNA plots, even if nontRNA_counts
    happens to be non-empty on this view."""
    df = pd.DataFrame({'count': [1, 2]}, index=['geneA', 'geneB'])
    adata = _make_adata(df)

    result_df, skip_message = resolve_nontrna_counts(adata, is_full_variant=False, feature_label='PCA plots')

    assert result_df is None
    assert skip_message is not None
    assert 'split' in skip_message.lower()


def test_feature_label_is_used_in_the_skip_message():
    adata = _make_adata(None)

    _, skip_message = resolve_nontrna_counts(adata, is_full_variant=True, feature_label='correlation matrices')

    assert 'correlation matrices' in skip_message
