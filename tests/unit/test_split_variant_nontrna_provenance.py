"""Regression test for roadmap.md's small-RNA/non-tRNA split-variant bug. merge_variant_into_adata()
previously stored each split variant's OWN independently-recomputed non-tRNA counts (from that
split's own length-filtered counting output) into uns['size_splits'][tag]['nontRNA_counts'], which
build_variant_view() then overlays as uns['nontRNA_counts'] whenever a split (--variant norm:o60 /
norm:u60) is rendered. A length cutoff cannot meaningfully partition non-tRNA reads the way it
partitions tRNAs (they're not being classified by that criterion at all), so one split's non-tRNA
counts could come out empty while the other's didn't -- exactly the observed "o60 got the small-RNA
subplots, u60 didn't" bug. Fix: a split variant's non-tRNA counts must always be the same value as
the object's existing full/unsplit variant, never a fresh per-split recomputation."""
import pandas as pd

from trnagraph.modules.adataBuild import _resolve_full_variant_nontrna_counts


def test_resolve_full_variant_nontrna_counts_prefers_the_existing_full_variant_value():
    full_value = pd.DataFrame({'count': [10, 20]}, index=['geneA', 'geneB'])
    split_specific_value = pd.DataFrame({'count': [0, 0]}, index=['geneA', 'geneB'])

    resolved = _resolve_full_variant_nontrna_counts(full_value, split_specific_value)

    pd.testing.assert_frame_equal(resolved, full_value)


def test_resolve_full_variant_nontrna_counts_falls_back_to_the_split_value_when_no_full_value_exists():
    split_specific_value = pd.DataFrame({'count': [5]}, index=['geneA'])

    resolved = _resolve_full_variant_nontrna_counts(None, split_specific_value)

    pd.testing.assert_frame_equal(resolved, split_specific_value)
