# Changes from tRAX

This document tracks statistical methodology differences, logic changes, and implementation improvements between tRNAgraph and tRAX (the legacy pipeline it replaces). It is a technical record, not a feature list — entries describe what changed and why, not how much better the result is. For validation status (which outputs match tRAX numerically and which differ), see `MigrationTesting/docs/validation_status.md`. For architectural background, see `MigrationTesting/docs/pipeline_differences.md`.

Entries are added as changes are made; resolved/superseded entries are not removed, since this is a historical record (unlike `roadmap.md`, which stays forward-looking only).

## Statistical Methodology Differences

### Read classification: deterministic sorting vs. non-deterministic set iteration

tRAX classifies a read that overlaps multiple genomic features (e.g. a pre-tRNA region overlapping a mature tRNA locus) via Python `set` iteration, which is hash-order dependent and effectively random from run to run. tRNAgraph sorts candidate features deterministically instead, so classification is reproducible.

Effect: ~0.1-0.5% of reads may be classified into a different category than a given tRAX run would have assigned. This concentrates in categories with the most feature overlap — pre-tRNA vs. mature tRNA is the primary case; anticodon/amino-acid count differences from this cause are on the order of single-count discrepancies. See `MigrationTesting/docs/counting_logic_changes.md` for the algorithm-level diff.

### Differential expression: PyDESeq2 vs. R DESeq2

tRAX's DESeq2 calls run through R via `Rscript` subprocesses. tRNAgraph uses PyDESeq2 (a Python reimplementation of the same negative-binomial GLM methodology) directly, with no R dependency anywhere in the pipeline. Results are statistically comparable but not bit-identical, due to implementation differences in dispersion estimation and optimization between the two packages.

### Primary size factors: tRNA/tRX-controlled vs. all-feature-controlled

tRAX's `SizeFactors.txt` is always computed by DESeq2's `estimateSizeFactors()` over all features (tRNAs plus any other GTF-annotated features). tRNAgraph's primary size factors (`adata.obs['deseq2_sizefactor']`, `adata.uns['deseq2_sizefactors_trna']`, on-disk `<exp>-SizeFactors.txt`) are computed using only tRNA/tRX features as the normalization reference instead, so non-tRNA feature abundance changes don't distort tRNA normalization. The all-feature-driven computation (tRAX's method) is still run and kept as `adata.uns['deseq2_sizefactors_allfeatures']` / `<exp>-allfeature_SizeFactors.txt` for direct tRAX-parity comparison.

### Trimming: fastp vs. cutadapt/SeqPrep

tRAX trims adapters with SeqPrep/cutadapt. tRNAgraph uses fastp. The two tools' adapter-detection and quality-trimming logic differ slightly at read edges, which propagates into small differences in downstream mapped-read counts. `MigrationTesting`'s test suite runs tRNAgraph against both trimmers (`vibrChol1_fastp`/`vibrChol1_cutadapt`) so this effect can be isolated from other differences during validation.

### Variance-stabilizing transform: VST vs. rlog

tRAX's R pipeline produces an rlog-transformed count matrix (`<exp>-rlog.txt` / `<exp>-<type>_rlog.txt`, via DESeq2's `rlog()`). tRNAgraph implements PyDESeq2's VST (`adata.layers['vst']`) instead of rlog; the two are alternative variance-stabilizing transforms serving the same purpose (DESeq2's own documentation recommends VST over rlog for anything but small sample counts, given rlog's higher computational cost). tRNAgraph does not produce an rlog-equivalent file.

### Low-coverage feature filtering: `--mincoverage` (hard row drop) vs. `--minfeaturereads` (VST-fit-only exclusion)

tRAX's `--mincoverage` (default 30) and tRNAgraph's original, identically-behaving equivalent both worked as a hard filter: a tRNA gene whose total raw read count (summed across all samples) fell below the threshold was dropped from the coverage output file entirely. Because tRNAgraph's `adata.obs` row universe is built directly from that same coverage file's index, this silently removed the gene from the *entire* AnnData object -- not just the coverage/Sprinzl-position matrix, but the gene-level columns volcano/heatmap plots read from too.

As of 2026-08-20, renamed to `--minfeaturereads` and re-scoped to affect only `adata.layers['vst']`. Every gene now always gets a full raw coverage row, a full normalized coverage row, and its own VST value regardless of this threshold. The threshold instead controls only which genes' counts are allowed to influence the VST dispersion-trend *fit* (`AnnDataBuilder._compute_vst_`'s `fit_mask`): a handful of noisy, very-low-count genes can otherwise distort the fitted trend curve for every other gene too, so they're excluded from the fit itself, then transformed using the resulting fit like everything else -- standard DESeq2 trend-shrinkage methodology (fit on a well-behaved subset, apply to all), not a compromise. Which genes were excluded from the fit is recorded on `adata.obs['vst_fit_excluded']`. `--volcutoff`/`--heatcutoff` remain the correct, independent knobs for low-count protection on the DE/plotting side; this flag is unrelated to those.

### Differential expression: manual t-test vs. PyDESeq2 GLM (log2FC path)

Through 2026-08-19, tRNAgraph's `adata.obs`-column-driven log2FC computation (used by heatmap/volcano plots and the `tools log2fc` command's default `group` comparison) computed log2 fold-change and significance via a two-sample t-test (`scipy.stats.ttest_ind_from_stats`) on precomputed per-group mean/std/count, with no multiple-testing correction. This was a tRNAgraph-internal implementation choice, not something inherited from tRAX (tRAX does not have an equivalent ad hoc per-`obs`-column log2FC path).

As of 2026-08-19, this path fits a PyDESeq2 `DeseqDataSet`/`DeseqStats` per comparison instead — the same statistical methodology tRNAgraph already uses for its build-time DESeq2 output, applied consistently to this on-demand path. The reported p-value changed from an uncorrected raw p-value to `padj` (Benjamini-Hochberg adjusted). See `toolsTG.py`'s `adataLog2FC.log2fc_df`.

## Outputs Not Produced by tRNAgraph

tRAX outputs with no tRNAgraph equivalent, because the underlying mechanism doesn't exist in tRNAgraph's architecture (see `MigrationTesting/docs/validation_status.md` for the full list and reasoning):

- `trimindex.txt` — tRAX's SeqPrep/cutadapt trim-run manifest bookkeeping; tRNAgraph's fastp-based trimmer uses a different manifest convention.
- `Rlog-<expname>.txt` — captured stdout/stderr of tRAX's `Rscript` subprocess calls; tRNAgraph has no R subprocess to capture.
- `positionmismatches.txt` — generated by an R script (`boxplotmismatches.R`); the underlying per-position mismatch data is stored in tRNAgraph's AnnData object instead (see below).

## Data Storage Changes

### Flat files vs. AnnData

tRAX writes each pipeline stage's output to a flat text file. tRNAgraph assembles the same information into a single AnnData object (`obs` = tRNA gene × sample observations, `var` = Sprinzl position × coverage-type features, `X`/`layers` = count matrices, `uns` = unstructured/aggregate results). Flat-file outputs matching tRAX's format are still written for `analyze build`'s primary DESeq2 runs (for validation and for users who want tRAX-compatible files directly), but downstream analysis in tRNAgraph (plotting, clustering, further DE) operates on the AnnData object, not the flat files.

Per-position mismatch/base-composition data (mismatch count, deletions, per-nucleotide coverage) that tRAX only exposes via `mismatches.txt` and R-generated `positionmismatches.txt` is stored at full per-position granularity in the AnnData object's coverage-type layers (`mismatchedbases`, `deletedbases`, `adenines`, `thymines`, `cytosines`, `guanines`, `deletions`); tRNAgraph's flat-file `mismatches.txt` output is kept for tRAX parity but is not the only place this data lives.

## CLI Flag Clarity Fixes

Renames/documentation fixes for a flag whose name or help text was ambiguous, carried over unchanged from tRAX, with no behavior change (default stays the same as tRAX's).

### `--uniqueonly` → `--filtermultimapped` (`analyze build`)

tRAX's `processsamples.py` has an identically-ambiguous `--uniqueonlycov` flag (default `False`, help text just "Show only unique coverage"), and tRNAgraph inherited both the name and the ambiguity unchanged. The problem: tRNAgraph already has a completely different, unrelated "uniqueness" concept that's always computed regardless of this flag — whether a read maps unambiguously to one specific tRNA transcript (as opposed to sharing an anticodon/amino acid with other tRNAs), tracked via `isuniqueaminomapping()`/`isuniqueacmapping()`/`isuniquetrnamapping()` and surfaced in the always-generated `results/unique/`/`graphs/unique/` output. The flag being renamed controls a second, orthogonal axis instead: whether a read maps to exactly one location in the *genome* at all (`issinglemapped()`) — when set, genomically-multi-mapped reads are dropped from every coverage column before any computation happens (main `coverage`, `readstarts`/`readends`, mismatch tracking, all of it), not filtered into a separate file. Renamed to `--filtermultimapped` to describe the actual mechanism (drop genomically-multi-mapped reads) without colliding with the unrelated tRNA-identity "unique" terminology already in use elsewhere in the output. Default unchanged (`False`, matching tRAX's own default/recommended behavior of including multi-mapped reads).

## Implementation Notes

Fixes made during tRNAgraph development, not differences from tRAX, but recorded here since they affect statistical output and are the kind of thing that could silently regress:

- **VST hang on larger sample counts (fixed 2026-08-19)**: PyDESeq2's `vst_fit()` only skips its own internal size-factor recomputation when both `obsm['size_factors']` is set and `self.logmeans` is not `None`; tRNAgraph was setting only the former, so its pre-computed tRNA-control size factors were silently discarded and PyDESeq2 recomputed its own via an `"iterative"` method whose cost scales non-linearly with sample count on zero-heavy data (fine at ~50 samples, unresponsive past ~100). Fixed in `adataBuild.py`'s `_compute_vst_`.
- **Cluster/UMAP coordinates stored in the wrong AnnData slot (fixed 2026-08-19)**: sample-level UMAP/HDBSCAN cluster coordinates were stored wholesale in `adata.uns` instead of `adata.obsm`, despite being per-observation-aligned data. Moved to `obsm` (namespaced per split variant where applicable). Fixed in `adataCluster.py`'s `adataCombine`.
