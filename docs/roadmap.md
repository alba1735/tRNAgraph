# Project Roadmap

This document outlines the planned features and current areas of focus for tRNAgraph development, organized into four remaining phases plus a standing testing effort. Phase 1 (correctness/data-integrity) is complete; Phase 2 (known issues/refactoring), Phase 3 (the Sprinzl position framework), Phase 4 (the multivariate/timecourse feature), and Phase 5 (legacy script migration/misincorporation analysis) can proceed in parallel, whenever someone is ready to pick any of them up — Phase 3 was split out of Phase 2 once the size of the Sprinzl-numbering overhaul became clear, and Phase 5 was split out separately since it's expected to grow with additional not-yet-incorporated analysis code. Resolved items are removed from this document once done rather than kept as a changelog — see git history for that record, or `changes_from_trax.md` for the subset of changes that affect statistical output specifically.

Phase 1 is complete: the VST hang, the `obs`/`obsm`/`obsp`/`uns` mapping audit, replacing manual log2FC with native DESeq2 output, and all 5 `MigrationTesting` follow-up items (all resolved as either the deterministic-sorting accepted divergence or architecture differences with no tRAX equivalent — see `MigrationTesting/docs/validation_status.md` and `changes_from_trax.md`). See "Implementation Notes" below before writing new PyDESeq2/AnnData code — each fix ran into a non-obvious gotcha worth not re-discovering.

## Implementation Notes (PyDESeq2 / AnnData gotchas)

Landmines hit while doing Phase 1's work, relevant to any future code that touches PyDESeq2, AnnData slots, or log2FC/DE -- most directly Phase 4's multi-factor DE engine, which will construct fresh `DeseqDataSet`s the same way these fixes did.

- **PyDESeq2 + zero-heavy tRNA count data hangs on its default size-factor method.** `DeseqDataSet`'s default `size_factors_fit_type='ratio'` silently falls back to its `'iterative'` method whenever every feature has at least one zero-count sample -- normal for tRNA coverage data -- and that method's cost (a `scipy.optimize` Powell search over one parameter per _sample_) blows up non-linearly with sample count: fine at ~50 samples, unresponsive past ~100. Any new code constructing a `DeseqDataSet` (Phase 4's multi-factor engine included) should pass `size_factors_fit_type='poscounts'` explicitly rather than relying on the default. See `adataBuild.py._compute_vst_` and `toolsTG.py.adataLog2FC.log2fc_df` for the fix in place.
- **Reusing externally-precomputed size factors needs more than `obsm['size_factors']`.** `DeseqDataSet.vst_fit()` only skips its own size-factor recompute when `self.logmeans` is also set (not just `obsm['size_factors']`) -- `logmeans` is otherwise only ever set inside `fit_size_factors()` itself, so setting `obsm['size_factors']` alone gets silently overwritten. See the workaround in `adataBuild.py._compute_vst_` if a future change needs to inject size factors again.
- **AnnData slot convention: per-obs-aligned data goes in `obsm`, not `uns`.** Anything with one row per (a subset of) `adata.obs_names`, reindexable onto the full obs axis, belongs in `obsm` -- `uns` is for data that genuinely doesn't conform to any axis (different cardinality/shape entirely). See `adataCluster.py.adataCombine`'s sample-vs-group split for a worked example of telling the two cases apart.
- **DE code must feed PyDESeq2 raw counts, never pre-normalized ones.** PyDESeq2 does its own internal normalization; feeding it an already-normalized (`_norm`) column double-normalizes and biases the fit. Every `adata.obs` readtype column has a raw/normalized pair (`_raw`/`_norm`); always derive the `_raw` column for anything DESeq2-based.
- **A "no results" edge case must still return the expected columns.** A DE/stats function that can legitimately produce zero result rows (e.g. nothing passes a cutoff) must still return the same column set callers index unconditionally downstream -- an empty result with the wrong (missing) columns caused a real `KeyError` in `plotsVolcano.py` during this work. Keep this in mind for Phase 4's Venn/PERMANOVA "no significant hits" cases.

## Phase 2 — Known Issues & Refactoring

### Refactoring Targets

- **Trimmer UMI Handling**: Confirmed real bug, not just an untested assumption — for `--umi3`, `toolsTrim.py`'s code does a bare `pass`, so no `--umi_loc` is passed to fastp at all, a silent no-op. Confirmed via fastp's own docs this can't be fixed with a fastp flag: every `--umi_loc` value (`index1/index2/read1/read2/per_index/per_read`) is 5'/head-anchored only, with no tail-anchored option. Fix: keep fastp's native `--umi_loc read1` for 5' UMIs (the default case, and the one that lets fastp fully replace tRAX's `umi_tools` usage); for `--umi3`, reintroduce a `umi_tools extract --3prime` post-processing step after fastp, mirroring tRAX's original approach (`umi_tools` is already a pinned conda dependency in `requirements.yaml`, so no new dependency is needed).
- **Logging**: Standardize logging across scripts for better traceability (currently `print()` is used in 37 of 42 module files; `toolsTestSuite.py` and `toolsTrim.py` use `logging`). Convert incrementally, only in files already being touched for other fixes in this phase — not as a dedicated 37-file sweep.
- **Warnings**: A full repo-wide catalog found exactly one warning suppression in the entire codebase: `adataGraph.py`'s blanket, undocumented, import-time suppression of a matplotlib "identical low and high ylims" `UserWarning`. Scope it to the specific plotting call that triggers it (a flat-line/zero-variance case) with a comment explaining why it's safe, rather than suppressing globally. Separately, the `FutureWarning`s during AnnData creation (likely from `.fillna(0.0)` without an explicit dtype at `adataBuild.py:1280`, and/or a non-`observed` `groupby`) should get a root-cause fix — explicit dtype/`observed=` parameters — not a new suppression.
- **`--lazy` flag**: Remove the flag. Skip-mapping-if-BAM-exists-with-header-check becomes the unconditional default behavior; add a new `--force-remap` flag as the explicit opt-out for a full remap.
- **`--nofrag` flag**: Remove entirely. It's real ported tRAX behavior (originally for TGIRT-only full-length libraries) but the `nofrag=True` collapsed-count path has zero test coverage or real invocation since migration, and doesn't fit workflows like OTTR-seq that produce mixed fragment/full-length reads. Its other tRAX-side half (QC-threshold recalibration for TGIRT mode) was ported into `toolsQC.py` but never wired to anything — delete that dead code too.
- **`--nosizefactors` flag**: Remove entirely. Confirmed currently broken: setting it crashes `analyze build` with an uncaught `FileNotFoundError` (`AnnDataBuilder.__init__` unconditionally reads size-factor/normalized-count files that this flag skips writing), and no fallback normalization exists anywhere in that path — this isn't just redundant, it's non-functional. DESeq2 size factors are always computed at build time going forward (the Phase 1 `poscounts` fix already resolved the performance concern that likely motivated skipping it); add a raw-vs-normalized display toggle at `graph` time instead.
- **`map` makes extra directories**: Both `preprocess map` and `analyze build` unconditionally create a `graphsdir` that's never written to by either command (`build` already has a comment marking this as dead legacy scaffolding) — remove the dead creation from both now.
- **`-maponly` flag**: Remove from `tools test`. Confirmed a pure convenience shortcut for stopping the demo pipeline after mapping, with no unique test-harness value (no extra assertions/validation) beyond what the already-separate `--trim`/`--makedb`/`--map` flags on `tools test` provide.
- **tqdm**: implement this for progress updates on long-running per-sample loops (mapping, trimming, counting).
- **Directory structure**: Narrow scope for now to the `graphsdir` removal above. A broader audit of directory creation across all modules (this item's original full ask) is real but open-ended — track as a separate future item once there's a concrete list of what each command actually needs.
- **Split-file variance**: Root cause identified. The Phase 1 deterministic-sort fix (`counting_logic_changes.md`) was only applied to the `trnaloci` candidate-set call site in `toolsCountReads.py`; two more call sites (`featurelist`, `otherseqlist`) still use unsorted, hash-order-dependent `RangeBin.getbin()` iteration with first-match-`break` semantics. Fix: apply the same deterministic-sort pattern to both. Validate via a same-input-run-twice reproducibility test in `tests/unit/` (not a `trna-test` audit against tRAX — tRAX has the identical non-determinism, so there's no single "correct" tRAX output to diff against for this one).

## Phase 3 — Sprinzl Position Framework

Split out of Phase 2 once its scope became clear: this isn't a one-line labeling fix but an overhaul of tRNAgraph's Sprinzl position-to-structural-region annotation system (`adata.var['location']` and the per-domain position lists it's derived from), which several other parts of the pipeline read from (CSV exports, coverage plots, the strand-orientation heuristic in `adataBuild.py`).

- **Cross-tool consistency audit (do first)**: Before finalizing any fix, review how tRNAscan-SE (whose output already feeds `preprocess makedb` directly), tDRnamer, and other Lowe Lab tRNA tools compute and label Sprinzl positions across domains/species. Check for absolute consistency with what tRNAgraph implements or plans to implement, and identify any existing position-mapping code, tables, or conventions that could be reused/shared rather than reimplemented from scratch — tRNAgraph's database is already built directly from tRNAscan-SE's own output, so drift between the two tools' position conventions is a real risk worth ruling out before investing in a from-scratch overhaul.
- **Acceptor stem labeling fix**: `adata.var['location']`'s `loc_dict` (`adataBuild.py`, ~lines 1029-1075) labels the acceptor stem asymmetrically — 5' side `range(-1,8)` vs. 3' side `range(66,77)` — incorrectly folding the discriminator base (73) and CCA tail (74-76) into the stem, plus a phantom `'0'` key that never matches (Sprinzl numbering jumps -1 to 1). Fix: correct the range to the canonical 1-7/66-72 pairing, give the discriminator base and CCA tail their own labels, drop the phantom key. Confirmed isolated to the acceptor stem — D-stem, anticodon-stem, and T-stem are all correctly symmetric already, vetted against canonical Sprinzl numbering (tRAX has no equivalent named-region labeling to diff against, so this is validated against the published numbering scheme instead, and against the cross-tool audit above). Also fix `docs/data_structure.md`'s documented location category names, which don't match any of the actual code strings (stale docs, independent of the range bug). Document the change and its rationale in `improvements_from_trax.md` once implemented.
- **e-series / variable-arm extension**: e-series variable-arm insertion codes are capped at e19 in both `loc_dict` and the underlying per-domain position lists (`toolsGetCoverage.py`), while canonical class-II tRNAs (e.g. tRNA-Leu, tRNA-Ser) need up to e27 — those variable arms are silently truncated from coverage entirely, a data-completeness gap rather than a labeling bug. Must be extended per-domain (Euk/Bact/Archaea already have separate position lists in `toolsGetCoverage.py`; this can't be a single universal range) and needs real eukaryotic/archaeal test data to validate against actual lab usage (mouse/human/archaea), not just `vibrChol1`. `MigrationTesting`'s `trna-test setup-data` currently only supports `--organism bact` — extending that is a prerequisite for validating this properly.

## Phase 4 — Multivariate / Timecourse Analysis

tRAX only handles one dataset at a time, run repeatedly to build up different views. tRNAgraph is designed for much more robust exploration of multidimensional experiments — the core motivating example is a cell-to-organoid timecourse with samples at multiple timepoints (or more generally, any multi-factor design: timepoint, cell type, drug treatment, etc.), and comparing fragment (`u<N>`, e.g. under 60nt) vs. full-length (`o<N>`, e.g. over 60nt) tRNA expression within that design.

This phase covers building the analysis machinery. **A separate, dedicated design pass is required before implementation begins** — this roadmap entry scopes the pieces and their dependencies, not the final algorithms. Two concrete deliverables are prerequisites for implementation and are tracked here but not written as part of this roadmap refresh:

- A use-cases documentation file, since not every dataset needs multivariate or timecourse analysis (single-condition experiments should be unaffected).
- Config support for declaring category order (e.g. correct chronological timepoint order), since order isn't always inferable from the data itself.

### 4a. Fragment vs. full-length comparison (standalone)

A baseline comparison that works for any single grouping, independent of study design — tRNAs are expressed as both full-length and fragments (e.g. in OTTR-seq data), so comparing frag vs. full expression across any group/timepoint/celltype is useful on its own.

- Venn diagrams of DE-hit overlap between `u<N>`/`o<N>` splits for a given comparison, using the `venn` PyPI package. A working prototype of the dict-of-sets → `venn()` pattern already exists in `test_organoids/organoid_analysis.ipynb` and should be ported/generalized rather than built from scratch.

### 4b. General multi-factor DE engine

The current DE machinery is single-factor only: `design_factors="condition"` is hardcoded in the PyDESeq2 build-time path (`adataBuild.py`), and every grouping consumer (PCA colors, volcano groups, `tools log2fc -g`, clustering's hardcoded `grpby='group'`) takes exactly one `adata.obs` column. Arbitrary extra metadata columns (e.g. `timepoint`, `celltype`) can already be added to `adata.obs`, but nothing lets them combine into a real multi-factor design.

- Build real multi-factor DESeq2 design formulas from N arbitrary `obs` columns (e.g. `~timepoint+celltype`), correctly modeling main effects and interactions rather than collapsing columns into a synthetic composite category.
- Extend the single-column grouping consumers (plotting, clustering, `log2fc`) to work with this multi-factor engine.
- Absorbs Phase 2's former "External Analysis Integration" item (running arbitrary DESeq2 analyses directly on an existing `.h5ad` as a `tools`-pattern command): the closest existing building block, `adataLog2FC`, has `design_factors` hardcoded to `'condition'` in both cases, so a narrowly-scoped Phase 2 version would either be trivial or duplicate this engine — tracked here instead of as a separate item.

### 4c. Trajectory / decoupling analysis

For ordered conditions (timecourse, dose-response), determine whether the frag:full expression ratio changes as the condition changes — e.g. does fragmentation stay proportionally stable across a timecourse, or does it decouple from full-length expression? Whether this applies to a timecourse or to a simpler two-condition comparison (e.g. drug vs. no drug), the underlying question is the same: does the split ratio respond to the condition.

- Model this as an interaction term (split × ordered-condition) on top of the 4b multi-factor engine, testable via DESeq2's likelihood ratio test. Exact algorithm/statistical framing is deferred to the separate design pass mentioned above.

### Multivariate statistics: PERMANOVA

- No multivariate-stats dependency currently exists in the project (no `statsmodels`, no `scikit-bio`). Add **`scikit-bio`** for PERMANOVA (a cited, maintained implementation was preferred over a custom implementation of the stats algorithm; dependency weight is not a concern at this project's scale).
- Default distance metric: **Bray-Curtis**, since the tRNA count data is zero-heavy and low-count for many features — Euclidean distance is distorted by zero-inflation in this kind of sparse compositional data.

### Sample/group ordering

No mechanism currently lets a user specify explicit category order anywhere in tRNAgraph — every plotting module falls back to `sorted()` (alphabetical), and the two existing config files (`--config` filter config, `--colormap`) are transient, per-`graph`-call JSON with no ordering concept.

- Bake explicit order into `adata.obs` as an ordered `pd.Categorical` at `analyze build` time, so there's a single source of truth reused both for plot legend/axis order and for DESeq2 reference-level ordering (rather than a transient flag that would need to be re-supplied on every plot/DE call, and wouldn't cover the statistics side at all). Ties into the Phase 1 `obs`/`obsm`/`obsp`/`uns` mapping audit.

## Phase 5 — Legacy Script Migration & Misincorporation Analysis

- **Legacy script cleanup** (low priority): none of the 12 `plotsLegacy*.py` files are wired into `cli.py` — they're fully dead/orphaned code today, reachable from nothing. Catalogued disposition:
  - Safe to delete now (fully superseded by a native `plots*.py` module): `plotsLegacyCoverage.py`, `plotsLegacyPCA.py`, `plotsLegacyVolcano.py`, `plotsLegacyLocusCoverage.py`, `plotsLegacyFeatureTypes.py`, `plotsLegacyFeatureTypesReal.py`, `plotsLegacyGeneFeatures.py`.
  - Needs verification before deciding (partial overlap with a native module): `plotsLegacyMismatchBoxplot.py`, `plotsLegacyScatter.py`.
  - Not yet ported (no native equivalent exists) — port to tRNAgraph's own style/conventions, then delete the legacy source: `plotsLegacyReadLength.py` (read-length distribution), `plotsLegacyMismatch.py` (mismatch-count-by-type), `plotsLegacyCCA.py` (3' end/CCA composition).
- **Misincorporation analysis** (not yet in this repo): additional misincorporation-analysis code exists already, written outside this repo, and is planned to be incorporated here. Given the overlap with the still-unported `plotsLegacyMismatch.py`/`plotsLegacyMismatchBoxplot.py` items above, scope both together when this work is picked up rather than porting the legacy mismatch plots first and redoing them once the new analysis lands.

## Testing

`tests/unit/` now holds real pytest coverage (`pytest` added as a `dev` optional-dependency in `pyproject.toml`), started alongside the three completed Phase 1 fixes above: `test_vst.py` (the VST hang + size-factor-preservation regression), `test_adata_cluster_slots.py` (the obsm/uns placement fix), and `test_log2fc.py` (the DESeq2-based log2FC replacement, including the empty-result column-set regression). `toolsTestSuite.py` remains a separate integration/demo pipeline runner (the `trnagraph tools test` end-to-end demo, now defaulting its output to `test_vibrChol1/` instead of `tests/`, so generated run output never collides with `tests/unit/`'s real test code), not a substitute for unit coverage.

- Keep building out unit test coverage, prioritizing correctness-critical and statistically risky code first.
- Phase 4's new modules (multi-factor DE, Venn comparison, PERMANOVA, trajectory analysis) should get test coverage as each lands, rather than retrofitting tests after the whole feature is "done" — this is new statistical code and the highest-risk work on this roadmap.
- The same as-it-lands bar applies to Phase 2's correctness-bearing items (merge conflicting-run-info enforcement, split-file variance reproducibility, metadata sample-name validation) — cosmetic/DX-only items (logging, tqdm, trim-plot styling) don't need new unit tests.
- Phase 3's Sprinzl fixes (acceptor-stem labeling, e-series extension) get the same bar, validated against canonical Sprinzl numbering and the tRNAscan-SE/tDRnamer cross-tool audit rather than a `trna-test` diff against tRAX, since tRAX has no equivalent labeling to diff against.

## Unvetted / future items (Ignore for now)

- Improve entire test suite's CI/CD coverage (GitHub Actions) — currently only `tests/unit/` is run, not `toolsTestSuite.py` or any of the legacy `trna-test` scripts. This will need to be implemented after everything else and after careful consideration of how to handle the large, generated `test_vibrChol1/` output directory (which is currently ignored in `.gitignore` and not checked into the repo).
  - Also incorporate Archea and Eukaryotic test data into the CI/CD pipeline, since the current `test_vibrChol1/` demo is bacterial-only.
- Add a Google Collab notebook (or other tutorial) for new users to run through the pipeline without installing anything locally.
  - We can use the one I created for tRAX as a guidline, but it will need to be updated for tRNAgraph's new CLI and new features etc.
- Go through the entire tRAXs issue backlog on github and see if any of those issues are still relevant to tRNAgraph, and if so, add them to this roadmap.
- Better terminal outputs in general that are more readable and informative, especially for long-running processes (e.g. mapping, trimming, counting).
- Better logs (markdown) for long-running processes (e.g. mapping, trimming, counting) that are more readable and informative.
- Automatic gtRNAdb db creation or a curated list of pre-built dbs for common organisms (e.g. human, mouse, E. coli, S. cerevisiae, etc.) to avoid the need for users to build their own dbs.
- Modification analysis of some kinda... This is dependent on sequencing type and species, so it will need to be a separate module that is only run when the user has the correct data type and organism. This is a long-term goal and will require more research and development.
- update command to cli so that users can update the software.
- Better mitochondrial tRNA handling, since mito tRNAs are often very different from nuclear tRNAs and may require special handling in the pipeline.
  - Probably an extra flag at build/map to perform additional mapping/analysis for mito tRNAs, and/or a separate mito tRNA database. This will require some consideration of how to handle mito tRNAs in the context of the rest of the pipeline, since they are often very different from nuclear tRNAs and may require special handling.
- Mouse/Rat have tons of SINEs that are tRNA-derived which caused huge problems in tRAX, and will likely cause problems in tRNAgraph as well. The work around was using the High Confidence tRNA set from GtRNAdb, but this is not a perfect solution.
  - Some SINE analysis might be interesting.
