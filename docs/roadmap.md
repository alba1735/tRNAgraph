# Project Roadmap

This document outlines the planned features and current areas of focus for tRNAgraph development.

## Planned Features

### Analysis & Statistics

- ~~**Alternative Normalization**: Add options to compute plots with RAW counts or alternative layers (normalization).~~ Resolved: `graph`/`analyze cluster`/`tools log2fc` now accept a `--variant <norm>:<tag>` flag (`<norm>` = `norm`/`raw`/`allfeatures`/`vst`) that selects the underlying layer a plot/DE call is built from, e.g. `--variant raw:full`. See `data_structure.md#split-variants---readlengthsplit`.
- **External Analysis Integration**: Allow new DESeq2 or other analyses to be run directly on the combined object as a tool.
- ~~**Differential Expression**: Add functionality to redo differential expression analysis, and/or add new layers and methods as a tool.~~ Largely resolved: `tools log2fc` already exists as a standalone DE-recompute command and is now `--variant`-aware, so it can redo DE for any normalization/split variant. A more general "run new DE methods" tool is still open.
- **Single-Object Size-Range Splits**: ~~Currently calculating for size range splits (over/under Nbp) creates separate adata objects.~~ Resolved: `analyze build --readlengthsplit N` now adds `u<N>`/`o<N>` variants to the _same_ object as layers/`obsm`/`uns` entries instead of writing separate `_uN.h5ad`/`_oN.h5ad` files, and a new `analyze addsplit` command adds further cutoffs to an existing object incrementally, without disturbing previously-added variants. On-disk `results_uN`/`graphs_uN` directories are unchanged. See `data_structure.md#split-variants---readlengthsplit`.

### Visualization Enhancements

- **Sorting**: Add functionality to sort combined coverage plots from most to least abundant.
- **Radar Plots**: Improve amino acid dictionary generation to be species-universal (currently uses a subset).
- **Styling**: Cleanup `trim` plots to match the general tRNAgraph aesthetic.
- **Correlation**: Add function to plot correlation with smallRNAs included.

### Data Validation

- **Metadata Check**: Add validation to ensure `sample` names in metadata match those in the samples file before writing the AnnData object.
- **Parameter Fallback**: Add checks for graphing parameters; default to `sample` if the requested column does not exist in `adata`.

## Known Issues & Refactoring Targets

- **Acceptor Stems**: Fix labeling issues with acceptor stems in `adata.var` (affects CSV exports).
- **Merge Logic**: Prevent merging of AnnData objects if they have conflicting tRNAgraph run info. Partially addressed: `analyze addsplit` now refuses (unless `--force`) to add a split variant whose explicitly-overridden `--database`/`--gtf` conflicts with the object's original build provenance, and `analyze merge` now warns if the two objects being merged have different `size_splits` variants present (which `ad.concat`'s `uns_merge='same'` would otherwise silently drop). Full conflicting-run-info enforcement across `analyze merge` is still open.
- **Trimmer**: Add log file output and quiet mode to `toolTrim.py`.
- **Trimmer UMI Handling**: Test that UMI handling in `toolTrim.py` works as expected.
- **logging**: Standardize logging across all scripts for better traceability. Improve log messages for clarity.
- **Warnings**: Address remaining `FutureWarnings` during AnnData creation.
- **Rename de_results**: Change `de_results` in `adata.uns` to `log2FC` for clarity.
- **Replace log2FC with deseq implementation**: Use DESeq2 results directly instead of calculating log2FC manually.
- **--lazy flag**: Decide to remove this flag since mapping is now seperate from processing.
- **--nofrag flag**: Currently unsure if this is needed as the plan is to replace fragment logic with a more general approach of spliting BAMs by size.
- **--nosizefactors flag**: Since both normalized and raw counts are stored, this flag may be redundant. Maybe implement into graphs instead to plot raw counts directly.
- **map makes extra directories**: both graph and results directories are created when only results is needed at this step.
- **-maponly flag**: Consider removing this flag since mapping is now seperate from processing.
- **tqdm**: implement this for progress updates.
- **Unit Tests**: Expand unit tests to cover more functions and edge cases.
- **dir Structure**: Standardize directory structure creation across all modules for consistency and update readme/docs accordingly.
- **heatmap.csv**: Redirect these to the results directory instead of graph directory.
- **validate if warnings in adataGraph are critical**: Currently suppressed, need to ensure they are not hiding important issues.
- **split files have small variances**: Currently unknown cause, need to investigate, very very small variances are not expected with raw layers.
