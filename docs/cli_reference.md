# Command Line Interface (CLI) Reference

This document provides a detailed reference for all command-line commands and flags available in tRNAgraph.

## Global Options

These options apply to most commands in the toolkit.

| Flag                  | Description                                                                         | Default |
| :-------------------- | :---------------------------------------------------------------------------------- | :------ |
| `-q`, `--quiet`       | Suppresses console output. A run's log is always persisted regardless -- see below. | `False` |
| `-v`, `--verbose`     | Enables detailed execution logs.                                                    | `False` |
| `--skip-env-check`    | Skips the environment validation checks (dependencies and versions).                | `False` |
| `--skip-update-check` | Skips the background check for a newer tRNAgraph release (see `update` below).      | `False` |

Every command except `tools test` (which keeps its own fixed `toolsTestSuite.log`, overwritten each run) always writes a timestamped log under `./.log/` (e.g. `.log/20260101_120000_trim.log`), unconditionally -- even under `--quiet`, which only suppresses the console, never the file. There is no `--log <path>` flag to redirect this; on success the log moves into the command's real output directory (next to whatever it produced); on a crash or premature exit, a warning is printed pointing at the log still sitting in `.log/`, and it's left there rather than moved to a possibly-incomplete destination. `.log/` is untracked by git.

---

## Preprocess Modules

### `makedb`

Builds a Bowtie2 index from gtRNAdb/tRNAScan-SE output and a reference genome. This index is required for the map step.

**Usage:**

```bash
trnagraph preprocess makedb -g <genome> -t <trnaout> -r <trnafa> -m <namemap> [options]
```

**Flags:**

- **`-g`, `--genome`** (Required): Path to the reference genome FASTA file. This is used to extract the exact sequences of the tRNAs including flanking regions.
- **`-t`, `--trnaout`** (Required): Path to the tRNAScan-SE output file. This file contains the genomic coordinates of predicted tRNAs.
- **`-r`, `--trnafa`** (Required): Path to the tRNA reference FASTA file (usually from gtRNAdb).
- **`-m`, `--namemap`** (Required): Path to the tRNA name mapping file. This maps the tRNAScan-SE IDs to standard tRNA names.
- **`-s`, `--orgmode`**: Organism mode used for tRNAScan-SE. This ensures the correct scoring model and structure are used. Organism mode (`euk`, `bact`, `arch`). Default: `euk`.
- **`-o`, `--output`**: Output prefix. Default: `db`.
- **`--addtrna`**: Path to a FASTA file containing additional tRNA sequences to include in the database (e.g., mitochondrial tRNAs if not in the main set).
- **`--addseqs`**: Path to a FASTA file containing additional non-tRNA sequences to include (e.g., spike-ins, rRNA).
- **`--forcecca`**: Forces the addition of a CCA tail to the 3' end of all tRNA sequences in the database.
- **`-n`, `--threads`**: Number of threads to use for building the index. Default: max_cores.

### `trim`

Trims adapters and processes UMIs using `fastp`.

**Usage:**

```bash
trnagraph preprocess trim -i <manifest> [options]
```

**Flags:**

- **`-i`, `--manifest`** (Required): Tab-delimited file (`OutputPrefix <tab> R1 <tab> R2`). `OutputPrefix` doubles as the sample name and where trimmed output is written — a bare name writes to `processed/trimmed/<name>_trimmed.fastq.gz` (or `_merged.fastq.gz` for paired-end); a name containing a directory writes there instead.
- **`-a1`, `--adapter1`**: Adapter sequence for R1 (Auto-detect if omitted).
- **`-a2`, `--adapter2`**: Adapter sequence for R2 (for paired-end data, Auto-detect if omitted).
- **`-l`, `--length`**: Minimum sequence length allowed after trimming. Reads shorter than this will be discarded. Default: `15`.
- **`-u`, `--umilength`**: Length of the Unique Molecular Identifier (UMI) in base pairs. Set to `0` to disable UMI extraction. Default: `0`.
- **`--umi3`**: Specifies that the UMI is located at the 3' end of the read. If not set, it is assumed to be at the 5' end. Default: `False` (5' end).
- **`-n`, `--threads`**: Number of threads to use for fastp. Default: max_cores.
- **`--colormap`**: Path to a JSON file specifying custom colors for the trim-stats plot's read-type bars, under a top-level `trimtype` key (see [Custom Colormaps](advanced_usage.md#custom-colormaps---colormap)). Falls back to the default palette if omitted or if the file has no `trimtype` key. Recognized bar categories: `Merged`/`Unmerged` (paired-end samples only -- fastp's merge step doesn't run on single-end input), `Trimmed` (single-end samples' filter-passing reads), `Discarded` (either type).

> [!NOTE]
> If R2_Path is omitted in the manifest file, the sample will be treated as single-end.

### `map`

Maps reads to the tRNA database using Bowtie2.

**Usage:**

```bash
trnagraph preprocess map -i <metadata> -d <database> -o <output> [options]
```

**Flags:**

- **`-i`, `--input`** (Required): Path to the metadata file.
- **`-d`, `--database`** (Required): Path/Name of the Bowtie2 index (tRNA database) created by `makedb`.
- **`-o`, `--output`** (Required): Experiment name. This is used to name the output directory and report files.
- **`--bamdir`**: Directory where BAM files will be stored. Default: `processed/bam`
- **`--force-remap`**: Force remapping even if a matching bam file already exists. Default (when omitted): skip mapping if a bam already exists (after a fastq/header consistency check).
- **`--local`**: Use Bowtie2 local alignment mode instead of end-to-end. Default: `False`.
- **`--minnontrnasize`**: Minimum read length required for a read to be assigned to a non-tRNA feature. Default: `20`
- **`--skipcheck`**: Skips the validation check that ensures FASTQ read names match BAM headers. Use with caution. Default: `False`
- **`-n`, `--threads`**: Number of threads to use for Bowtie2 mapping. Default: `8`

> [!IMPORTANT]
> The number of threads used for mapping is highly system dependent. Bowtie2 can be memory intensive, and using too many threads can cause the system to run out of memory or lose performance because of overhead. It is recommended to use a number of threads that is appropriate for your system's available memory and CPU cores although between 8-10 has been commonly used on high performance machines as an optimal range and starting point.

---

## Analyze Modules

### `build`

Constructs the AnnData object from mapping outputs.

**Usage:**

```bash
trnagraph analyze build -i <metadata> -d <database> -o <out_dir> [options]
```

**Flags:**

- **`-i`, `--input`** (Required): Path to the metadata file.
- **`-d`, `--database`** (Required): Path/Name of the Bowtie2 index (tRNA database) created by `makedb`.
- **`-o`, `--output`** (Required): Output directory. The AnnData file will be created inside this directory with the name `<basename>.h5ad` alongside a results/ directory with additional output files.
- **`--gtf`**: Path to an Ensembl GTF file. Used to annotate non-tRNA features.
- **`--pairs`**: File listing sample pairs for direct comparison.
- **`--bed`**: List of additional BED files to define custom features.
- **`--maxmismatches`**: Maximum mismatches allowed per read. This is one consistent read-level filter applied identically everywhere reads get counted from BAM files -- affects the X matrix, uns aggregate counts (type/amino/anticodon), and coverage data the same way, not a different mismatch-filtered subset for different outputs. Default: `None` (No limit)
- **`--minfeaturereads`**: Minimum total raw read count a tRNA gene needs to count toward the VST fit (see [Data Structure: Layers](data_structure.md#layers) for exactly what this does and doesn't affect). Default: `30`
- **`--minnontrnasize`**: Minimum read length for non-tRNA features. Default: `20`
- **`--hub`**: Generates a UCSC Genome Browser track hub for the data. Default: `False`
- **`--hubonly`**: Only generates the track hub and skips AnnData generation. Default: `False`
- **`--filterother`**: Dumps reads counted in the `other` type category (the `other` row in `typecounts.txt`/`adata.uns['type_counts']`) to a separate BAM file for manual inspection -- i.e. reads matching no tRNA, bed, or GTF-annotated feature (being classified `other` and matching no other feature are the same condition; there's no separate/narrower `other` tag). Default: `False`
- **`--bamdir`**: Custom directory to look for BAM files if they are not in the default location. Default: `processed/bam`
- **`--filtermultimapped`**: Drops genomically multi-mapped reads (a read aligning to more than one location in the genome) from the entire coverage build before any column is computed. Unrelated to the separate, always-computed tRNA-identity uniqueness in `results/unique/`/`graphs/unique/`, which this flag does not affect. Default: `False`
- **`--dispfittype`**: Dispersion fit type for DESeq2. parametric`(default, standard for DESeq2),`mean`(robust for small sample sizes). Default:`parametric`
- **`-c`, `--readlengthsplit`**: Read length cutoff for splitting. When specified, adds `u<N>`/`o<N>` split variants to the _same_ output `.h5ad` as additional layers/obsm/uns entries (see [Data Structure: Split Variants](data_structure.md#split-variants---readlengthsplit)) — no separate `_uN.h5ad`/`_oN.h5ad` files are written. Further cutoffs can be added later to an existing object via `analyze addsplit`. Default: `None` (disabled)
- **`--overwritebams`**: Force overwrite of existing BAM files during map/split. Default: `False`
- **`--savesplitbams`**: The `u<N>`/`o<N>` split BAM files generated by `--readlengthsplit` (under `--bamdir`) are scratch files used only to compute that variant's coverage/counts — by default they're deleted once merged into the AnnData object. Pass this flag to keep them on disk instead. Default: `False`
- **`--vst`**: Sets the Variance Stabilizing Transformation (VST) computation method. Options: `vst` (native PyDESeq2 VST), `log1p` (np.log1p + StandardScaler, Default), or `none` (disable VST computation).
- **`-n`, `--threads`**: Number of threads to use for processing. Default: `8`

> [!NOTE]
> PyDeseq2 will automatically switch to `mean` dispersion fitting if the number of samples is too small for `parametric` to be stable.

> [!NOTE]
> DESeq2 size factors are always computed twice for the main feature matrix: once using only tRNA/tRX features as the normalization reference (the **default**, backing `adata.X`, `adata.layers["raw"]`, and `adata.obs["deseq2_sizefactor"]`), and once using all features as the reference (a secondary set kept for comparison). Both are stored in `adata.uns["deseq2_sizefactors_trna"]` / `adata.uns["deseq2_sizefactors_allfeatures"]`, and the all-feature-normalized layer is available as `adata.layers["norm_allfeatures"]`. If a GTF file is not provided no non-tRNA features are counted and the two sizefactor sets should be identical.

### `addsplit`

Adds an additional read-length split variant (an under/over cutoff pair) to an _existing_ AnnData object, without disturbing any variant already present — e.g. build with `-c 60`, then later add `-c 50` to see `u50`/`o50` alongside the existing `u60`/`o60`. Uses the same computation as `build`'s `--readlengthsplit`, so both paths produce identical results for the same cutoff/data. See [Data Structure: Split Variants](data_structure.md#split-variants---readlengthsplit).

**Usage:**

```bash
trnagraph analyze addsplit -i <input.h5ad> -c <cutoff> [options]
```

**Flags:**

- **`-i`, `--anndata`** (Required): Path to the existing AnnData file to add a split variant to.
- **`-c`, `--readlengthsplit`** (Required): New read length cutoff to add (generates `u<N>`/`o<N>` variants).
- **`--metadata`**: Metadata file. Default: recovered from provenance (see below).
- **`--bamdir`**: Original (unsplit) BAM directory. Default: recovered from provenance.
- **`-d`, `--database`**: Override tRNA database. Default: recovered from provenance.
- **`--gtf`**: Override GTF path. Default: recovered from provenance.
- **`--dispfittype`**: Override DESeq2 dispersion fit type. Default: recovered from provenance.
- **`--vst`**: VST strategy for this split's `vst` layer. Default: recovered from provenance.
- **`--minfeaturereads`**: Override the minimum total raw read count a tRNA gene needs for this split's VST fit. Default: recovered from provenance, or `30`.
- **`--overwritebams`**: Force overwrite of existing split BAM files. Default: `False`.
- **`--savesplitbams`**: Keep the `u<N>`/`o<N>` split BAM files (under `--bamdir`) instead of deleting them once merged into the AnnData object. Default: `False`.
- **`-o`, `--output`**: Output `.h5ad` path. Default: overwrite the input file in place.
- **`-w`, `--overwrite`**: Overwrite this cutoff's `u<N>`/`o<N>` data if already present in the object. Default: `False`.
- **`--force`**: Proceed even if an explicitly-overridden `--database`/`--gtf` conflicts with the object's original build provenance. Default: `False`.
- **`-n`, `--threads`**: Number of threads to use. Default: `8`.

> [!NOTE]
> "Recovered from provenance" means read back from `trnagraphruninfo`, the build history tRNAgraph stores inside the h5ad itself (`adata.uns["trnagraphruninfo"]["flags"]`) — see [Data Structure: Run Information](data_structure.md#run-information). Any flag above falls back to that recorded value unless you pass it explicitly here.

> [!NOTE]
> Each split variant gets its own independently-fit DESeq2 size factors, computed from that variant's own BAMs — never a shared/global normalization.

### `cluster`

Performs clustering (UMAP, HDBSCAN) on an existing AnnData object.

> [!NOTE]
> The database object can be clustered using the `cluster` function. The following code will cluster the database object using [UMAP](https://umap-learn.readthedocs.io/en/latest/index.html) and cluster using [HDBSCAN](https://hdbscan.readthedocs.io/en/latest/index.html). The default parameters used in tRNAgraph work well on ARMseq, DM-tRNAseq, and OTTRseq data; however, each dataset is different and may require fine-tuning of the parameters to yield the best results. The following code will cluster the database object and save the result to a new database object.

**Usage:**

```bash
trnagraph analyze cluster -i <input.h5ad> [options]
```

**Flags:**

- **`-i`, `--input`** (Required): Path to the anndata file to be clustered.
- **`-o`, `--output`**: Path to save the clustered anndata file. Default: `trnagraph.cluster.h5ad`.
- **`-w`, `--overwrite`**: Overwrites existing cluster information in the file if present. Default: `False`.
- **`-r`, `--randomstate`**: Seed for random number generation to ensure reproducible UMAP results. Default: `None` (Random)
- **`-t`, `--readcutoff`**: Minimum read count for a tRNA to be included in clustering. Default: `20`
- **`--coveragetype`**: List of coverage features to use for clustering. Default: `['uniquecoverage', 'readstarts', 'readends', 'mismatchedbases', 'deletions']`
- **`--variant`**: Select which `<norm>:<tag>` variant to cluster, e.g. `norm:u60` (see [Data Structure: Split Variants](data_structure.md#split-variants---readlengthsplit) for the syntax). Default: `norm:full`. Cluster results for a split variant are stored namespaced under that variant (`adata.uns['size_splits'][tag]`/`adata.obsm['size_split_<tag>']`) so they never overwrite another variant's results.
- **UMAP Parameters:**
  - **`-c1`, `--ncomponentsmp`**: Components for sample clustering (Default: 2).
  - **`-c2`, `--ncomponentgrp`**: Components for group clustering (Default: 2).
  - **`-l1`, `--neighborclusmp`**: Neighbors for sample clustering (Default: 150).
  - **`-l2`, `--neighborclusgrp`**: Neighbors for group clustering (Default: 40).
  - **`-m`, `--mindist`**: Minimum distance for UMAP (Default: 0.1).
- **HDBSCAN Parameters:**
  - **`-d1`, `--hdbscanminsampsmp`**: Min samples for sample clustering (Default: 6).
  - **`-d2`, `--hdbscanminsampgrp`**: Min samples for group clustering (Default: 3).
  - **`-b1`, `--hdbscanminclusmp`**: Min cluster size for sample clustering (Default: 30).
  - **`-b2`, `--hdbscanminclugrp`**: Min cluster size for group clustering (Default: 10).

> [!IMPORTANT]
> Clustering is performed across the `uniquecoverage`, `readstarts`, `readends`, `mismatchedbases`, and `deletions` categories of the AnnData object. When performing clustering, verifying that it is reproducible and that the results reflect the data is important. This can be done by running the clustering multiple times and comparing the results. The clustering is also performed on `sample` and `group` observations. In the case of samples, every set of reads for every single tRNA is used for clustering. In the case of groups, the mean of the reads is taken for each tRNA across the read categories and then used for clustering. This is done to reduce the number of samples used for clustering and to reduce the noise in the clustering. The results will be saved in the `obs` attribute of the database object as `sample_cluster\umap1\umap2` and `group_cluster\umap1\umap2`, respectively. Clusters annotated as `-1` are considered noise and are not included in the clustering. Plotting of the clustering is done as well for convenience.

> [!WARNING]
> When working with downstream analysis of the cluster groups, it is important to note that reads that are dropped via the `--readcutoff` flag will not be included in the clustering however, they are still present in the AnnData object. This means that your object can contain NaN values in the clustering columns. Depending on your use case, you may want to filter these out before performing any analysis.

### `graph`

Generates a wide variety of visualizations from the AnnData file.

**Usage:**

```bash
trnagraph graph -i <input.h5ad> -o <output_dir> [options]
```

**General Flags:**

- **`-i`, `--input`** (Required): Input AnnData file.
- **`-o`, `--output`**: Output directory. Default: `figures`.
- **`-g`, `--graphtypes`**: List of graphs to generate (`all`, `cluster`, `correlation`, `count`, `coverage`, `heatmap`, `logo`, `pca`, `radar`, `volcano`). Default: `all`.
- **`--config`**: JSON configuration file for filtering.
- **`--colormap`**: JSON file for custom colors.
- **`--regen_uns`**: Force regeneration of calculated stats.
- **`--variant`**: Select which `<norm>:<tag>` variant to plot, e.g. `raw:full`, `norm:u60`, `allfeatures:o60` (see [Data Structure: Split Variants](data_structure.md#split-variants---readlengthsplit) for the syntax). Default: `norm:full`. Plots for a non-default variant are written to a `<norm>_<tag>/` subfolder under each graph type's output directory, so different `--variant` runs into the same `-o` don't overwrite each other.

**Cluster Plot Options:**

- **`--clustergrp`**: Grouping variable. Default: `amino`.
- **`--clusterlabels`**: Custom labels for clusters.
- **`--clusteroverview`**: Generate overview plot. Default: `False`.
- **`--clustermask`**: Mask unclustered points. Default: `False`.

**Coverage Plot Options:**

- **`--covgrp`**: Grouping variable. Default: `group`.
- **`--covtype`**: Coverage type (e.g., `uniquecoverage`).
- **`--covmethod`**: Combination method (`mean`).
- **`--combinedpdfonly`**: Skip individual tRNA PDFs. Default: `False`.

**Heatmap Options:**

- **`--heatgrp`**: Grouping variable. Default: `group`.
- **`--diffrts`**: Read types for differential analysis (shared with volcano).
- **`--heatcutoff`**: Read count cutoff. Default: `80`.
- **`--heatsubplots`**: Also save each individual comparison's heatmap as its own PDF, in an `individual/` subfolder next to the combined multi-page PDFs (which are unaffected). Default: `False`.

**Volcano Options:**

- **`--volgrp`**: Grouping variable used both to define the pairwise group comparisons and to look up per-group colors in `--colormap` (the same `<obs_column>: {<value>: <color>}` shape used by PCA, keyed on `--volgrp`'s value, e.g. `"group"`). Default: `group`.
- **`--diffrts`**: Read types to generate per-readtype tRNA volcano plots for, shared with heatmap.
- **`--volcutoff`**: Read count cutoff. Default: `80`.
- **`--vollabels`**: Number of top significant markers to label on each plot, ranked by `|log2FC| * -log10(p-value)`. Omit to label every significant marker (default); pass `0` to disable labels entirely.

> [!NOTE]
> Two extra subplots are generated automatically (in combined plots) whenever `adata.uns['nontRNA_counts']` is non-empty (i.e., `--gtf` was used at `analyze build`): a non-tRNA-only plot and a combined tRNA + non-tRNA plot. No additional flags are needed, and both are skipped with a log message if non-tRNA counts are unavailable. These use a different DESeq2 normalization than the per-readtype plots — see [Data Structure: Graphing Notes](data_structure.md#6-graphing-notes).

**PCA Options:**

- **`--pcamarkers`**: Marker style grouping. Default: `sample`.
- **`--pcacolors`**: Color grouping. Default: `group`.
- **`--pcareadtypes`**: Read types to generate per-readtype tRNA PCA plots for (`tRNA_<pcamarkers>_by_<pcacolors>_<readtype>_*`). Default: `total_unique`, `total`.

> [!NOTE]
> Two extra plots are generated automatically whenever `adata.uns['nontRNA_counts']` is non-empty (i.e., `--gtf` was used at `analyze build`): a non-tRNA-only plot and a combined tRNA + non-tRNA plot. No additional flags are needed, and both are skipped with a log message if non-tRNA counts are unavailable. These use a different DESeq2 normalization than the per-readtype plots — see [Data Structure: Graphing Notes](data_structure.md#6-graphing-notes).

**Radar Options:**

- **`--radargrp`**: Grouping variable. Default: `group`.
- **`--radarscaled`**: Scale axes to 100%. Default: `False`.

**SeqLogo Options:**

- **`--logogrp`**: Grouping variable. Default: `amino`.
- **`--logopseudocount`**: Pseudocount added. Default: `20`.
- **`--logosize`**: Sequence size preset. Default: `noloop`.
- **`--ccatail`**: Keep CCA tail. Default: `True`.
- **`--pseudogenes`**: Keep pseudo-tRNAs. Default: `True`.

---

## Tools Commands

The `tools` module contains utility functions for specific analysis tasks or testing.

### `log2fc`

Computes log2 fold change for specified groups and read types.

**Usage:**

```bash
trnagraph tools log2fc -i <input.h5ad> [options]
```

- **`-i`, `--anndata`** (Required): Input file.
- **`-g`, `--group`**: Grouping variable from obs. Default: `group`.
- **`-r`, `--readtypes`**: List of read types to analyze.
- **`-x`, `--cutoff`**: Read count cutoff(s). Default: `80`.
- **`-c`, `--config`**: Config file for filtering.
- **`--variant`**: Select which `<norm>:<tag>` variant to compute log2fc for, e.g. `norm:u60` (see [Data Structure: Split Variants](data_structure.md#split-variants---readlengthsplit) for the syntax). Default: `norm:full`. Results for a split variant are written to `adata.uns['size_splits'][tag]['log2FC']` instead of the top-level `adata.uns['log2FC']`.

### `csv`

Exports the AnnData object to a set of CSV files (obs, var, X).

**Usage:**

```bash
trnagraph tools csv -i <input.h5ad> -o <output_dir>
```

- **`-i`, `--anndata`** (Required): Input file.
- **`-o`, `--output`**: Output directory. Default: `csv`.

### `merge`

Merges two existing AnnData objects into a single file. This is useful for combining datasets processed separately.

**Usage:**

```bash
trnagraph tools merge -i1 <file1.h5ad> -i2 <file2.h5ad> [options]
```

**Flags:**

- **`-i1`, `--anndata1`** (Required): Path to the first AnnData file.
- **`-i2`, `--anndata2`** (Required): Path to the second AnnData file.
- **`-o`, `--output`**: Path for the merged output file. Default: `trnagraph.merge.h5ad`.
- **`--dropno`**: Drop non-tRNA genes that are not present in both objects. Default: `False`.
- **`--droprna`**: Drop RNA categories that are not present in both objects. Default: `False`.
- **`--force`**: Proceed even if the two objects' build provenance (database/gtf) conflicts. Default: `False` (refuses the merge on a conflict).

### `test`

Runs the internal automated test suite.

**Usage:**

```bash
trnagraph tools test [options]
```

**General Options:**

- **`--all`**: Run all tests including split analysis, forcing a clean workspace and full redownload.
- **`--skip-download`**: Skip the metadata/fastq/tRNA/genome download steps and run everything else. Downloads are already skipped by default when the target files are already present; this forces the skip regardless.
- **`--cleanrun`**: Clean up test files after completion.
- **`-d`, `--directory`**: Directory to run tests in.

**Step Flags (run specific steps):**

| Flag            | Description                                                    |
| :-------------- | :------------------------------------------------------------- |
| `--metadata`    | Download metadata                                              |
| `--fastq`       | Download FASTQ files                                           |
| `--trna`        | Download tRNA sequences                                        |
| `--genome`      | Download reference genome                                      |
| `--trim`        | Run adapter trimming                                           |
| `--makedb`      | Create tRNA database                                           |
| `--map`         | Run read mapping                                               |
| `--build`       | Build AnnData object (no split analysis)                       |
| `--split-build` | Build AnnData with read length split (generates u60/o60 files) |
| `--cluster`     | Run clustering algorithms                                      |
| `--graph`       | Generate visualization plots (main h5ad only)                  |
| `--split-graph` | Generate plots for split h5ad files                            |
| `--hubonly`     | Generate UCSC track hubs only                                  |

---

## Update

### `update`

Updates this git checkout to the latest source and re-syncs the local environment. Refuses to run if the checkout has uncommitted local changes to tracked files, so nothing is at risk of being lost or silently merged over.

**Usage:**

```bash
trnagraph update [options]
```

**Options:**

- **`--branch <name>`**: Update to this branch instead of `main` (e.g. `dev`). Mutually exclusive with `--tag`. Refused if the branch's own reported version is older than what's currently installed (see below).
- **`--tag <version>`**: Check out this release tag instead of a branch (e.g. `v1.9.0` or `1.9.0`). This leaves the checkout in git's standard "detached HEAD" state for a tag checkout -- a printed message explains what that means and how to get back to a branch afterward. Mutually exclusive with `--branch`. Refused if the requested version is older than v1.9.0, the version `update` was itself introduced in.
