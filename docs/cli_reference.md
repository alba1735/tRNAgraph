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

- **`--config`**: JSON configuration file whose `flags.makedb` block pins this command's options, so one file can carry a whole run. A flag typed on the command line always beats the file. See [Run Configuration](advanced_usage.md#run-configuration---config).

### `trim`

Trims adapters and processes UMIs using `fastp`.

**Usage:**

```bash
trnagraph preprocess trim -i <manifest> [options]
```

**Flags:**

- **`-i`, `--input`** (Required): Tab-delimited manifest file (`OutputPrefix <tab> R1 <tab> R2`). `OutputPrefix` doubles as the sample name and where trimmed output is written — a bare name writes to `processed/trimmed/<name>_trimmed.fastq.gz` (or `_merged.fastq.gz` for paired-end); a name containing a directory writes there instead.
- **`-a1`, `--adapter1`**: Adapter sequence for R1 (Auto-detect if omitted).
- **`-a2`, `--adapter2`**: Adapter sequence for R2 (for paired-end data, Auto-detect if omitted).
- **`-l`, `--length`**: Minimum sequence length allowed after trimming. Reads shorter than this will be discarded. Default: `15`.
- **`-u`, `--umilength`**: Length of the Unique Molecular Identifier (UMI) in base pairs. Set to `0` to disable UMI extraction. Default: `0`. The UMI is moved into the read name, separated by an underscore (`@READ_AAAACAC`) — the same shape `umi_tools` produces, so a dataset trimmed here can be deduplicated later with [`preprocess map --dedup`](#map) whichever UMI path it took.
- **`--umi3`**: Specifies that the UMI is located at the 3' end of the read. If not set, it is assumed to be at the 5' end. Default: `False` (5' end). fastp has no tail-anchored UMI option, so this case runs fastp with no UMI flags and applies `umi_tools extract --3prime` to its output afterward.
- **`-n`, `--threads`**: Number of threads to use for fastp. Default: max_cores.
- **`--style`**: Path to a JSON style file. Only its `colors.trimtype` block is read here (see [Style Files](advanced_usage.md#style-files---style)), but it is the same file `analyze graph` takes, so one file can style the whole pipeline. Falls back to the default palette if omitted or if the file has no `trimtype` key. Recognized bar categories: `Merged`/`Unmerged` (paired-end samples only -- fastp's merge step doesn't run on single-end input), `Trimmed` (single-end samples' filter-passing reads), `Discarded` (either type).

> [!NOTE]
> If R2_Path is omitted in the manifest file, the sample will be treated as single-end.

- **`--config`**: JSON configuration file whose `flags.trim` block pins this command's options, so one file can carry a whole run. A flag typed on the command line always beats the file. See [Run Configuration](advanced_usage.md#run-configuration---config).

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
- **`--dedup`**: Deduplicate mapped reads by UMI using `umi_tools dedup`, as a separate phase after mapping completes. The deduplicated reads take over the sample's ordinary `<sample>.bam` name, so `analyze build`/`--bamdir` need no knowledge that it happened. Default: `False`.
- **`--keep-prededup`**: Keep the pre-deduplication BAM as `<sample>.prededup.bam` instead of discarding it. Useful for comparing deduplicated against non-deduplicated output without paying for a second mapping run. Default: `False`.
- **`--dedup-method`**: `umi_tools dedup --method` to use: `unique`, `percentile`, `cluster`, `adjacency` or `directional`. Default: `directional`. Note that `directional` merges UMIs within edit distance 1, which becomes progressively more aggressive as per-position UMI density rises; on a deeply sequenced library with a short UMI it can compress counts for the most abundant features substantially. `unique` does no such merging and is the comparison worth running if that matters for your data.
- **`-n`, `--threads`**: Number of threads to use for Bowtie2 mapping. Default: `8`

> [!IMPORTANT]
> `--dedup` **refuses to run** if it cannot find a UMI in the BAM's read names, rather than warning and continuing. Given UMI-less reads, `umi_tools dedup` does not fail — it falls back to collapsing reads that share an alignment position, and for short, deeply-covered tRNA transcripts those reads are overwhelmingly genuine molecules rather than PCR duplicates. That would delete real signal and leave nothing in the output to show it had happened. Extract UMIs at trimming time with [`preprocess trim -u/--umilength`](#trim) first.

> [!NOTE]
> The separator between the read name and its UMI is detected from the BAM rather than assumed, because `umi_tools dedup` defaults to `_` and would silently mis-parse anything else. tRNAgraph's own trimming produces `_` on both UMI paths, but externally-trimmed BAMs vary.

> [!NOTE]
> Deduplication runs samples concurrently. The job count is the smallest of half `--threads`, what total system memory affords, and the number of samples. Memory is the binding constraint rather than cores: `umi_tools` is single-threaded, so samples are the only unit of parallelism, but each job holds a position's reads and their UMI network in memory. The per-job estimate is eight times the largest BAM's on-disk size, since real usage scales with read count and UMI length rather than being a fixed figure. On a nine-sample human dataset this takes the phase from ~2.4 hours serial to roughly the length of its slowest single sample.

> [!NOTE]
> These flags are a convenience wrapper, not a replacement for `umi_tools`. It is already installed in the project's conda environment, so for a protocol this wrapper does not cover — paired UMIs, a cell barcode, a non-standard read-name layout, or any other `umi_tools dedup` option — run `umi_tools` directly against the BAM directory and then continue with `analyze build` as normal. A per-sample `umi_tools` log is written next to each BAM, and `results/<exp>-dedupinfo.txt` records the method, the detected separator and the `umi_tools` version used.

> [!TIP]
> `results/<exp>-dedupstats.txt` is the file to read to decide whether deduplication behaved, and whether a sample should be dropped. Per sample: `input_reads`, `output_reads`, `retained_pct`, `reads_per_molecule`, `positions`, `reads_per_position`, `mean_umis_per_position`, `max_umis_per_position`. Every value comes from `umi_tools`' own end-of-run log, so it costs a text parse rather than another pass over the BAMs.
>
> How to read it:
>
> - **`reads_per_molecule`** is the duplication level. A sample well above its group-mates was amplified harder, and its raw counts carry proportionally more PCR noise.
> - **`positions` and `reads_per_position`** together separate a deeply-sequenced library from a low-complexity one. Few positions with low reads each means the library itself was limited, which deduplication cannot fix and more sequencing will not either.
> - **`max_umis_per_position`** read against your UMI's own ceiling (4ⁿ for an n-base UMI; 16,384 for 7 nt) says whether the tag space is saturating. As it approaches that ceiling, deduplicated counts at the deepest features measure UMI diversity rather than molecule count, and dynamic range compresses.
>
> Pair it with `results/<exp>-replicatecorrelation.txt` from `analyze build`: deduplication should _widen_ the gap between within-group and between-group correlation, because it removes a technical artifact. A run where that gap narrows is a warning that it removed signal instead.

> [!IMPORTANT]
> The number of threads used for mapping is highly system dependent. Bowtie2 can be memory intensive, and using too many threads can cause the system to run out of memory or lose performance because of overhead. It is recommended to use a number of threads that is appropriate for your system's available memory and CPU cores although between 8-10 has been commonly used on high performance machines as an optimal range and starting point.

---

- **`--config`**: JSON configuration file whose `flags.map` block pins this command's options, so one file can carry a whole run. A flag typed on the command line always beats the file. See [Run Configuration](advanced_usage.md#run-configuration---config).

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
- **`--dispfittype`**: Dispersion fit type for DESeq2. parametric`(default, standard for DESeq2),`mean`(robust for small sample sizes). Default:`parametric`
- **`-c`, `--readlengthsplit`**: Read length cutoff for splitting. When specified, adds `u<N>`/`o<N>` split variants to the _same_ output `.h5ad` as additional layers/obsm/uns entries (see [Data Structure: Split Variants](data_structure.md#split-variants---readlengthsplit)) — no separate `_uN.h5ad`/`_oN.h5ad` files are written. Further cutoffs can be added later to an existing object via `analyze addsplit`. Default: `None` (disabled)
- **`--overwritebams`**: Force overwrite of existing BAM files during map/split. Default: `False`
- **`--savesplitbams`**: The `u<N>`/`o<N>` split BAM files generated by `--readlengthsplit` (under `--bamdir`) are scratch files used only to compute that variant's coverage/counts — by default they're deleted once merged into the AnnData object. Pass this flag to keep them on disk instead. Default: `False`
- **`--vst`**: Sets the Variance Stabilizing Transformation (VST) computation method. Options: `vst` (native PyDESeq2 VST), `log1p` (np.log1p + StandardScaler, Default), or `none` (disable VST computation).
- **`-n`, `--threads`**: Number of threads to use for processing. Default: `8`

> [!TIP]
> Every build writes `results/<exp>-replicatecorrelation.txt`, a quality-control table asking whether samples agree with their own group more than with other groups. The `separation` value at the top (mean within-group r² minus mean between-group r²) says whether your grouping is visible in the data at all; the per-sample rows identify which sample disagrees with its replicates if one does. A `n_within_pairs` of 0 means every sample is its own group -- see the `trim_metadata.tsv` warning under [`trim`](#trim). The same tables are stored in `adata.uns['replicate_correlation*']`.

> [!NOTE]
> PyDeseq2 will automatically switch to `mean` dispersion fitting if the number of samples is too small for `parametric` to be stable.

> [!NOTE]
> DESeq2 size factors are always computed twice for the main feature matrix: once using only tRNA/tRX features as the normalization reference (the **default**, backing `adata.X`, `adata.layers["raw"]`, and `adata.obs["deseq2_sizefactor"]`), and once using all features as the reference (a secondary set kept for comparison). Both are stored in `adata.uns["deseq2_sizefactors_trna"]` / `adata.uns["deseq2_sizefactors_allfeatures"]`, and the all-feature-normalized layer is available as `adata.layers["norm_allfeatures"]`. If a GTF file is not provided no non-tRNA features are counted and the two sizefactor sets should be identical. This applies to the complete variant only: a read-length split variant excludes non-tRNA features entirely, so it runs no all-feature pass and `--variant allfeatures:<tag>` is not available for it.

- **`--config`**: JSON configuration file whose `flags.build` block pins this command's options, so one file can carry a whole run. A flag typed on the command line always beats the file. See [Run Configuration](advanced_usage.md#run-configuration---config).

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

- **`--config`**: JSON configuration file whose `flags.addsplit` block pins this command's options, so one file can carry a whole run. A flag typed on the command line always beats the file. See [Run Configuration](advanced_usage.md#run-configuration---config).

### `order`

Writes a `--config` file's declared category order into an _existing_ AnnData object, as an ordered categorical on `adata.obs`. Order is otherwise alphabetical everywhere, which is wrong for any set of labels whose meaningful sequence is not its alphabet — `D0`, `D7`, `D14`, `D35` sorts to `D0`, `D14`, `D35`, `D7`.

The order has two consumers, which is why it is stored on the object rather than passed per call: plot legend and axis order, and the DESeq2 reference level, which is the **first** listed level of each column. `build` applies the same block automatically when given a `--config`; this command exists so an object built before the order was declared can gain it without a rebuild.

**Usage:**

```bash
trnagraph analyze order -i <input.h5ad> --config <config.json>
```

**Flags:**

- **`-i`, `--anndata`** (Required): Path to the existing AnnData file to order.
- **`--config`** (Required): JSON file whose top-level `order` block declares each column's category order.
- **`-o`, `--output`**: Output path. Default: overwrite the input file in place.
- **`-q`, `--quiet`**: Suppress output to stdout. Default: `False`.

**Config block:**

```json
{
  "name": "organoid_timecourse",
  "order": { "timepoint": ["Day 0", "Day 35", "Day 70"] }
}
```

`order` is a **top-level** key, a sibling of the `obs`/`var` filters rather than a member of `flags` — it describes what the experiment is, not how one command ran.

Every level present in the data must be listed. A level left out is an error naming both what is missing and what was declared, rather than being appended silently, since an incomplete declaration would no longer describe the column's order. A declared level that is absent from the data is kept, so a filtered subset does not change the ordering.

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
- **`-e`, `--variancethreshold`**: Variance threshold used for feature selection before clustering — features varying less than this across the dataset are dropped. Default: `0.1`.
- **`--clusterobsexperimental`**: **Experimental.** Names `adata.obs` columns to copy into `adata.var`/`adata.X` so they participate in clustering alongside the coverage features. Repeatable. Default: none.
- **`-n`, `--threads`**: Number of threads to use. Default: `0` (all available cores). Always passed to HDBSCAN's `core_dist_n_jobs`, and to UMAP's `n_jobs` only when no `--randomstate` seed is set — UMAP forces `n_jobs` to 1 when seeded, for reproducibility.
- **UMAP Parameters:**
  - **`-c1`, `--ncomponentsmp`**: Components for sample clustering (Default: 2).
  - **`-c2`, `--ncomponentgrp`**: Components for group clustering (Default: 2).
  - **`-l1`, `--neighborclusmp`**: Neighbors for sample clustering (Default: 150).
  - **`-l2`, `--neighborclusgrp`**: Neighbors for group clustering (Default: 40).
  - **`-n1`, `--neighborstdsmp`**: Neighbors used for the UMAP projection _plot_ of samples, separate from `-l1`, which governs the clustering itself (Default: 75).
  - **`-n2`, `--neighborstdgrp`**: Neighbors used for the UMAP projection _plot_ of groups, separate from `-l2` (Default: 20).
  - **`-m`, `--mindist`**: Minimum distance for UMAP (Default: 0.1).
  - **`-us`, `--umapstatsmetrics`**: Distance metric for UMAP (Default: `euclidean`).
- **HDBSCAN Parameters:**
  - **`-d1`, `--hdbscanminsampsmp`**: Min samples for sample clustering (Default: 6).
  - **`-d2`, `--hdbscanminsampgrp`**: Min samples for group clustering (Default: 3).
  - **`-b1`, `--hdbscanminclusmp`**: Min cluster size for sample clustering (Default: 30).
  - **`-b2`, `--hdbscanminclugrp`**: Min cluster size for group clustering (Default: 10).
  - **`-uh`, `--hdbstatsmetrics`**: Distance metric for HDBSCAN (Default: `euclidean`).

> [!IMPORTANT]
> Clustering is performed across the `uniquecoverage`, `readstarts`, `readends`, `mismatchedbases`, and `deletions` categories of the AnnData object. When performing clustering, verifying that it is reproducible and that the results reflect the data is important. This can be done by running the clustering multiple times and comparing the results. The clustering is also performed on `sample` and `group` observations. In the case of samples, every set of reads for every single tRNA is used for clustering. In the case of groups, the mean of the reads is taken for each tRNA across the read categories and then used for clustering. This is done to reduce the number of samples used for clustering and to reduce the noise in the clustering. The results will be saved in the `obs` attribute of the database object as `sample_cluster\umap1\umap2` and `group_cluster\umap1\umap2`, respectively. Clusters annotated as `-1` are considered noise and are not included in the clustering. Plotting of the clustering is done as well for convenience.

> [!WARNING]
> When working with downstream analysis of the cluster groups, it is important to note that reads that are dropped via the `--readcutoff` flag will not be included in the clustering however, they are still present in the AnnData object. This means that your object can contain NaN values in the clustering columns. Depending on your use case, you may want to filter these out before performing any analysis.

- **`--config`**: JSON configuration file whose `flags.cluster` block pins this command's options, so one file can carry a whole run. A flag typed on the command line always beats the file. See [Run Configuration](advanced_usage.md#run-configuration---config).

### `graph`

Generates a wide variety of visualizations from the AnnData file.

**Usage:**

```bash
trnagraph graph -i <input.h5ad> -o <output_dir> [options]
```

**General Flags:**

- **`-i`, `--input`** (Required): Input AnnData file.
- **`-o`, `--output`**: Output directory. Default: `figures`.
- **`-g`, `--graphtypes`**: List of graphs to generate (`all`, `cluster`, `compare`, `correlation`, `count`, `coverage`, `heatmap`, `logo`, `mismatch`, `pca`, `radar`, `venn`, `volcano`). Default: `all`. `all` now **unions** with anything else named, so `-g all -g venn` gives you both — it previously replaced the list and dropped the extra silently. Repeatable rather than space-separated: `-g volcano -g pca`, not `-g volcano pca` (the latter parses as one graph type plus two stray positional arguments and errors out).
- **`-n`, `--threads`**: Number of threads to use. Default: `0` (all available cores).
- **`--config`**: JSON configuration file for filtering, and optionally a `flags.graph` block pinning most `graph` options (grouping columns, cutoffs, readtypes, `--variant`, `--allreads`, ...). The same file carries a block per command, so one file can drive a whole run. A flag typed on the command line always beats the file. `--input`/`--output`/`--config`/`--style`/`--threads`/`--quiet`/`--verbose` are not settable, and `--format` belongs to `--style`. See [Run Configuration](advanced_usage.md#run-configuration---config).
- **`--style`**: JSON style file carrying the color palette and presentation settings (figure size, marker/font/line size, dpi, alpha, output format). See [Style Files](advanced_usage.md#style-files---style).
- **`--format`**: Output image format for every plot: `pdf`, `svg` or `png`. Overrides a `format` set in `--style`. Multi-page combined outputs stay PDF. Default: `pdf`.
- **`--regen_uns`**: Force regeneration of calculated stats.
- **`--variant`**: Select which `<norm>:<tag>` variant to plot, e.g. `raw:full`, `norm:u60`, `allfeatures:full` (see [Data Structure: Split Variants](data_structure.md#split-variants---readlengthsplit) for the syntax). Default: `norm:full`. Plots for a non-default variant are written to a `<norm>_<tag>/` subfolder under each graph type's output directory, so different `--variant` runs into the same `-o` don't overwrite each other.
- **`--allreads`**: Plot every graph type from all reads instead of unique reads. Default: `False` (unique). "Unique" here means **transcript-specific** — a read that aligned to exactly one mature tRNA transcript — not genome-level uniqueness; the genome MAPQ filter is separate and always applied. This is one option for the whole command, so two plots of the same dataset never rest on different denominators. Plots are written to an `allreads/` subfolder (nested inside the `--variant` folder when both are used) for every graph type **except coverage and PCA**, whose filenames already record the basis: coverage's carry `--covtype`'s category, and PCA's carry the readtype label (`total` against `total_unique`). Adding a directory segment there would write byte-identical files into a second place — PCA in particular always plots both bases on its comparative overview, so at default settings the two runs produced exactly the same files. Readtypes with no transcript-specific counts (`antisense`, and the pre-tRNA categories) fall back to all reads with a warning.

> [!NOTE]
> The PCA and volcano _combined overview_ pages always show both read bases side by side, whatever `--allreads` is set to. That is deliberate: it is the only place you can see how much transcript-level multi-mapping actually moves your data, and a labelled comparison is not the same thing as two plots silently disagreeing.

> [!NOTE]
> `compare` and `venn` are not part of `all`'s fixed list, but are folded in automatically when their prerequisites are satisfied — two different columns for `--comparegrp1`/`--comparegrp2`, and a `multivariate` block for `venn`. When a prerequisite is missing the type is left out and the run logs why, so a missing figure is never confused with an empty result. Naming one explicitly always includes it, so the run fails with an instruction rather than quietly skipping. `all` expands to `cluster`, `correlation`, `count`, `coverage`, `heatmap`, `logo`, `mismatch`, `pca`, `radar` and `volcano`; each of the other two is only produced when you ask for it by name (`-g compare`, `-g venn`). It cannot produce anything at default settings — it needs two _different_ `obs` columns, and the defaults leave `--comparegrp1` and `--comparegrp2` both set to `group` — so it depends on metadata beyond what a minimal experiment carries and is only meaningful when reached for on purpose. Including it in `all` would abort every ordinary run, since naming the plot is what makes the two flags' shared `group` default a contradiction. `venn` is excluded for a related but distinct reason: it needs a `multivariate` block in `--config` declaring the grouping and thresholds its sets are built from, and set-overlap figures produced from cutoffs nobody chose invite wrong conclusions. See **Compare Options** and **Venn Options** below.

> [!IMPORTANT]
> Every option below that names something inside the object -- a grouping column
> (`--covgrp`, `--volgrp`, `--heatgrp`, `--clustergrp`, `--corrgroup`, `--radargrp`,
> `--logogrp`, `--pcamarkers`, `--pcacolors`, `--covobs`, `--comparegrp1`/`--comparegrp2`,
> `--clusterlabels`), a coverage type (`--covtype`), or a readtype (`--diffrts`,
> `--pcareadtypes`) -- is checked against the object before any plotting begins. An
> unrecognised value **aborts the run**, reporting every bad label at once with its near
> matches. It is not substituted for a default: a mistyped grouping column would otherwise
> produce a complete, plausible set of figures grouped by a column you did not ask for, with
> nothing in the output to say so. Run [`trnagraph tools info`](#info) to list the valid
> values for any of them.

**Cluster Plot Options:**

- **`--clustergrp`**: Grouping variable. Default: `amino`.
- **`--clusterlabels`**: Custom labels for clusters.
- **`--clusteroverview`**: Generate overview plot. Default: `False`. Forced on when `-g all` is used.
- **`--clustermask`**: Mask unclustered points. Default: `False`.
- **`--clusternumeric`**: Treat the `--clustergrp` category as numeric rather than categorical, which changes how it is coloured and ordered. Default: `False`.

**Compare Options:**

Only used by `-g compare`, which `all` does not include.

- **`--comparegrp1`**: AnnData `obs` column drawn as the coloured **series** within each figure. Default: `group`.
- **`--comparegrp2`**: AnnData `obs` column the fold change is taken **between**. Default: `group`.

> [!IMPORTANT]
> The direction is the opposite of what the names suggest, and this reference previously described it backwards. The fold change is taken **between `--comparegrp2` values, within each `--comparegrp1` value**. So `--comparegrp2` is the axis being compared, and `--comparegrp1` is the series the bars are coloured by. One figure is written per pair of `--comparegrp2` values, per count grouping (`amino` and `iso`).

> Fold changes come from the same PyDESeq2 negative-binomial fit as the heatmap and volcano, taken over the `amino` or `iso` axis, with BH-adjusted p-values. Each `--comparegrp1` value is fitted separately on its own subset rather than as an interaction term. Because the fit re-derives size factors from the raw counts, a dataset in which every feature moves together is read as a library-size difference and reported as no change.
>
> The two must name **different** columns: with one column for both, there is no comparison to make, and the run now aborts saying so rather than failing inside pandas with `Grouper for 'group' not 1-dimensional`. Because both flags default to `group`, reaching this plot always means setting at least one of them. A comparison also needs a `--comparegrp2` value shared across every value of `--comparegrp1`; where none exists — always the case when `--comparegrp1` is a per-observation-unique column such as `sample` — that count grouping is skipped with a warning rather than failing. A count grouping is likewise skipped, naming both columns, when a `--comparegrp1` value has only one sample per `--comparegrp2` value: an unreplicated design leaves no residual variation to estimate dispersion from, so no fold change can be fitted.

**Venn Options:**

Only used by `-g venn`. Set-overlap diagrams over feature populations: which tRNAs are detected as one species, the other, or both. Configured entirely through `--config`'s `multivariate` block, not through flags.

- **`multivariate.grouping`**: `obs` column the analysis is taken over. Default: `group`.
- **`multivariate.presence_cutoff`**: mean normalized count, per group, for a feature to count as present in a circle. Default: `20`.
- **`multivariate.venn`**: optional list of complex diagrams, each with its circles enumerated.

Two diagrams — fragment vs full-length, and 5' vs 3' — are drawn automatically whenever the data supports them, with nothing to declare. See [Multivariate Analyses](advanced_usage.md#multivariate-analyses-multivariate) for the block's full shape, how circles are declared, and how diagram size decides the layout.

**Correlation Options:**

- **`--corrmethod`**: Correlation method. Default: `pearson`.
- **`--corrgroup`**: Grouping variable to generate correlation matrices for. Default: `sample`.
- **`--corrmask`**: Hide one half of each matrix: `none` (default), `upper` or `lower`. A correlation matrix is symmetric, so half of it restates the other half. The **diagonal is kept** either way — it is R² = 1 by construction and carries no information, but it anchors each row to its axis label, and dropping it turns the triangle into a stair-step. Applies to every matrix the module writes (per-readtype, non-tRNA, and combined tRNA + non-tRNA). A masked run appends the mask name to its filenames, so a masked and an unmasked run of the same data sit side by side rather than one overwriting the other; `none` keeps the original filenames unchanged.

**Coverage Plot Options:**

- **`--covgrp`**: Grouping variable. Default: `group`. Also the columns of the read-specificity grid (see below), which refuses to draw more than 24 of them — that check runs before any plotting, so a too-wide column costs nothing.
- **`--covtype`**: Which coverage category to plot. Reads are binned by how specifically they could be assigned, into four mutually exclusive categories that sum to total coverage:

  | Alias                   | `adata.var` value        | Meaning                                 |
  | ----------------------- | ------------------------ | --------------------------------------- |
  | `unique` / `transcript` | `uniquecoverage`         | Assigned to exactly one tRNA transcript |
  | `isodecoder`            | `multitrnacoverage`      | One anticodon, several transcripts      |
  | `isotype`               | `multianticodoncoverage` | One amino acid, several anticodons      |
  | `notamino`              | `multiaminocoverage`     | Several amino acids                     |
  | `total`                 | `coverage`               | The sum of all four                     |

  Any other `adata.var` coverage value (`readstarts`, `readends`, `mismatchedbases`, `deletions`, …) is also accepted. Default: `unique`, or `total` under `--allreads`; an explicit value is always honored. Each category is written to its own subfolder named for its alias, so separate `--covtype` runs into the same `-o` never overwrite each other.

- **`--covmethod`**: Combination method (`mean`).
- **`--covobs`**: The basis each individual coverage plot is drawn for — i.e. what one plot represents. Default: `trna` (one plot per tRNA transcript).
- **`--covgap`**: Include alignment gap positions in coverage plots. Default: `False`, since gap columns carry no reads and stretch the x-axis.
- **`--combinedpdfonly`**: Skip individual tRNA PDFs. Default: `False`. Covers the read-specificity individuals too, which are one file per tRNA **per group** and so the larger of the two sets.

> [!NOTE]
> The four categories stacked together are written under `coverage/specificity/`, which shows what fraction of each position's signal is transcript-specific versus only isodecoder-, isotype- or amino-level assignable — reading that from the per-category folders would mean comparing four PDFs by eye. Two views are produced:
>
> - `combined_<covobs>_specificity_by_<covgrp>_<covmethod>.pdf` — a **grid**, one row per `--covobs` value and one column per `--covgrp` value, with the panels in a row sharing a y-axis so the groups can be read against each other. Paginated at 4 rows × 8 columns; a row's scale covers every group, including those on a later column-page, so the pages stay comparable.
> - `specificity/<covobs>/` — the same plot per group at full size, one file per `--covobs` value per group, with quiet ones sorted into `low_coverage/` by the same ceiling the individual coverage plots use.
>
> Both are the same under `--allreads`, which selects a category rather than changing the partition. With a single group the grid is one column wide, which is what this plot looked like before it gained a group axis.

**Heatmap Options:**

- **`--heatgrp`**: Grouping variable. Default: `group`.
- **`--diffrts`**: Read types for differential analysis (shared with volcano). Bare readtypes only (`total`, `wholecounts`, `fiveprime`, `threeprime`, `other`) — the read basis comes from `--allreads`, so a value carrying a `_unique` suffix is rejected. Default: all five.
- **`--heatcutoff`**: Read count cutoff. Default: `80`.
- **`--heatbound`**: How many features to show from each end of the ranking — the heatmap is bounded to the top and bottom N counts rather than rendering every feature. Default: `25`.

> [!NOTE]
> On the **combined** heatmap (`<grp>_combine_..._heatmap.pdf`), each row label carries a marker showing which read type it came from — ● 5′, ■ 3′, ▲ other, ◆ whole, ✕ pre-tRNA — with a legend below the panels. That file stacks every `--diffrts` read type except `total`, so one tRNA appears once per read type under the same name; sorting is by effect size, so the copies scatter. The per-read-type heatmaps carry no markers, since every row there is the same read type.

- **`--heatorient`**: Heatmap layout: `vertical` (default) or `horizontal`. The default puts features on rows and comparisons on columns, with the log2FC and p-value panels side by side. `horizontal` transposes that — features to columns, comparisons to rows — and stacks the two panels, which fits a landscape page once the feature count grows. A horizontal run appends the orientation to its **PDF** filenames so both layouts can sit side by side; the CSV written alongside is the _data_, so it is neither renamed nor transposed and two runs of one analysis still produce exports that diff against each other.
- **`--heatsubplots`**: Also save each individual comparison's heatmap as its own PDF, in an `individual/` subfolder next to the combined multi-page PDFs (which are unaffected). Default: `False`.

**Volcano Options:**

- **`--volgrp`**: Grouping variable used both to define the pairwise group comparisons and to look up per-group colors in `--style`'s `colors` block (the same `<obs_column>: {<value>: <color>}` shape used by PCA, keyed on `--volgrp`'s value, e.g. `"group"`). Default: `group`.
- **`--diffrts`**: Read types to generate per-readtype tRNA volcano plots for, shared with heatmap. See the heatmap entry above for accepted values.
- **`--volcutoff`**: Read count cutoff. Default: `80`.
- **`--vollabels`**: Number of top significant markers to label on each plot, ranked by `|log2FC| * -log10(p-value)`. Default: `100` (labeling every significant marker has unbounded cost on large datasets); pass `0` to disable labels entirely, or any other value to label exactly that many.
- **`--volxlim`**: Force the x-axis half-width to this log2 fold change. By default the axis is capped at the 95th percentile of `|log2FC|` whenever the largest value exceeds 1.5x that percentile, and points beyond the cap are drawn as triangles at the boundary — so one extreme feature cannot compress every other one. Nothing is ever dropped, only pinned to the edge.

**Differential Expression Options:**

- **`--shrink`**: How the log2 fold changes are shrunk: `apeGLM` (default) or `none`. p-values are unaffected either way. Shrinking costs one DESeq2 fit per distinct baseline group instead of one overall, so `--shrink none` is faster for iteration. The method is part of the cached `uns['log2FC']` key, so changing it recomputes rather than serving a frame produced by the other setting. See [Coming from tRAX](improvements_from_trax.md) for how this compares to tRAX.

> [!NOTE]
> Two extra subplots are generated automatically (in combined plots) whenever `adata.uns['nontRNA_counts']` is non-empty (i.e., `--gtf` was used at `analyze build`): a non-tRNA-only plot and a combined tRNA + non-tRNA plot. No additional flags are needed, and both are skipped with a log message if non-tRNA counts are unavailable. These use a different DESeq2 normalization than the per-readtype plots — see [Data Structure: Graphing Notes](data_structure.md#6-graphing-notes).

**PCA Options:**

- **`--pcamarkers`**: Marker style grouping. Default: `sample`.
- **`--pcacolors`**: Color grouping. Default: `group`.
- **`--pcareadtypes`**: Read types to generate per-readtype tRNA PCA plots for (`tRNA_<pcamarkers>_by_<pcacolors>_<readtype>_*`). Bare readtypes only — the read basis comes from `--allreads`. Default: `total`. Pass `all` for every readtype that exists in both bases. The combined overview always adds both bases of `total` regardless.

> [!NOTE]
> Two extra plots are generated automatically whenever `adata.uns['nontRNA_counts']` is non-empty (i.e., `--gtf` was used at `analyze build`): a non-tRNA-only plot and a combined tRNA + non-tRNA plot. No additional flags are needed, and both are skipped with a log message if non-tRNA counts are unavailable. These use a different DESeq2 normalization than the per-readtype plots — see [Data Structure: Graphing Notes](data_structure.md#6-graphing-notes).

**Radar Options:**

- **`--radargrp`**: Grouping variable. Default: `group`.
- **`--radarscaled`**: Scale axes to 100%. Default: `False`.
- **`--radarmethod`**: Aggregation method(s) used to combine samples within a group. Repeatable; pass `all` for every available method. Default: `mean`.

**SeqLogo Options:**

- **`--logogrp`**: Grouping variable. Default: `amino`.
- **`--logopseudocount`**: Pseudocount added. Default: `20`.
- **`--logosize`**: Sequence size preset. Default: `noloop`.
- **`--ccatail`**: Keep CCA tail. Default: `True`.
- **`--pseudogenes`**: Keep pseudo-tRNAs. Default: `True`.
- **`--logomanualgrp`**: An explicit list of tRNAs to build one seqlogo from, instead of grouping by the `--logogrp` `obs` column. Repeatable. Default: unset (group by `--logogrp`).
- **`--logomanualname`**: Output filename for the `--logomanualgrp` logo. Default: unset, in which case the file is timestamped. Ignored unless `--logomanualgrp` is given.
- **`--logornamode`**: Render the logo as RNA (U) rather than DNA (T). Default: `False`.

**Mismatch Options:**

- **`--mismatchpseudocount`**: Pseudocount added to coverage when computing per-position misincorporation rates, damping positions with almost no reads behind them. Default: `10`. Graph-time only — the build-time `results/mismatch/` outputs keep tRAX's own constants so they stay directly comparable to a tRAX run.

## Tools Commands

The `tools` module contains utility functions for specific analysis tasks or testing.

### `info`

Reports what an AnnData object contains: every `obs` and `var` column with its dtype, its
number of unique values and those values, the `uns` keys, the `layers`/`obsm`/`obsp` slots, and
the `--variant` strings the object can actually resolve. This is the command to run when you
need to know what to type after a grouping option such as `--covgrp`, `--volgrp` or `--covtype`
-- those vocabularies live inside the `.h5ad` and nowhere else.

**Usage:**

```bash
trnagraph tools info -i <input.h5ad> [options]
```

**Options:**

- **`-i`, `--input`** (Required): Input `.h5ad` file.
- **`--column`**: Print one `obs` or `var` column's unique values in full, instead of the whole
  object. In the default report each column's values are capped at 20 and each individual value
  is elided past 32 characters, because a column such as `trna` (hundreds of values) or
  `refseq` (75-90nt sequences) would otherwise bury every column you would actually type.
- **`--json`**: Emit the report as JSON rather than text, for scripting.

> [!NOTE]
> Columns that hold measurements rather than labels -- a numeric column with more than 20
> distinct values, such as `nreads_total_norm` -- are reported as a range instead of an
> enumeration. A low-cardinality numeric column (an ordered timepoint or dose) keeps its values,
> since that is a label and is exactly the sort of column a grouping option gets pointed at.

### `log2fc`

Computes log2 fold change for specified groups and read types.

**Usage:**

```bash
trnagraph tools log2fc -i <input.h5ad> [options]
```

- **`-i`, `--anndata`** (Required): Input file.
- **`-g`, `--group`**: Grouping variable from obs. Default: `group`. A column that does not exist in the object aborts the run (see the note under [`graph`](#graph)).
- **`-r`, `--readtypes`**: List of read types to analyze.
- **`-x`, `--cutoff`**: Read count cutoff(s). Default: `80`.
- **`--config`**: JSON configuration file whose `flags.log2fc` block pins this command's options; its `name` also namespaces the computed log2FC cache. The `-c` short form was removed, since `-c` means `--readlengthsplit` on `analyze build`. See [Run Configuration](advanced_usage.md#run-configuration---config).
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

### `template`

Writes a blank, fully-enumerated JSON config template into the current directory. Every key the file format accepts is present and set to `null`, so the template is a valid no-op until you fill something in -- it documents the format, which JSON itself cannot do. (The config template's `name` is the one exception: it is a required field that becomes an output subdirectory, so it ships as `"default"`.)

**Usage:**

```bash
trnagraph tools template [options]
```

**Options:**

- **`--style`**: Write the `--style` template (colors, gradients, categorical palette and presentation settings).
- **`--config`**: Write the `--config` template (data filters and pinned graph options). With no selector, every available template is written.
- **`-o, --output`**: Directory to write the template(s) into. Default: `.`.
- **`--overwrite`**: Replace an existing template file instead of refusing. Default: `False` -- the expected workflow is emit, edit, then re-emit after an upgrade, so an edited file is never replaced silently.

### `test`

Runs the internal automated test suite.

**Usage:**

```bash
trnagraph tools test [options]
```

**General Options:**

- **`--all`**: Run all tests including split analysis, forcing a clean workspace and full redownload. Refuses to run if the workspace contains files the suite did not create, naming them -- `--all` recursively deletes the workspace's contents, so it only runs in a directory the suite owns. There is no override; move the files or point `-d` elsewhere.
- **`--skip-download`**: Skip the metadata/fastq/tRNA/genome download steps and run everything else. Downloads are already skipped by default when the target files are already present; this forces the skip regardless.
- **`--cleanrun`**: Clean up test files after completion.
- **`-d`, `--directory`**: Directory to run tests in.
- **`--log`**: Disables the live progress panel, so the run prints plain sequential lines instead. Use it when the panel's redrawing is unhelpful — under `nohup`, in CI, or when piping to a file. It takes no value: the suite always writes its own fixed `toolsTestSuite.log`, and each `trnagraph` invocation it makes writes a timestamped log under `.log/` regardless of this flag.

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
| `--merge`       | Exercise `tools merge` on the objects the suite built          |
| `--hubonly`     | Generate UCSC track hubs only                                  |

---

## Update

### `update`

Updates this git checkout to the latest source and re-syncs the local environment. Refuses to run if the checkout has uncommitted local changes to tracked files, so nothing is at risk of being lost or silently merged over.

Requires an editable (`pip install -e .`) install: the command updates a working tree with git and re-syncs conda from `requirements.yaml`, and neither exists in a built distribution. From a non-editable install it stops immediately with an explanatory error, before running any git command -- it will not act on a directory that merely happens to sit inside some other repository.

**Usage:**

```bash
trnagraph update [options]
```

**Options:**

- **`--branch <name>`**: Update to this branch instead of the one currently checked out (e.g. `dev`). Once given, this branch becomes what a later plain `trnagraph update` (no `--branch`) continues tracking -- git's own branch state is the memory, not a separately stored preference. Mutually exclusive with `--tag`. Refused if the branch's own reported version is older than what's currently installed (see below). If run with no `--branch` from a detached HEAD (e.g. right after a `--tag` checkout), there's no "current branch" to default to, so `--branch` must be given explicitly.
- **`--tag <version>`**: Check out this release tag instead of a branch (e.g. `v1.9.0` or `1.9.0`). This leaves the checkout in git's standard "detached HEAD" state for a tag checkout -- a printed message explains what that means and how to get back to a branch afterward. Mutually exclusive with `--branch`. Refused if the requested version is older than v1.9.0, the version `update` was itself introduced in.

`trnagraph --version` (and the header of every command's log) also shows which channel this checkout is running: `stable` (the `main` branch), `prerelease` (the `dev` branch), or `nightly @ <commit hash>` (any other branch, or a detached HEAD not sitting exactly on a release tag -- which instead just shows that tag).
