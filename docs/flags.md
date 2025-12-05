# tRNAgraph Flag Reference

This document provides a detailed reference for all command-line flags available in tRNAgraph.

## Global Options

These options can be used with any command.

- **`--log <file>`**

  - **Description**: Redirects all output logging to the specified file instead of standard output. Useful for keeping records of long-running jobs.
  - **Default**: `None` (Output to stdout)

- **`-q`, `--quiet`**

  - **Description**: Suppresses all standard output (stdout). Errors and warnings will still be printed to stderr.
  - **Default**: `False`

- **`-v`, `--verbose`**
  - **Description**: Enables verbose output for commands that support it (e.g., `trim`). This will print detailed command execution logs.
  - **Default**: `False`

---

## Preprocess Commands

The `preprocess` module handles the initial processing of raw sequencing data, including database creation, trimming, and mapping.

### `makedb`

Builds a Bowtie2 index from gtRNAdb/tRNAScan-SE output and a reference genome. This index is required for the `map` step.

**Usage:**

```bash
python trnagraph.py preprocess makedb -g <genome> -t <trnaout> -r <trnafa> -m <namemap> [options]
```

#### Flags

- **`-g`, `--genome`** (Required)

  - **Description**: Path to the reference genome FASTA file. This is used to extract the exact sequences of the tRNAs including flanking regions.

- **`-t`, `--trnaout`** (Required)

  - **Description**: Path to the tRNAScan-SE output file. This file contains the genomic coordinates of predicted tRNAs.

- **`-r`, `--trnafa`** (Required)

  - **Description**: Path to the tRNA reference FASTA file (usually from gtRNAdb).

- **`-m`, `--namemap`** (Required)

  - **Description**: Path to the tRNA name mapping file. This maps the tRNAScan-SE IDs to standard tRNA names.

- **`-s`, `--orgmode`**

  - **Description**: Organism mode used for tRNAScan-SE. This ensures the correct scoring model and structure are used.
  - **Options**: `euk` (Eukaryotic), `bact` (Bacterial), `arch` (Archaeal).
  - **Default**: `euk`

- **`-o`, `--output`**

  - **Description**: Output prefix for the generated Bowtie2 index files.
  - **Default**: `db` (Creates files like `db.1.bt2`, `db.rev.1.bt2`, etc.)

- **`-n`, `--threads`**

  - **Description**: Number of threads to use for building the index.
  - **Default**: `0` (Uses all available CPU cores).
  - **Note:** The number of threads used for mapping is highly system dependent. Bowtie2 can be memory intensive, and using too many threads can cause the system to run out of memory or lose performance because of overhead. It is recommended to use a number of threads that is appropriate for your system's available memory and CPU cores although between 8-10 has been commonly used on high performance machines as an optimal range and starting point.

- **`--addtrna`**

  - **Description**: Path to a FASTA file containing additional tRNA sequences to include in the database (e.g., mitochondrial tRNAs if not in the main set).
  - **Default**: `None`

- **`--addseqs`**

  - **Description**: Path to a FASTA file containing additional non-tRNA sequences to include (e.g., spike-ins, rRNA).
  - **Default**: `None`

- **`--forcecca`**
  - **Description**: Forces the addition of a CCA tail to the 3' end of all tRNA sequences in the database.
  - **Default**: `False`

---

### `trim`

Trims adapters, merges paired-end reads, and extracts UMIs using `fastp`. This prepares the FASTQ files for mapping.

**Usage:**

```bash
python trnagraph.py preprocess trim -n <name> -i <manifest> [options]
```

#### Flags

- **`-n`, `--name`** (Required)

  - **Description**: Name of the run. This is used to prefix output filenames and directories.

- **`-i`, `--manifest`** (Required)

  - **Description**: A tab-delimited file describing the samples.
  - **Format**: `SampleName <tab> R1_Path [<tab> R2_Path]`
  - **Note**: If R2_Path is omitted, the sample is treated as single-end.

- **`-a1`, `--adapter1`**

  - **Description**: Adapter sequence for Read 1. If not specified, `fastp` will attempt to auto-detect it.
  - **Default**: `None` (Auto-detect)

- **`-a2`, `--adapter2`**

  - **Description**: Adapter sequence for Read 2 (for paired-end data).
  - **Default**: `None` (Auto-detect)

- **`-l`, `--length`**

  - **Description**: Minimum sequence length allowed after trimming. Reads shorter than this will be discarded.
  - **Default**: `15`

- **`-u`, `--umilength`**

  - **Description**: Length of the Unique Molecular Identifier (UMI) in base pairs. Set to `0` to disable UMI extraction.
  - **Default**: `0`

- **`--umi3`**

  - **Description**: Specifies that the UMI is located at the 3' end of the read. If not set, it is assumed to be at the 5' end.
  - **Default**: `False` (5' end)

- **`-t`, `--threads`**
  - **Description**: Total number of threads to use for `fastp`.
  - **Default**: `0` (Uses all available CPU cores).

---

### `map`

Maps reads to the tRNA database using Bowtie2 and generates coverage files.

**Usage:**

```bash
python trnagraph.py preprocess map -n <name> -d <database> -s <samples> [options]
```

#### Flags

- **`-n`, `--name`** (Required)

  - **Description**: Experiment name. This is used to name the output directory and report files.

- **`-d`, `--database`** (Required)

  - **Description**: Path/Name of the Bowtie2 index (tRNA database) created by `makedb`.

- **`-s`, `--samples`** (Required)

  - **Description**: Path to the sample file (usually the manifest file used in `trim` or a similar format).

- **`-t`, `--threads`**

  - **Description**: Number of threads to use for Bowtie2 mapping.
  - **Default**: `8`
  - **Note**: Increasing threads significantly speeds up mapping.

- **`--bamdir`**

  - **Description**: Directory where BAM files will be stored.
  - **Default**: `bam/<name>`

- **`--lazy`**

  - **Description**: If set, skips the mapping step if the expected BAM files already exist. Useful for re-running downstream steps without re-mapping.
  - **Default**: `False`

- **`--minnontrnasize`**

  - **Description**: Minimum read length required for a read to be assigned to a non-tRNA feature.
  - **Default**: `20`

- **`--local`**

  - **Description**: Configures Bowtie2 to use local alignment mode (`--local`) instead of end-to-end.
  - **Default**: `False`

- **`--skipcheck`**
  - **Description**: Skips the validation check that ensures FASTQ read names match BAM headers. Use with caution.
  - **Default**: `False`

---

## Core Commands

### `build`

Builds an AnnData object (`.h5ad`) from the output of `preprocess map`. This object is the core data structure for all subsequent analysis in tRNAgraph.

**Usage:**

```bash
python trnagraph.py build -i <metadata> -o <output_directory> -d <database> [options]
```

#### Flags

- **`-i`, `--input`, `--metadata`** (Required)

  - **Description**: Path to the metadata file describing the samples (groups, conditions, etc.).
  - **Note**: Can also be the sample file if no extra metadata is needed, but a header row is recommended.

- **`-o`, `--output`** (Required)

  - **Description**: Output directory. The AnnData file will be created inside this directory with the name `<basename>.h5ad` alongside a results/ directory with additional output files.
  - **Default**: `h5ad`

- **`-d`, `--database`** (Required)

  - **Description**: Path/Name of the Bowtie2 index (tRNA database) created by `makedb`.

- **`-t`, `--threads`**

  - **Description**: Number of threads to use for processing.
  - **Default**: `8`

- **`--gtf`**

  - **Description**: Path to an Ensembl GTF file. Used to annotate non-tRNA features.
  - **Default**: `None`

- **`--pairs`**

  - **Description**: File listing sample pairs for direct comparison.
  - **Default**: `None`

- **`--bed`**

  - **Description**: List of additional BED files to define custom features.
  - **Default**: `None`

- **`--nofrag`**

  - **Description**: Omits fragment determination logic. Useful for specific protocols like TGIRT-seq where fragmentation patterns differ.
  - **Default**: `False`

- **`--nosizefactors`**

  - **Description**: Disables the use of DESeq2 size factors for normalization.
  - **Default**: `False`

- **`--maxmismatches`**

  - **Description**: Maximum number of allowed mismatches for a read to be counted.
  - **Default**: `None` (No limit)

- **`--mincoverage`**

  - **Description**: Minimum read count required for a transcript to be included in coverage plots.
  - **Default**: `None`

- **`--minnontrnasize`**

  - **Description**: Minimum read length for non-tRNA features.
  - **Default**: `20`

- **`--hub`**

  - **Description**: Generates a UCSC Genome Browser track hub for the data.
  - **Default**: `False`

- **`--hubonly`**

  - **Description**: Only generates the track hub and skips AnnData generation.
  - **Default**: `False`

- **`--dumpother`**

  - **Description**: Includes 'other' features (non-tRNA) when counting gene types.
  - **Default**: `False`

- **`--bamdir`**

  - **Description**: Custom directory to look for BAM files if they are not in the default location.
  - **Default**: `None`

- **`--uniqueonly`**

  - **Description**: Restricts analysis to uniquely mapped reads only.
  - **Default**: `False`

- **`--dispfittype`**
  - **Description**: Dispersion fit type for DESeq2.
  - **Options**: `parametric` (default, standard for DESeq2), `mean` (robust for small sample sizes).
  - **Default**: `parametric`
  - **Note**: If the number of samples is small, this will automatically resolve to `mean` to ensure stability.

---

### `cluster`

Performs clustering (UMAP, HDBSCAN) on an existing AnnData object.

> [!NOTE]
> The database object can be clustered using the `cluster` function. The following code will cluster the database object using [UMAP](https://umap-learn.readthedocs.io/en/latest/index.html) and cluster using [HDBSCAN](https://hdbscan.readthedocs.io/en/latest/index.html). The default parameters used in tRNAgraph work well on ARMseq, DM-tRNAseq, and OTTRseq data; however, each dataset is different and may require fine-tuning of the parameters to yield the best results. The following code will cluster the database object and save the result to a new database object:

**Usage:**

```bash
python trnagraph.py cluster -i <input.h5ad> -o <output.h5ad> [options]
```

#### Flags

- **`-i`, `--anndata`** (Required)

  - **Description**: Path to the input AnnData (`.h5ad`) file.

- **`-o`, `--output`**

  - **Description**: Path for the output AnnData file containing cluster information.
  - **Default**: `trnagraph.cluster.h5ad`

- **`-w`, `--overwrite`**

  - **Description**: Overwrites existing cluster information in the file if present.
  - **Default**: `False`

- **`-r`, `--randomstate`**

  - **Description**: Seed for random number generation to ensure reproducible UMAP results.
  - **Default**: `None` (Random)

- **`-t`, `--readcutoff`**

  - **Description**: Minimum read count for a tRNA to be included in clustering.
  - **Default**: `20`

- **`--coveragetype`**

  - **Description**: List of coverage features to use for clustering.
  - **Default**: `['uniquecoverage', 'readstarts', 'readends', 'mismatchedbases', 'deletions']`

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

> [!NOTE]
> Clustering is performed across the `uniquecoverage`, `readstarts`, `readends`, `mismatchedbases`, and `deletions` categories of the AnnData object. When performing clustering, verifying that it is reproducible and that the results reflect the data is important. This can be done by running the clustering multiple times and comparing the results. The clustering is also performed on `sample` and `group` observations. In the case of samples, every set of reads for every single tRNA is used for clustering. In the case of groups, the mean of the reads is taken for each tRNA across the read categories and then used for clustering. This is done to reduce the number of samples used for clustering and to reduce the noise in the clustering. The results will be saved in the `obs` attribute of the database object as `sample_cluster\umap1\umap2` and `group_cluster\umap1\umap2`, respectively. Clusters annotated as `-1` are considered noise and are not included in the clustering. Plotting of the clustering is done as well for convenience.

> [!WARNING]
> When working with downstream analysis of the cluster groups, it is important to note that reads that are dropped via the `--readcutoff` flag will not be included in the clustering however, they are still present in the AnnData object. This means that your object can contain NaN values in the clustering columns. Depending on your use case, you may want to filter these out before performing any analysis.

---

### `merge`

Merges two AnnData objects into a single file.

**Usage:**

```bash
python trnagraph.py merge -i1 <file1.h5ad> -i2 <file2.h5ad> -o <merged.h5ad>
```

#### Flags

- **`-i1`, `--anndata1`** (Required): First AnnData file.
- **`-i2`, `--anndata2`** (Required): Second AnnData file.
- **`-o`, `--output`**: Output file path (Default: `trnagraph.merge.h5ad`).
- **`--dropno`**: Drop non-tRNA genes that are not present in both objects (Default: `False`).
- **`--droprna`**: Drop RNA categories that are not present in both objects (Default: `False`).

---

### `graph`

Generates a wide variety of visualizations from the AnnData object.

**Usage:**

```bash
python trnagraph.py graph -i <input.h5ad> -o <output_dir> -g <graph_type> [options]
```

#### Flags

- **`-i`, `--anndata`** (Required)

  - **Description**: Path to the input AnnData file.

- **`-o`, `--output`**

  - **Description**: Directory where figures will be saved.
  - **Default**: `figures`

- **`-g`, `--graphtypes`**

  - **Description**: List of graph types to generate.
  - **Options**: `all`, `bar`, `heatmap`, `pca`, `radar`, `logo`, `volcano`, `coverage`, `count`, `correlation`.
  - **Default**: `all`

- **`-n`, `--threads`**

  - **Description**: Number of threads to use for parallel graph generation.
  - **Default**: `0` (All available)

- **`--config`**

  - **Description**: Path to a JSON configuration file for filtering data (observations/variables).
  - **Default**: `None`

- **`--colormap`**

  - **Description**: Path to a JSON file defining custom colors for groups/features.
  - **Default**: `None`

- **`--regen_uns`**
  - **Description**: Forces regeneration of `uns` (unstructured) data like log2fc even if it exists.
  - **Default**: `False`

#### Graph-Specific Options

**Bar Plots:**

- **`--barcol`**: Grouping for individual bar stacks (Default: `group`).
- **`--bargrp`**: Grouping for bar columns (Default: `amino`).
- **`--barsubgrp`**: Secondary grouping for subplots.
- **`--barsort`**: Column to sort bars by.

**Cluster Plots:**

- **`--clustergrp`**: Grouping variable (Default: `amino`).
- **`--clusteroverview`**: Generate overview plot (Default: `False`).
- **`--clustermask`**: Mask unclustered points (Default: `False`).

**Coverage Plots:**

- **`--covgrp`**: Grouping variable (Default: `group`).
- **`--covtype`**: Coverage type to plot (Default: `uniquecoverage`).
- **`--covmethod`**: Method for combining coverage (Default: `mean`).
- **`--combinedpdfonly`**: Only generate combined PDFs, skip individual ones.

**Heatmap / Volcano:**

- **`--heatgrp` / `--volgrp`**: Grouping variable (Default: `group`).
- **`--heatcutoff` / `--volcutoff`**: Read count cutoff (Default: `80`).
- **`--diffrts`**: Read types to use (Default: `wholecounts_unique`, etc.).

**PCA:**

- **`--pcamarkers`**: Marker style grouping (Default: `sample`).
- **`--pcacolors`**: Color grouping (Default: `group`).

**Radar:**

- **`--radargrp`**: Grouping variable (Default: `group`).
- **`--radarscaled`**: Scale axes to 100% (Default: `False`).

**SeqLogo:**

- **`--logogrp`**: Grouping variable (Default: `amino`).
- **`--logopseudocount`**: Pseudocount added (Default: `20`).
- **`--logosize`**: Sequence size preset (Default: `noloop`).
- **`--ccatail`**: Keep CCA tail (Default: `True`).

---

## Tools Commands

### `log2fc`

Computes log2 fold change for specified groups and read types.

- **`-i`, `--anndata`**: Input file.
- **`-g`, `--group`**: Grouping variable (Default: `group`).
- **`-r`, `--readtypes`**: Read types to analyze.
- **`-x`, `--cutoff`**: Read count cutoff (Default: `80`).

### `csv`

Exports the AnnData object to a set of CSV files.

- **`-i`, `--anndata`**: Input file.
- **`-o`, `--output`**: Output directory (Default: `csv`).

### `test`

Runs the internal test suite.

- **`--all`**: Run all tests.
- **`--cleanrun`**: Remove test files after completion.
- **`-d`, `--directory`**: Directory to run tests in.
