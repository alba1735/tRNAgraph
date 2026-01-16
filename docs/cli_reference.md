# Command Line Interface (CLI) Reference

This document provides a detailed reference for all command-line commands and flags available in tRNAgraph.

## Global Options

These options apply to most commands in the toolkit.

| Flag              | Description                                                    | Default         |
| :---------------- | :------------------------------------------------------------- | :-------------- |
| `--log <file>`    | Redirects all output logging to the specified file.            | `None` (Stdout) |
| `-q`, `--quiet`   | Suppresses standard output (stdout). Errors are still printed. | `False`         |
| `-v`, `--verbose` | Enables detailed execution logs.                               | `False`         |
| `--skip-env-check`| Skips the environment validation checks (dependencies and versions). | `False` |

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

- **`-i`, `--manifest`** (Required): Tab-delimited file (`SampleName <tab> R1 <tab> R2`).
- **`-a1`, `--adapter1`**: Adapter sequence for R1 (Auto-detect if omitted).
- **`-a2`, `--adapter2`**: Adapter sequence for R2 (for paired-end data, Auto-detect if omitted).
- **`-l`, `--length`**: Minimum sequence length allowed after trimming. Reads shorter than this will be discarded. Default: `15`.
- **`-u`, `--umilength`**: Length of the Unique Molecular Identifier (UMI) in base pairs. Set to `0` to disable UMI extraction. Default: `0`.
- **`--umi3`**: Specifies that the UMI is located at the 3' end of the read. If not set, it is assumed to be at the 5' end. Default: `False` (5' end).
- **`-n`, `--threads`**: Number of threads to use for fastp. Default: max_cores.

> [!NOTE]
> If R2_Path is omitted in the manifest file, the sample will be treated as single-end.

### `split`

Splits BAM files based on read length. This is useful for separating tRNAs into different size fractions (fragments vs full-length) (e.g., <60bp and >=60bp) for separate analysis. This function generates new metadata files for the split datasets in addition to the split BAM files.

**Usage:**

```bash
trnagraph preprocess split -i <metadata> [options]
```

**Flags:**

- **`-i`, `--input`** (Required): Tab-delimited metadata file.
- **`-c`, `--cutoff`**: Read length cutoff for splitting. Reads with length < cutoff go to one file, >= cutoff to another. Default: `60`.
- **`--bamdir`**: Directory containing input BAM files. Default: current directory.
- **`-n`, `--threads`**: Number of threads to use. Default: `0` (1 thread per sample).

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
- **`--bamdir`**: Directory where BAM files will be stored. Default: `bam/<name>`
- **`--lazy`**: Skip mapping if BAMs exist. Default: `False`.
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
- **`--nofrag`**: Omits fragment determination logic. Useful for specific protocols like TGIRT-seq where fragmentation patterns differ.
- **`--nosizefactors`**: Disables the use of DESeq2 size factors for normalization.
- **`--maxmismatches`**: Maximum number of allowed mismatches for a read to be counted. Default: `None` (No limit)
- **`--mincoverage`**: Minimum read count required for a transcript to be included in coverage plots. Default: `None`
- **`--minnontrnasize`**: Minimum read length for non-tRNA features. Default: `20`
- **`--hub`**: Generates a UCSC Genome Browser track hub for the data. Default: `False`
- **`--hubonly`**: Only generates the track hub and skips AnnData generation. Default: `False`
- **`--dumpother`**: Includes 'other' features (non-tRNA) when counting gene types. Default: `False`
- **`--bamdir`**: Custom directory to look for BAM files if they are not in the default location. Default: `None`
- **`--uniqueonly`**: Restricts analysis to uniquely mapped reads only. Default: `False`
- **`--dispfittype`**: Dispersion fit type for DESeq2. parametric`(default, standard for DESeq2),`mean`(robust for small sample sizes). Default:`parametric`
- **`-n`, `--threads`**: Number of threads to use for processing. Default: `8`

> [!NOTE]
> PyDeseq2 will automatically switch to `mean` dispersion fitting if the number of samples is too small for `parametric` to be stable.

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

### `merge`

Merges two existing AnnData objects into a single file. This is useful for combining datasets processed separately.

**Usage:**

```bash
trnagraph analyze merge -i1 <file1.h5ad> -i2 <file2.h5ad> [options]
```

**Flags:**

- **`-i1`, `--anndata1`** (Required): Path to the first AnnData file.
- **`-i2`, `--anndata2`** (Required): Path to the second AnnData file.
- **`-o`, `--output`**: Path for the merged output file. Default: `trnagraph.merge.h5ad`.
- **`--dropno`**: Drop non-tRNA genes that are not present in both objects. Default: `False`.
- **`--droprna`**: Drop RNA categories that are not present in both objects. Default: `False`.

---

### `graph`

Generates a wide variety of visualizations from the AnnData file.

**Usage:**

```bash
trnagraph graph -i <input.h5ad> -o <output_dir> [options]
```

**General Flags:**

- **`-i`, `--input`** (Required): Input AnnData file.
- **`-o`, `--output`**: Output directory. Default: `figures`.
- **`-g`, `--graphtypes`**: List of graphs to generate (e.g., `all`, `bar`, `heatmap`, `pca`). Default: `all`.
- **`--config`**: JSON configuration file for filtering.
- **`--colormap`**: JSON file for custom colors.
- **`--regen_uns`**: Force regeneration of calculated stats.

**Bar Plot Options:**

- **`--barcol`**: Column for bar stack components. Default: `group`.
- **`--bargrp`**: Column for bar stacking groups. Default: `amino`.
- **`--barsubgrp`**: Column for subplots.
- **`--barsort`**: Column to sort bars by.

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

**Heatmap / Volcano Options:**

- **`--heatgrp` / `--volgrp`**: Grouping variable. Default: `group`.
- **`--diffrts`**: Read types for differential analysis.
- **`--heatcutoff` / `--volcutoff`**: Read count cutoff. Default: `80`.
- **`--heatsubplots`**: Generate subplots per comparison. Default: `False`.

**PCA Options:**

- **`--pcamarkers`**: Marker style grouping. Default: `sample`.
- **`--pcacolors`**: Color grouping. Default: `group`.

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

### `csv`

Exports the AnnData object to a set of CSV files (obs, var, X).

**Usage:**

```bash
trnagraph tools csv -i <input.h5ad> -o <output_dir>
```

- **`-i`, `--anndata`** (Required): Input file.
- **`-o`, `--output`**: Output directory. Default: `csv`.

### `test`

Runs the internal automated test suite.

**Usage:**

```bash
trnagraph tools test [options]
```

- **`--all`**: Run all tests.
- **`--cleanrun`**: Clean up test files after completion.
- **`-d`, `--directory`**: Directory to run tests in.
- **Step Flags**: `--metadata`, `--fastq`, `--trna`, `--genome`, `--trim`, `--makedb`, `--map`, `--build`, `--cluster`, `--graph`.
