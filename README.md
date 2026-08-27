# tRNAgraph

![tRNAgraph Logo](docs/images/logo.png)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/590619343.svg)](https://doi.org/10.5281/zenodo.14669314)
[![Build Status](https://github.com/alba1735/tRNAgraph/actions/workflows/release.yml/badge.svg)](https://github.com/alba1735/tRNAgraph/actions)
![Python](https://img.shields.io/badge/python-3.12-blue.svg)

tRNAgraph is a comprehensive toolkit for analyzing tRNA-seq data. Built upon the foundation of [tRAX](https://github.com/UCSC-LoweLab/tRAX), it generates [AnnData](https://anndata.readthedocs.io/en/latest/index.html) objects, allowing for high-dimensional visualization, clustering, and differential expression analysis.

## Overview

tRNAgraph bridges the gap between raw alignment data and biological insight. It allows users to:

1. **Preprocess** raw FASTQ files (Trim, Map, Index).
2. **Analyze** data by building a structured database (AnnData) and performing clustering.
3. **Graph** results using a suite of built-in visualization tools (Heatmaps, PCA, Coverage Plots, etc.).

```mermaid
flowchart LR
  %% Input Nodes
  I1[/"manifest.tsv"/]
  I2[/"FASTQ Files"/]
  I3[/"Genomic References"/]
  I4[/"metadata.tsv"/]

  %% Output Nodes
  O1[Visualizations]

  subgraph Step1 [1. Preprocess]
    P1[Trim Reads]
    P2[Make Database]
    P3[Map Reads]
  end

  I3 --> P2
  I2 & I1 --> P1

  D1[("Bowtie2 Index")]
  D2[("Trimmed Reads")]
  P2 --> D1
  P1 --> D2
  D1 & D2 --> P3

  D3[("BAM Directory")]
  P3 --> D3

  subgraph Step2 [2. Analyze]
    B1[Build AnnData]
    C1[Cluster Data]
  end

  D3 & I4 --> B1

  D4[("tRNAgraph.h5ad")]
  D5[("Results Directory")]
  B1 -->|creates| D4
  B1 -->|creates| D5
  C1 -.->|updates in place| D4

  subgraph Step3 [3. Graph]
    V1[Generate Graphs]
  end

  D4 --> V1 --> O1

  %% Styling
  classDef input fill:#000000,stroke:#01579b,stroke-width:2px;
  classDef process fill:#000000,stroke:#2e7d32,stroke-width:2px;
  classDef storage fill:#000000,stroke:#ef6c00,stroke-width:2px;
  classDef output fill:#000000,stroke:#7b1fa2,stroke-width:2px;

  class I1,I2,I3,I4 input;
  class P1,P2,P3,B1,C1,V1 process;
  class D1,D2,D3,D4,D5 storage;
  class O1 output;
```

## Documentation

- **[Installation & Quick Start](#installation)**: Get up and running.
- **[Tutorial](docs/tutorial.md)**: Three worked examples, from a reproducible yeast run to a full human OTTR-seq analysis.
- **[CLI Reference](docs/cli_reference.md)**: Detailed documentation for all commands and flags.
- **[Data Structure](docs/data_structure.md)**: Details on the AnnData object, observations, and variables.
- **[Advanced Usage](docs/advanced_usage.md)**: Python API, configuration files, and downstream analysis.
- **[Test Suite](docs/testSuite.md)**: How to run the automated validation pipeline.
- **[Roadmap](docs/roadmap.md)**: Future features and planned improvements.

## Installation

Dependencies can be installed using conda/mamba, and the package itself is installed via pip:

```bash
# 1. Create the environment with non-Python dependencies (bowtie2, etc.)
conda env create -f requirements.yaml
conda activate trnagraph

# 2. Install tRNAgraph in editable mode
pip install -e .
```

> [!NOTE]
> The editable (`-e`) install is the supported mode. A plain `pip install .` gives you a working
> CLI, but `trnagraph update` and the `requirements.yaml` environment check both need the source
> checkout and will report that it is unavailable rather than acting on a guessed location.

## Quick Start

### 1. Prepare Inputs

You need two tab-delimited files to begin:

**manifest.tsv** (For trimming and merging reads):
_Format: OutputPrefix `<tab>` R1_Path `<tab>` R2_Path (optional)_

> [!NOTE]
> If single-end reads are used, only provide the R1 column.

> [!NOTE]
> OutputPrefix doubles as the sample name **and** where the trimmed output is written: a bare name (as below) writes to `processed/trimmed/<name>_trimmed.fastq.gz`; a name containing a directory (e.g. `some/dir/SampleA`) writes there instead.

```text
SampleA   <path_to_fastq>/VC_24h_1_R1.fastq.gz   <path_to_fastq>/VC_24h_1_R2.fastq.gz
SampleB   <path_to_fastq>/VC_24h_2_R1.fastq.gz   <path_to_fastq>/VC_24h_2_R2.fastq.gz
```

> [!TIP]
> A `trim_metadata.tsv` template is automatically generated in the output directory after trimming. While this simplifies metadata creation, it defaults to setting `group` names identical to `sample` names. **You must update this file** with your actual experimental groups to ensure accurate normalization and downstream analysis.

**metadata.tsv** (For mapping reads and building the database):
_Format: Must contain `fastq`, `sample` and `group` columns. Add other metadata columns as needed._

> [!TIP]
> More metadata columns allow for richer analysis and visualization.

```text
fastq	sample	group	treatment
processed/trimmed/SampleA_merged.fastq.gz	SampleA	Control	None
processed/trimmed/SampleB_merged.fastq.gz	SampleB	Treated	DrugX
```

### 2. Preprocess Data

Run the integrated wrapper for trimming (fastp), tRNA database generation (bowtie2), mapping (bowtie2):

```bash
# Create database
trnagraph preprocess makedb -g genome.fa -t trnascan.out -r gtrna.fa -m namemap.txt

# Trim reads (writes to processed/trimmed/ by default -- see the manifest note above)
trnagraph preprocess trim -i manifest.tsv

# Map reads (-i takes metadata.tsv, not the trim manifest; -o names this experiment)
trnagraph preprocess map -i metadata.tsv -d database -o experiment1
```

### 3. Build Database

Convert coverage files into a tRNAgraph AnnData object, attaching your metadata:

```bash
trnagraph analyze build -i metadata.tsv -o build_output -d db_name
```

### 4. Generate Graphs

Create a standard suite of visualizations. Convention is to write graphs into a `graphs/` subfolder alongside `results/`, inside the same output directory `build` used, so everything for one experiment stays together:

```bash
trnagraph graph -i build_output/build_output.h5ad -o build_output/graphs -g all
```

## Output Structure

A standard run of the full pipeline will generate the following directory structure with default settings:

```text
project_root/
├── db/                              # Generated by 'preprocess makedb' (Default: db/)
│   ├── database.1.bt2
│   ├── database.fa
│   └── ...
├── processed/                       # Generated by 'preprocess'
│   ├── trimmed/                     # Generated by 'trim' (Default: processed/trimmed)
│   │   ├── SampleA_trimmed.fastq.gz
│   │   ├── ...
│   │   ├── trim_feature_types.pdf   # Example QC plot
│   │   └── trim_stats.csv
│   │
│   └── bam/                         # Generated by 'map' (Default: processed/bam)
│       ├── SampleA.bam
│       └── ...
└── <output_dir>/                    # Passed to 'build -o' -- name it anything you like
    ├── <output_dir>.h5ad            # The main AnnData database object
    ├── results/                     # Fixed subfolder name -- per-run text output files
    │   ├── <output_dir>-runinfo.txt        # Build provenance for this run
    │   ├── <output_dir>-readcounts.txt
    │   ├── ...                             # normalizedreadcounts/typecounts/aminocounts/etc.
    │   ├── trna/                           # tRNA-only-matrix DESeq2 output
    │   ├── allfeature/                     # All-feature-controlled DESeq2 output (tRAX parity)
    │   └── unique/                         # Uniquely-mapped-read count/DE output
    └── graphs/                      # Wherever you pass to 'graph -o' -- convention is
        ├── bar/                     # <output_dir>/graphs, alongside results/, as shown here
        ├── cluster/                 # UMAP/HDBSCAN plots
        ├── correlation/             # Correlation plots
        ├── count/                   # Count summary plots
        ├── coverage/                # Coverage profiles per tRNA
        ├── heatmap/                 # Differential expression heatmaps
        ├── logo/                    # Sequence logo plots
        ├── pca/                     # PCA plots
        ├── radar/                   # Radar plots
        └── volcano/                 # Volcano plots for DE
```

> [!NOTE]
> Passing `--readlengthsplit <N>` to `analyze build` adds two more variants -- `u<N>` (fragments under the cutoff) and `o<N>` (over it) -- computed alongside the default/full-length one and merged into the same `.h5ad`. When (and only when) a split is requested, `results/`'s per-run text output nests into `results/complete/`, `results/u<N>/`, and `results/o<N>/` subfolders instead of writing straight into `results/`, so the three variants don't collide; a plain build with no split keeps writing directly to `results/` as shown above. `graph -o` isn't affected by this directly (its output directory is always whatever you pass to `-o`), but the convention this project's own demo pipeline (`tools test`) follows -- and the one shown above -- is to mirror the same nesting under `graphs/` (`graphs/complete/`, `graphs/u<N>/`, `graphs/o<N>/`) when graphing more than one variant into the same output directory.

See the [CLI Reference](docs/cli_reference.md) and [Data Structure](docs/data_structure.md) docs for the full list of per-run output files and what each contains.

## License

tRNAgraph is licensed under the [GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) license.
