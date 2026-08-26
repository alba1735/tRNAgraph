# tRNAgraph Test Suite

This document outlines the automated test suite for the tRNAgraph project. The suite validates the functionality of the entire pipeline, from raw data acquisition to high-dimensional visualization. It is designed to ensure reproducibility and regression testing by running a complete "end-to-end" workflow.

## Overview

The test suite automates the following workflow:

1.  **Metadata Retrieval**: Downloads and filters experiment metadata from SRA.
2.  **Data Acquisition**: Downloads FASTQ files and subsamples them for speed.
3.  **Reference Preparation**: Fetches _Vibrio cholerae_ genome/annotations and prepares tRNA sequences.
4.  **Preprocessing**: Trims adapters, builds the database, and maps reads.
5.  **Database Construction**: compiles mapping results into an AnnData object.
6.  **Analysis & Visualization**: Performs clustering and generates a full suite of graphs.

## Usage

Invoke the test suite via the CLI:

```bash
python trnagraph.py tools test [options]
```

### Options

- `--all`: Run the complete pipeline including split analysis (cleans existing test data first, forcing a full redownload).
- `--skip-download`: Skip the metadata/FASTQ/tRNA/genome download steps and run everything else. Downloads are already skipped by default when their target files are present (the slowest step to redundantly repeat); this forces the skip regardless.
- `--cleanrun`: Clean up all generated files after the run completes.
- `--directory`: Specify a custom working directory. A relative path is resolved against the directory the command is run from. (Default: `test_vibrChol1/` under the current working directory.)

**Step-specific flags:**

- `--metadata`: Download metadata.
- `--fastq`: Download FASTQ files.
- `--trna`: Download tRNA sequences.
- `--genome`: Download reference genome.
- `--trim`: Run adapter trimming.
- `--makedb`: Create tRNA database.
- `--map`: Run read mapping.
- `--build`: Build the AnnData object (no split analysis).
- `--split-build`: Build with read length split enabled (creates main + u60/o60 h5ad files).
- `--cluster`: Run clustering algorithms.
- `--graph`: Generate visualization plots (main h5ad only).
- `--split-graph`: Generate plots for split h5ad files (u60/o60).
- `--hubonly`: Generate UCSC track hubs without building the full database.

## Dataset Details

The suite uses data from:

- **Study**: Kimura et. al, 2020 ["Comparative tRNA sequencing and RNA mass spectrometry for surveying tRNA modifications"](https://www.nature.com/articles/s41589-020-0558-1)
- **Organism**: _Vibrio cholerae_
- **Accessions**: SRP254278 (Filtered subset)
- **Optimization**: FASTQ downloads are subsampled to the first 100,000 reads (`-X 100000`) to ensure rapid execution.

## Pipeline Breakdown

### 1. Data & Metadata Acquisition

- **Metadata**: Uses `pysradb` to fetch metadata for SRP254278. It specifically filters out non-_Vibrio_ entries (e.g., _E. coli_ controls) and specific gene knockouts to create a clean test dataset.
- **FASTQ**: Uses `prefetch` and `fastq-dump` to retrieve raw sequencing data based on the filtered metadata.

### 2. Reference Preparation

- **Genome**: Downloads the _Vibrio cholerae_ genome (GCF_000006745.1) from NCBI.
  - Converts GFF annotations to GTF using `gffread`.
  - Standardizes chromosome names (e.g., changing NCBI RefSeq IDs to `chrI`, `chrII`) in both FASTA and GTF files to match gtRNAdb conventions.
- **tRNAs**: Extracts pre-packaged tRNA sequences from the local `assets` directory (simulating a gtRNAdb download).

### 3. Preprocessing

- **Trimming**: Runs `preprocess trim` using `fastp` with specific adapter sequences to clean raw reads.
- **Database Generation**: Runs `preprocess makedb` to create a Bowtie2 index using the standardized genome and tRNA sequences (`-s bact` mode).
- **Mapping**: Runs `preprocess map` to align trimmed reads to the generated index, producing BAM files in `processed/vibrChol1/bam`.

### 4. Database Construction (Build)

Runs `trnagraph.py build` to aggregate the coverage data.

- **Inputs**: Uses the generated BAM files, the downloaded GTF for non-tRNA features, and a pairs file (`config/vibrChol1.pair.txt`).
- **Configuration**: Runs with the default multi-mapped-read filtering (tRAX-parity behavior; no flag needed). The `--split-build` test step re-runs build with `--readlengthsplit 60`, which splits BAM files by read length internally (into temporary `u60`/`o60` BAM files that are deleted once merged into the AnnData object) rather than as a separate preprocessing command.

### 5. Post-Processing & Visualization

- **Clustering**: Runs `trnagraph.py cluster` to perform dimensionality reduction (UMAP) and density-based clustering on the AnnData object.
- **Graphing**: Runs `trnagraph.py graph` to generate the full suite of visualizations (Heatmaps, PCA, Coverage, etc.) in the `vibrChol1/graphs` directory.
