# tRNAgraph Test Suite

This document outlines the automated test suite for the tRNAgraph project. The test suite, implemented in `toolsTestSuite.py`, is designed to validate the functionality of the entire tRNAgraph pipeline, from data acquisition to read mapping. It is based on the original tRAX tutorial notebook but has been integrated directly into the package to ensure legacy feature support and regression testing.

## Overview

The test suite automates the following steps:

1. **Metadata Retrieval**: Downloads experiment metadata from SRA.
2. **Data Acquisition**: Downloads FASTQ files for specific accessions.
3. **Reference Preparation**: Downloads and prepares tRNA sequences and the reference genome (Vibrio cholerae).
4. **Preprocessing**: Trims adapters from raw reads.
5. **Database Creation**: Builds a Bowtie2 index for the tRNA database.
6. **Alignment**: Maps reads to the generated database.

## Usage

The test suite is typically invoked via the `trnagraph.py` CLI.

```bash
python trnagraph.py tools test [options]
```

### Options

- `--all`: Run the complete pipeline from scratch (cleans existing test data first).
- `--metadata`: Only download metadata.
- `--fastq`: Only download FASTQ files.
- `--trna`: Only download tRNA sequences.
- `--genome`: Only download reference genome.
- `--trim`: Only run adapter trimming.
- `--makedb`: Only create the tRNA database.
- `--map`: Only map reads.
- `--cleanrun`: Clean up all generated files after the run completes.

## Dataset Details

The test suite uses RNAseq data from ["Comparative tRNA sequencing and RNA mass spectrometry for surveying tRNA modifications"](https://www.nature.com/articles/s41589-020-0558-1) by Kimura et. al, 2020 (GEO: [GSE147614](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147614)).

### Samples

The suite downloads a subset of samples (Vibrio cholerae) to ensure quick execution while still providing comprehensive coverage of the pipeline's features.

Accessions used:

- SRR11431928 - SRR11431937

**Note**: For testing purposes, the pipeline subsamples reads to the first 100,000 reads (`-X 100000`) to speed up execution.

## Pipeline Steps

### 1. Metadata & Fastq Download

Uses `pysradb` to fetch metadata and `fastq-dump` (via `sra-tools`) to download raw sequencing data.

### 2. Reference Genome & tRNA Database

- Downloads _Vibrio cholerae_ tRNA sequences from [gtRNAdb](http://gtrnadb.ucsc.edu/).
- Downloads the reference genome from NCBI.
- Converts GFF annotations to GTF format using `gffread`.
- Standardizes chromosome names (e.g., changing NCBI accessions to `chrI`, `chrII`) to match gtRNAdb conventions.

### 3. Trimming

Uses the internal `preprocess trim` command (wrapping `fastp` or `cutadapt`) to remove adapters.

### 4. Database Generation

Runs `preprocess makedb` to generate a Bowtie2 index from the fetched tRNA sequences and genome.

### 5. Mapping

Runs `preprocess map` to align the trimmed reads to the custom tRNA database, generating BAM files and coverage statistics.
