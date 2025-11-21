# tRNAgraph Test Suite

This document outlines the steps taken to set up and run the test suite for the tRNAgraph project. The test suite is designed to validate the functionality of various tools and components within the tRNAgraph package.

## Dataset

RNAseq data is taken from ["Comparative tRNA sequencing and RNA mass spectrometry for surveying tRNA modifications"](https://www.nature.com/articles/s41589-020-0558-1) by Kimura et. al, 2020. This can be accessed via the GEO accession code [GSE147614](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147614) or its corresponding SRA value SRP254278.

[pysradb](https://github.com/saketkc/pysradb) is a tool to quickly pull SRA metadata for an experiment if given an SRA/ENA/GEO/Etc accession. We will use this to grab the metadata for the experiment and save it as `SRP254278.tsv`.

We only plan to download some of the samples from this experiment. Since the later samples are singleton, we will exclude them and only download the accessions listed below:

```
SRR11431928
SRR11431929
SRR11431930
SRR11431931
SRR11431932
SRR11431933
SRR11431934
SRR11431935
SRR11431936
SRR11431937
```

For this test, we will use `-X 100000` to subsample the fastq reads to the first 100,000 total reads so that this pipeline can run quickly. In typical use cases, you should **NOT** use this flag to subset your samples except for testing purposes.

## Creating tRNA Database

A bowtie2 tRNA database is required to run alignment. The files for this can be provided via [gtRNAdb](http://gtrnadb.ucsc.edu/index.html) or generated via [tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/). In the case of this test, we will use the [Vibr chol](https://gtrnadb.org/genomes/bacteria/Vibr_chol_O1_biovar_El_Tor_N16961/) database from gtRNAdb.

We will download and extract the tRNA names with `curl` and the reference genome for Vibr Chol. Sometimes, the nomenclature for the reference genome and gtRNAdb don't align to correct for this, we will also use a `sed` command to change the Chromosome names to match gtRNAdb.

Vibro Chol downloaded from NCBI only has a `gff` this needs to be converted into a `gtf` via `gffread`, and then the chromosome names will also need to be converted via `sed` commands.
