# Features to add

## General

* Add options to compute with RAW counts and alternative normalization methods
* Argparse updates [stackoverflow thread](https://stackoverflow.com/questions/7498595/python-argparse-add-argument-to-multiple-subparsers)
* Add multiprocessing to speed up the analysis via the graph generation
* Add a command to run both create and build in one command, skipping build if it already exists
* Add combination PDF plots to each graph type as needed

## trnagraph.py

* Add a class to concatenate multiple adata objects together
  * This will be useful for cross analysis of multiple experiments, different species, etc.
  * Allow for new DeSeq2 analysis to be run on the combined object
* Add validation that sample names match samples so it can actually write the adatafile

## pca_tools.py

* Add a function to plot the PCA with smallRNAs included
* `tsne` - Generates a tSNE plot of the tRNA coverage for each sample. -->
* `umap` - Generates a UMAP plot of the tRNA coverage for each sample. -->

## correlation_tools.py

* Add a function to plot the PCA with smallRNAs included

## All plots

* Add a function to check if the parameter exisits in adata columns and if not default to samples

## coverage_tools.py

* Move low coverage <20 reads to a seperate low coverage folder based on max coverage?

## radar_tools.py

* Add total amino dict not just the subset one I am using now also generate a list of all the amino acids from the adata object so it works universal to all species.
