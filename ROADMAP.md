# Features to add

## General

* Add options to compute with RAW counts and alternative normalization methods

## trnagraph.py

* Clean up the uniqefeat column and sorting function
* Add a function to get observations from a file rather than inherit them from the trax sample names
* Add a flag to allow for simple sample/group names
* Add a class to concatenate multiple adata objects together
  * This will be useful for cross analysis of multiple experiments, different species, etc.
  * Allow for new DeSeq2 analysis to be run on the combined object

## pca_tools.py

* Add a function to plot the PCA with smallRNAs included

## correlation_tools.py

* Add a function to plot the PCA with smallRNAs included

## All plots
* Add a function to check if the parameter exisits in adata columns and if not default to samples