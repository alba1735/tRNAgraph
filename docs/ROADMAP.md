# Features to add

## General

* Add options to compute with RAW counts and alternative normalization methods
* Add a command to run both create and build in one command, skipping build if it already exists
* Add combination PDF plots to each graph type as needed

## trnagraph.py

* Allow for new DeSeq2 analysis to be run on the combined object
* Add validation that sample names match samples so it can actually write the adatafile

## pca_tools.py

* Add a function to plot the PCA with smallRNAs included

## correlation_tools.py

* Add a function to plot the correlation with smallRNAs

## coverage_tools.py

* add a function to sort the combined coverage from most to least abundant

## All plots

* Add a function to check if the parameter exisits in adata columns and if not default to samples

## radar_tools.py

* Add total amino dict not just the subset one I am using now also generate a list of all the amino acids from the adata object so it works universal to all species.
