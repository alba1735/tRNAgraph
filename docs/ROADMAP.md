# Features to add

## General

* Add options to compute with RAW counts and alternative normalization methods

## trnagraph.py

* Allow for new DeSeq2 or other analysis to be run on the combined object
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

## Priority

* Fix the acceptor stems not being labeled in adata.var and other missing info with those positions
  * Shows up in CVS exports
* Move all flag validation in adata obs or var into the trnagraph.py script from the graphing scripts
* Prevent merging of adata objects with different trax/trnagraph runinfo


<!-- # Verify adata is valid for chosen coverage group or obs
#if adata.obs[self.coverage_grp].isna().any():
#    raise ValueError('Coverage group contains NaN values.\nThis most likely means that you forgot to include samples in your metadata file that are present in your trax directory.\n' \
#                     'Try adding the samples to your metadata file and rebuilding the AnnData object or selecting a different coverage group.')
# if self.coverage_obs:
#     if adata.obs[self.coverage_obs].isna().any().any():
#         raise ValueError('Coverage obs contains NaN values.\nThis most likely means that you forgot to include samples in your metadata file that are present in your trax directory.\n' \
#                         'Try adding the samples to your metadata file and rebuilding the AnnData object or selecting a different coverage obs.') -->