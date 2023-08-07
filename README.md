# tRNAgraph

tRNAgraph is a tool for analyzing tRNA-seq data generated from tRAX. It can be used to create an AnnData database object from a tRAX coverage file or to analyze an existing database object and generate expanded visualizations. The database object can also be used to perform further analysis beyond the scope of what tRAX can do.

## About

[tRAX](https://github.com/UCSC-LoweLab/tRAX) is a tool often used for analyzing tRNA-seq data. While it generates a comprehensive set of results, it does not provide a way to visualize specific meta-data associated with a particular experiment. tRNAgraph is a tool that can be used to create a database object from a tRAX coverage file containing various experimental conditions not captured by tRAX. The database object can then be used to generate a variety of visualizations, including heatmaps, coverage plots, PCA plots, and more that are more specific to the experimental conditions of interest.

## Installation

Dependencies can be installed using conda:

```bash
conda env create -f environment.yml
```

## Input files

tRNAgraph will work with any tRAX (Python3 Version) generated coverage file however imputing the meta-data associated with the samples will increase the utility of the tool.

### Metadata

You will need to impute the meta-data associated with the samples if you want to generate graphs based on specific experimental conditions. To do this you can provide a whitespace/.csv/.tsv file (`-m/--metadatafile`) containing the sample names, sample groups, and any meta-data associated with the samples. If you wish to run the tool without providing a metadata file, you can instead provide the [tRAX samples file](http://trna.ucsc.edu/tRAX/#step-3-analyze-sequencing-data-for-gene-expression) used to generate your tRAX run. An example meta-data file is shown below:

```tsv
sample1 sampleGroup1 celltypeX treatmentA condition1
sample2 sampleGroup1 celltypeX treatmentA condition2
sample3 sampleGroup2 celltypeY treatmentB condition1
```

It is recommended that either a list of observations (`-l/--observationslist`) or a file containing a list of observations (`-f/--observationsfile`) be provided. The observations should match the categories found in the tRAX generation file, not including the `sample` and `group` columns. The observations can be provided in the form of a list or a file containing a list of observations. If a file is provided, each observation should be tab separated.

Example list of observations:

```tsv
celltype treatment condition
```

If no list of observations is provided then all observations will be annotated automatically `obs_#` where `#` is the ordered observation.

## Usage

### Activating the environment

```bash
conda activate trnagraph
```

### Building a database object

The tRNAgraph.py script can be used to build a database object from a tRAX directory with the following command:

```bash
python trnagraph.py build -i <tRAX_directory> -s <tRAX_samples_file> -o <output_file> -l <list_of_observations>
```

* `-i` or `--input` is the path to the tRAX directory you want to build an AnnData database object from.
* `-m` or `--metadata` is the path to the tRAX generation file containing sample names, sample groups, and fastq paths.
* `-o` or `--output` is the path to the output file. The output file should be a `.h5ad` file. By default, the output file will be named `trnagraph.h5ad` if no path is provided.
* Observations are the metadata categories to include in the database object. This can be a list `-l/--observationslist` or a tab separated file containing a list of observations `-f/--observationsfile`. If no list is provided any metadata will be annotated `obs_#` where `#` is the ordered observation.

### Generating visualizations

The tRNAgraph.py script can be used to generate a variety of visualizations from a database object with the following command:

```bash
python trnagraph.py graph -i <input_database> -o <output_directory> -g <graph_type> <graph_parameters>
```

* `-i` or `--input` is the path to the database object you want to generate visualizations from.
* `-o` or `--output` is the path to the output directory. By default, the output directory will be named `trnagraph` if no path is provided.
* `-g` or `--graph` is the type of graph to generate. The following graphs can be generated:
  * `coverage` - Generates coverage plots.
  * `heatmap` - Generates heatmaps of the differential tRNA expression.
  * `pca` - Generates PCA plots.
  * `correlation` - Generates correlation plots of the tRNA coverage.
  * `volcano` - Generates volcano plots of differential tRNA expression.
  * `radar` - Generates radar plots of the tRNA coverage.
  * `all` - Generates all of the above plots.
* `--config` is an optional flag to the path to a JSON file containing additonal graph parameters. [See configuration section below for more details.](#configuration)
* `-n` or `--threads` is the number of threads to use for generating the graphs. By default, the number of threads will be set to 1. This is mostly useful for generating coverage plots as they can take a long time to generate.
* `--log` is wether to output a log of the shell commands used to generate the graphs. By default, the log will not be output.

The following parameters are optional and can be used to customize the graphs:

#### Coverage

* `--coveragegrp` is the observation to use for grouping the coverage plots. If no observation is provided, the coverage plots will be grouped by sample group.
* `--coverageobs` is the observation to use for coloring the coverage plots. If no observation is provided, the coverage plots will be colored by sample group. Multiple observations can be provided.
* `--coveragetype` is the type of tRNA coverage to plot, by default it will plot unique coverage. All posible coverage types match the [tRAX coverage types](http://trna.ucsc.edu/tRAX/outputs/#abundance-of-trna-tdrs-and-other-genes).
* `--coveragegap` is whether to include tRNA gaps in the coverage plots. By default, tRNA gaps will be skipped in the coverage plots.
* `--combineonly` is whether to only plot the combined coverage plots. By default, the coverage plots will be plotted for each individual tRNA as well as the combined coverage plots.

#### Heatmap

* `--heatgrp` is the observation to use for grouping the heatmap plots. If no observation is provided, the heatmap plots will be grouped by sample group.
* `--heatrts` is the read types to use for the heatmap plots. By default, the heatmap plots will be generated using the total unique normalized reads.
* `--heatcutoff` is the cutoff for reads to include in the heatmap plots. By default, the heatmap plots will discard anything with less than 80 reads.
* `--heatbound` is the range to bound the heatmap plots. By default, the heatmap plots will be bounded by the top 25 and the bottom 25 values in the data for each comparative column.

#### PCA

* `--pcamarkers` is the observation to use for choosing which markers to populate the pca plot. By default, the samples will be used as markers.
* `--pcacolors` is the observation to use for coloring the PCA plot. If no observation is provided, the PCA plot will be colored by samples.

#### Correlation

* `--corrmethod` is the method to use for calculating the correlation. The following methods can be used:
  * `pearson` - Pearson correlation coefficient.
  * `spearman` - Spearman rank correlation.
  * `kendall` - Kendall Tau correlation coefficient.

#### Volcano

* `--volgrp` is the observation to use for grouping the volcano plots. If no observation is provided, the volcano plots will be grouped by sample group.
* `--volrt` is the read type to use for the volcano plots. By default, the volcano plots will be generated using the total unique normalized reads.
* `--volcutoff` is the cutoff for reads to include in the volcano plots. By default, the volcano plots will discard anything with less than 80 reads.

#### Radar

* `--radargrp` is the observation to use for grouping the radar plots. If no observation is provided, the radar plots will be grouped by sample group.

## Database Variables

The database object aggregates all the information from the tRAX directory into a single object allowing for easy calls to the data. In addition to being able to use the following variables as flags for figure generation, they can be used for further analysis and data manipulation independent of tRNAgraph.

### Observations

The observations are the metadata categories that are used to group and color the plots derived from the provided metadata. The observations are stored in the `obs` attribute of the database object as a pandas dataframe. The following observations are automatically added to the database object:

* `trna` - The tRNA name.
* `iso` - The tRNA isotype group.
* `amino` - The tRNA amino acid group.
* `sample` - The sample name from tRAX.
* `group` - The sample group tRAX.
* Any metadata categories provided in the observations list/file.
* The uniquely mapped reads are the reads that map to a single tRNA via alignment and filtering in tRAX and can be broken down into the following categories:
  * `nreads_whole_unique_raw` - The raw number of uniquely mapped whole-counts in the sample.
  * `nreads_whole_unique_norm` - The normalized number of uniquely mapped whole-counts in the sample.
  * `nreads_fiveprime_unique_raw` - The raw number of uniquely mapped 5' end counts in the sample.
  * `nreads_fiveprime_unique_norm` - The normalized number of uniquely mapped 5' end counts in the sample.
  * `nreads_threeprime_unique_raw` - The raw number of uniquely mapped 3' end counts in the sample.
  * `nreads_threeprime_unique_norm` - The normalized number of uniquely mapped 3' end counts in the sample.
  * `nreads_other_unique_raw` - The raw number of uniquely mapped other counts in the sample.
  * `nreads_other_unique_norm` - The normalized number of uniquely mapped other counts in the sample.
  * `nreads_total_unique_raw` - The raw number of uniquely mapped total counts in the sample. This is the sum of the above categories.
  * `nreads_total_unique_norm` - The normalized number of uniquely mapped total counts in the sample. This is the sum of the above categories.
* The multi-mapped reads are the reads that map via alignment but ambiguous reads are randomly assigned and can be broken down into the following categories:
  * `nreads wholecounts_raw` - The raw number of whole-counts in the sample.
  * `nreads wholecounts_norm` - The normalized number of whole-counts in the sample.
  * `nreads fiveprime_raw` - The raw number of 5' end counts in the sample.
  * `nreads fiveprime_norm` - The normalized number of 5' end counts in the sample.
  * `nreads threeprime_raw` - The raw number of 3' end counts in the sample.
  * `nreads threeprime_norm` - The normalized number of 3' end counts in the sample.
  * `nreads other_raw` - The raw number of other counts in the sample.
  * `nreads other_norm` - The normalized number of other counts in the sample.
  * `nreads total_raw` - The raw number of total counts in the sample. This is the sum of the above categories.
  * `nreads total_norm` - The normalized number of total counts in the sample. This is the sum of the above categories.
  * The following additional categories are also found under multi-mapped reads:
  * `nreads_wholeprecounts_raw` - The raw number of whole-precounts in the sample.
  * `nreads_wholeprecounts_norm` - The normalized number of whole-precounts in the sample.
  * `nreads_partialprecounts_raw` - The raw number of partial-precounts in the sample.
  * `nreads_partialprecounts_norm` - The normalized number of partial-precounts in the sample.
  * `nreads_trailercounts_raw` - The raw number of trailer-counts in the sample.
  * `nreads_trailercounts_norm` - The normalized number of trailer-counts in the sample.

### Variables

The variables are the metadata categories that are used to filter the read coverage. The variables are stored in the `var` attribute of the database object as a pandas dataframe. The following variables are automatically added to the database object:

* `gap` - Whether a position is a gap in canonical Sprinzl positions. These gaps are skipped in the coverage plots so that they can be easier to interpret when comparing different tRNAs.
* `positions` - The canonical Sprinzl positions.
* `coverage` - The coverage type matching coverage types found in the tRAX coverage file.

Since all coverage types are stored in the database object it is useful to specify which coverage type you want to use if you plan on using the database object for further analysis.

### Downstream Analysis and Filtering

Using the database object for further analysis is easy and follows the same syntax as using the AnnData object. For example, the following code will filter the database object to only include the coverage type `unique` and drop the `gap` positions:

```python
import anndata as ad

adata = ad.read_h5ad("tRNAgraph.h5ad")
adata = adata[adata.var["coverage"] == "unique"]
adata = adata[adata.var["gap"] == False]
```

If you wanted to filter the database object further to only include samples from group `A` and tRNA `tRNA-Ala-AGC-1` you could use the following code:

```python
adata = adata[adata.obs["group"] == "A"]
adata = adata[adata.obs["trna"] == "tRNA-Ala-AGC-1"]
```

The resulting table can be called using `adata.X`.

### Configuration

A JSON file can be used for complicated filtering and grouping of the data as well as using custom colormaps. Any of these categories can be left blank and they will be skipped.

* `name` - Is a name for the filtering configuration and will be saved as a subfolder in the output directory.
* `obs` and `var` -  Are conditions to filter on in the AnnData observation and variable categories respectively. The values can be a single value or a list of values. If a list of values is provided, the data will be filtered to include only the values in the list. If no values are provided, the data will not be filtered on that category.
* `colormap` - It is a dictionary of dictionaries that allows custom colormaps to be used for the observations. The first level of the dictionary is the observation category, and the second level is the value of color. The value can be a hex color code or an RGB tuple value. If no colormap is provided, the default colormap will be used. The colormap will only be used if the observation for the colormap is selected, generating the plot. For example, if coverage plots are being generated, the colormap will only be used if the `--coveragegrp` flag matches an existing colormap. The JSON file should be formatted as follows:

```json
{
    "name": "name",
    "obs": {
        "observation_1": ["value1", "value2"],
        "observation_2": ["value1"]
    },
    "var": {
        "variable_1": ["value1", "value2"]
    },
    "colormap": {
        "observation_1": {
            "value1": "#000000",
            "value2": "#FFFFFF"
        }
    }
}
```