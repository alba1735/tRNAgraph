# tRNAgraph

tRNAgraph is a tool for analyzing tRNA-seq data generated from tRAX. It can be used to create a AnnData database object from a tRAX coverage file, or to analyze an existing database object and generate expanded visulizations. The database object can also be used to perform further analysis beyond the scope of what tRAX can do.

## About

[tRAX](https://github.com/UCSC-LoweLab/tRAX) is a tool often used for analyzing tRNA-seq data. While it is generates a comprehensive set of results, it does not provide a way to visualize specfic meta-data associated with a particular experiment. tRNAgraph is a tool that can be used to create a database object from a tRAX coverage file containing various experimental conditions not captured by tRAX. The database object can then be used to generate a variety of visualizations, including heatmaps, coverage plots, pca plots, and more that are more specific to the experimental conditions of interest.

## Installation

Dependencies can be installed using conda:

```bash
conda env create -f environment.yml
```

## Input files

tRNAgraph will work with any tRAX generated coverage file however imputing the meta-data associated with the samples will increase the utility of the tool. You can choose how you wish to impute the meta-data associated with the samples by either automatically generating the meta-data from the sample names if you follow a specific format, or you can provide a meta-data file. By default, tRNAgraph will automatically generate the meta-data from the sample names regardless of how you name the samples.

### Automatic meta-data generation

If you want to generate the metadata automatically all you need to do is seperate each experimental condition or observation by an underscore (`_`) with in the tRAX [samples file](http://trna.ucsc.edu/tRAX/#step-3-analyze-sequencing-data-for-gene-expression) and incorporate them into the sample and group names.

```txt
<sample_name>_<condition-1>...<condition-n> <sample_group>_<condition-1>...<condition-n> <Fastq>
```

Including the experimental conditions in the sample names and samples groups found in coulmns 1 and 2 of the tRAX generation file the resulting database object can be used to generate more specific visualizations without a limit on the number of conditions you want to include. 

Example tRAX generation file:

```txt
sample1_celltypeX_treatmentA_condition1 celltypeX_treatmentA_condition1 /path/to/fastq
sample2_celltypeX_treatmentA_condition2 celltypeX_treatmentA_condition2 /path/to/fastq
sample3_celltypeY_treatmentB_condition1 celltypeY_treatmentB_condition1 /path/to/fastq
...
```

If you are using this automatic method to get the most utility out of tRNAgraph, it is recommended that either a list of observations (`-l/--observationslist`) or a file containing a list of observations (`-f/--observationsfile`) be provided. The observations should match the cateogories found in the tRAX generation file. The observations can be provided in the form of a list or a file containing a list of observations. If a file is provided, each observation should be tab seperated.

Example list of observations:

```txt
celltype treatment condition
```

If no list of observations is provided then all observations will be annotated `obs_#` where `#` is the ordered observation.

### Manual meta-data imputation

If you want to manually impute the meta-data associated with the samples you can provide a tab seperated file (`-m/--metadatafile`) containing the sample names, sample groups, and any meta-data associated with the samples. The file should be formatted as follows:

```tsv
sample1 sampleGroup1 celltypeX treatmentA condition1
sample2 sampleGroup1 celltypeX treatmentA condition2
sample3 sampleGroup2 celltypeY treatmentB condition1
```

It is recommended to also use the `-l/--observationslist` or `-f/--observationsfile` options when using this method.

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

* `-i` or `--input` is the path to the tRAX directory you want to build a AnnData database object from.
* `-s` or `--samples` is the path to the tRAX generation file containing sample names, sample groups, and fastq paths.
* `-o` or `--output` is the path to the output file. The output file should be a `.h5ad` file. By default, the output file will be named `trnagraph.h5ad` if no path is provided.
* Observations are the metadata to include in the database object. This can be a list of observations `-l/--observationslist` or a tab seperated file containing a list of observations `-f/--observationsfile`. If no list is provided any metadata will be annotated `obs_#` where `#` is the ordered observation.

### Generating visualizations

The tRNAgraph.py script can be used to generate a variety of visualizations from a database object with the following command:

```bash
python trnagraph.py graph -i <input_database> -o <output_directory> -g <graph_type> <graph_parameters>
```

* `-i` or `--input` is the path to the database object you want to generate visualizations from.
* `-o` or `--output` is the path to the output directory. By default, the output directory will be named `trnagraph` if no path is provided.
* `-g` or `--graph` is the type of graph to generate. The following graphs can be generated:
  * `coverage` - Generates a coverage plots.
  * `heatmap` - Generates a heatmaps of the differantial tRNA expression.
  * `pca` - Generates PCA plots.
  <!-- * `tsne` - Generates a tSNE plot of the tRNA coverage for each sample. -->
  <!-- * `umap` - Generates a UMAP plot of the tRNA coverage for each sample. -->
  * `volcano` - Generates a volcano plot of differantial tRNA expression.
  * `radar` - Generates a radar plot of the tRNA coverage.
  * `all` - Generates all of the above graphs.

The following parameters are optional and can be used to customize the graphs:

<!-- * `-c` or `--color` is the observation to use for coloring the graph. If no observation is provided, the graph will be colored by sample group. -->
#### Coverage

* `--coveragegrp` is the observation to use for grouping the coverage plots. If no observation is provided, the coverage plots will be grouped by sample group.
* `--coverageobs` is the observation to use for coloring the coverage plots. If no observation is provided, the coverage plots will be colored by sample group. Multiple observations can be provided.
* `--coveragetype` is the type of tRNA coverage to plot, by default it will plot unique coverage. All posible coverage types match the [tRAX coverage types](http://trna.ucsc.edu/tRAX/outputs/#abundance-of-trna-tdrs-and-other-genes).
* `--coveragegap` is wether to include tRNA gaps in the coverage plots. By default, tRNA gaps will be skipped in the coverage plots.
* `--coveragefill` is wether to fill the area under the coverage plots `fill` or include a confidence interval `ci`. By default, a confidence interval will be included in the coverage plots.

#### PCA

* `--pcamarkers` is the observation to use for choosing which markers to populate the pca plot. By default the samples will be used as markers.
* `--pcacolors` is the observation to use for coloring the PCA plot. If no observation is provided, the PCA plot will be colored by samples.

#### Correlation

* `--corrmethod` is the method to use for calculating the correlation. The following methods can be used:
  * `pearson` - Pearson correlation coefficient.
  * `spearman` - Spearman rank correlation.
  * `kendall` - Kendall Tau correlation coefficient.
