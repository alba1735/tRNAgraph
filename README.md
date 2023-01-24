# tRNAgraph
tRNAgraph is a tool for analyzing tRNA-seq data generated from tRAX. It can be used to create a AnnData database object from a tRAX coverage file, or to analyze an existing database object and generate expanded visulizations. The database object can also be used to perform further analysis beyond the scope of what tRAX can do.

## About
[tRAX](https://github.com/UCSC-LoweLab/tRAX) is a tool often used for analyzing tRNA-seq data. While it is generates a comprehensive set of results, it does not provide a way to visualize specfic meta-data associated with a particular experiment. tRNAgraph is a tool that can be used to create a database object from a tRAX coverage file containing various experimental conditions not captured by tRAX. The database object can then be used to generate a variety of visualizations, including heatmaps, coverage plots, pca plots, and more that are more specific to the experimental conditions of interest.

## Installation
Dependencies can be installed using conda:

```
conda env create -f environment.yml
```

## Usage
### Input files
While tRNAgraph will work with any tRAX generated coverage file, it is designed to work with sample names that follow a [specific format](http://trna.ucsc.edu/tRAX/#step-3-analyze-sequencing-data-for-gene-expression). By including the experimental conditions in the sample names and samples groups found in coulmns 1 and 2 of the tRAX generation file the resulting database object can be used to generate more specific visualizations. The sample names within the tRAX generation file should be in the format of:

```
<sample_name>_<condition-1>...<condition-n> <sample_group>_<condition-1>...<condition-n> <Fastq>
```

Each experimental condition should be sperated by an underscore (`_`) and an unlimited amount of conditions can be included.

To get the most utility out of tRNAgraph, it is recommended that either a list of observations or a file containing a list of observations be provided. The observations should match the cateogories found in the tRAX generation file. The observations can be provided in the form of a list or a file containing a list of observations. If a file is provided, each observation should be tab seperated.

Example tRAX generation file:

```
sample1_celltypeX_treatmentA_condition1 celltypeX_treatmentA_condition1 /path/to/fastq
sample2_celltypeX_treatmentA_condition2 celltypeX_treatmentA_condition2 /path/to/fastq
sample3_celltypeY_treatmentB_condition1 celltypeY_treatmentB_condition1 /path/to/fastq
...
```

Example list of observations:

```
celltype treatment condition
```

The following observations will automatically be generated as well and so using the following should be avoided:
`trna`, `iso`, `amino`, `sample`, `deseq2_sizefactor`, and `refseq`.


### Activating the environment

```
conda activate tRNAgraph
```

### Building a database object

The tRNAgraph.py script can be used to build a database object from a tRAX directory with the following command:
```bash
python tRNAgraph.py build -i <tRAX_directory> -o <output_file> -l <list_of_observations>
```

* `-i` or `--input` is the path to the tRAX directory you want to build a AnnData database object from.
* `-o` or `--output` is the path to the output file. The output file should be a `.h5ad` file. By default, the output file will be named `trnagraph.h5ad` if no path is provided.
* Observations are the metadata to include in the database object. This can be a list of observations `-l/--observationslist` or a tab seperated file containing a list of observations `-f/--observationsfile`. If no list is provided any metadata will be annotated `obs_#` where `#` is the ordered observation.

### Generating visualizations

The tRNAgraph.py script can be used to generate a variety of visualizations from a database object with the following command:
```bash
python tRNAgraph.py graph -i <input_database> -o <output_directory> -g <graph_type>
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
#### PCA

* `--pcamarkers` is the observation to use for choosing which markers to populate the pca plot. By default the samples will be used as markers.
* `--pcacolors` is the observation to use for coloring the PCA plot. If no observation is provided, the PCA plot will be colored by samples.