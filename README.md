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

### Activating the environment
```
conda activate tRNAgraph
```

### Usage
```
trnagraph.py [-h] {build,graph} ...
    build:
        -i TRAX_COVERAGE, --trax_coverage TRAX_COVERAGE
                                tRAX coverage file
        -l OBSERVATION_LIST, --observationslist OBSERVATION_LIST
                                list of observations to include in the database
        -f OBSERVATION_FILE, --observationsfile OBSERVATION_FILE
                                file containing a list of observations to include in the
                                database
        -o OUTPUT, --output OUTPUT
                                output file name (default: tRNAgraph.h5ad)
    graph:
        -i DATABASE, --database DATABASE
                                database object
        -g GRAPH_TYPE, --graphtype GRAPH_TYPE
                                graph type
        -a GRAPH_ARGUMENTS, --grapharguments GRAPH_ARGUMENTS
                                graph arguments
        -o OUTPUT, --output OUTPUT
                                output file name (default: tRNAgraph/)
```
