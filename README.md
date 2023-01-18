# tRNAGraph
tRNAgraph is a tool for analyzing tRNA-seq data generated from tRAX. It can be used to create a AnnData database object from a trax coverage file, or to analyze an existing database object and generate expanded visulizations. The database object can also be used to perform further analysis beyond the scope of what tRAX can do.

## Installation
Dependencies can be installed using conda:
```
conda env create -f environment.yml
```

## Usage
### Activating the environment
```
conda activate tRNAgraph
```

### Usage
```
trnagraph.py [-h] {build,graph} ...
    build:
        -tc TRAX_COVERAGE, --trax_coverage TRAX_COVERAGE
                                trax coverage file
        -obs OBSERVATION_LIST, --observation_list OBSERVATION_LIST
                                list of observations to include in the database
        -obsf OBSERVATION_FILE, --observation_file OBSERVATION_FILE
                                file containing a list of observations to include in the
                                database
        -o OUTPUT, --output OUTPUT
                                output file name (default: tRNAgraph.h5ad)
    graph:
        -db DATABASE, --database DATABASE
                                database object
        -g GRAPH_TYPE, --graph_type GRAPH_TYPE
                                graph type
        -ga GRAPH_ARGUMENTS, --graph_arguments GRAPH_ARGUMENTS
                                graph arguments
        -o OUTPUT, --output OUTPUT
                                output file name (default: tRNAgraph/)
```
