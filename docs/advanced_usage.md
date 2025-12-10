# Advanced Usage and Configuration

## Workflow Overview

The following diagram illustrates the comprehensive workflow of tRNAgraph, detailing the flow of data from raw inputs through preprocessing, analysis, and visualization.

```mermaid
graph TD
    %% Classes
    classDef input fill:#000000,stroke:#01579b,stroke-width:2px;
    classDef process fill:#000000,stroke:#2e7d32,stroke-width:2px;
    classDef storage fill:#000000,stroke:#ef6c00,stroke-width:2px;
    classDef output fill:#000000,stroke:#7b1fa2,stroke-width:2px;

    subgraph Inputs
        M[Manifest File]
        MD[Metadata File]
        FQ[FASTQ Files]
        G[Genome FASTA]
        T[tRNA Scan Output]
    end

    subgraph Preprocess
        direction TB
        DB[Make Database]
        TR[Trim Reads]
        MP[Map Reads]
        SP[Split BAMs]

        G & T --> DB
        FQ & M --> TR
        TR --> MP
        DB --> MP
        MP --> SP
    end

    subgraph Analyze
        direction TB
        B[Build AnnData]
        C[Cluster Data]
        MG[Merge Datasets]

        MP & MD --> B
        SP --> B
        B --> C
        B --> MG
    end

    subgraph Tools
        direction TB
        L2FC[Log2 Fold Change]
        CSV[Export CSV]

        B --> L2FC
        B --> CSV
    end

    subgraph Graph
        direction TB
        V[Generate Graphs]

        C --> V
        B --> V
        L2FC --> V
    end

    subgraph Outputs
        H5[AnnData .h5ad]
        H5C[Clustered .h5ad]
        FIG[Figures PDF/PNG]
        CSVO[CSV Files]
    end

    B --> H5
    C --> H5C
    V --> FIG
    CSV --> CSVO
    L2FC -.-> H5

    %% Apply Classes
    class M,MD,FQ,G,T input;
    class DB,TR,MP,SP,B,C,MG,L2FC,CSV,V process;
    class H5,H5C storage;
    class FIG,CSVO output;
```

## Workflow Steps

### 1. Preprocess

The preprocessing module handles raw data preparation.

- **[makedb](cli_reference.md#makedb)**: Creates a Bowtie2 index from a reference genome and tRNA predictions.
- **[trim](cli_reference.md#trim)**: Uses `fastp` to remove adapters and process UMIs from raw FASTQ files.
- **[map](cli_reference.md#map)**: Aligns trimmed reads to the tRNA database using Bowtie2.
- **[split](cli_reference.md#split)**: Splits BAM files based on read length (e.g., <60bp and >=60bp).

### 2. Analyze

The analysis module builds and refines the core database.

- **[build](cli_reference.md#build)**: Aggregates alignment data (BAMs) and metadata into a structured AnnData object (`.h5ad`).
- **[cluster](cli_reference.md#cluster)**: Performs dimensionality reduction (UMAP) and density-based clustering (HDBSCAN) on the dataset.
- **[merge](cli_reference.md#merge)**: Combines multiple AnnData objects into a single dataset.

### 3. Graph

- **[graph](cli_reference.md#graph)**: Generates a comprehensive suite of visualizations (Heatmaps, PCA, Coverage Plots, etc.) from the AnnData object.

### 4. Tools

Utility functions for specific data operations.

- **[log2fc](cli_reference.md#log2fc)**: Calculates Log2 Fold Change statistics for differential expression analysis.
- **[csv](cli_reference.md#csv)**: Exports the internal data structures (obs, var, X) to CSV format for external use.

## Configuration Files

You can control filtering and coloring in `trnagraph.py graph` using JSON files.

### Filter Configuration (`--config`)

Used to filter data before plotting.

```json
{
  "name": "analysis_subset",
  "obs": {
    "treatment": ["treatment_A"],
    "celltype": ["HEK293"]
  },
  "obs_r": {
    "amino": ["Und", "Sup"]
  },
  "var": {
    "coverage": ["unique"]
  }
}
```

- `name`: Name for the filtering configuration. Will be used to create a subfolder in the output directory.
- `obs`: Include only these values.
- `obs_r`: **Reverse** filter (Exclude these values).

### Custom Colormaps (`--colormap`)

Define custom colors for groups or features. Supports Hex, RGB, or Matplotlib names.

```json
{
  "group": {
    "Control": "lightgrey",
    "Treated": "#FF5733"
  },
  "amino": {
    "Ala": "royalblue",
    "Gly": "#FF9896"
  }
}
```

> [!NOTE]
> Some plots default to using group as the default category for plotting making a colormap with this name will override the default colormap in those cases.

## Python API & Downstream Analysis

Because tRNAgraph produces standard AnnData objects, you can easily load them into Python for custom analysis using `scanpy` or `pandas`.

### Loading and Basic Filtering

```python
import anndata as ad

# Load the database
adata = ad.read_h5ad("results/tRNAgraph.h5ad")

# Filter for specific coverage type and remove gaps
adata = adata[adata.var["coverage"] == "unique"]
adata = adata[adata.var["gap"] == False]

# Filter by Sample Group
adata_subset = adata[adata.obs["group"] == "Treatment_A"]
```

### Accessing Raw Data

```python
# Access normalized data
norm_data = adata.X

# Access raw counts
raw_data = adata.layers["raw"]
```

## Testing

A test suite is included to validate installations and ensure functionality. See the [Test Suite Documentation](testSuite.md) for details on running tests and interpreting results.
