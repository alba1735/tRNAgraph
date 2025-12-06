# Advanced Usage and Configuration

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
