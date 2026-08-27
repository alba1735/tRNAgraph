import importlib
from typing import Any

class LazyLoader:
    """Lazy loader for internal tRNAgraph submodules."""
    def __init__(self, module_name: str):
        self.module_name = module_name
        self.module = None

    def __getattr__(self, name: str) -> Any:
        if self.module is None:
            self.module = importlib.import_module(f'.{self.module_name}', package='trnagraph.modules')
        return getattr(self.module, name)


class ExternalLazyLoader:
    """Lazy loader for external packages (not part of tRNAgraph)."""
    def __init__(self, module_name: str):
        self.module_name = module_name
        self.module = None

    def __getattr__(self, name: str) -> Any:
        if self.module is None:
            self.module = importlib.import_module(self.module_name)
        return getattr(self.module, name)


# Plot modules

plotsCount = LazyLoader('plotsCount')
plotsCluster = LazyLoader('plotsCluster')
plotsCompare = LazyLoader('plotsCompare')
plotsCorrelation = LazyLoader('plotsCorrelation')
plotsCoverage = LazyLoader('plotsCoverage')
plotsHeatmap = LazyLoader('plotsHeatmap')
plotsMismatch = LazyLoader('plotsMismatch')
plotsSeqlogo = LazyLoader('plotsSeqlogo')
plotsPca = LazyLoader('plotsPca')
plotsRadar = LazyLoader('plotsRadar')
plotsVolcano = LazyLoader('plotsVolcano')

# Legacy plot modules
plotsTrimmingStats = LazyLoader('plotsTrimmingStats')

# External cluster modules (umap-learn and hdbscan packages)
umap = ExternalLazyLoader('umap')
hdbscan = ExternalLazyLoader('hdbscan')
anndata = ExternalLazyLoader('anndata')
matplotlib = ExternalLazyLoader('matplotlib')

# Tools modules
toolsMap = LazyLoader('toolsMap')
toolsDedup = LazyLoader('toolsDedup')
toolsTDatabase = LazyLoader('toolsTDatabase')
toolsTrim = LazyLoader('toolsTrim')
toolsTG = LazyLoader('toolsTG')
toolsCountReads = LazyLoader('toolsCountReads')
toolsGetCoverage = LazyLoader('toolsGetCoverage')
toolsTrackHub = LazyLoader('toolsTrackHub')
toolsTestSuite = LazyLoader('toolsTestSuite')
toolsUpdate = LazyLoader('toolsUpdate')

# Adata modules
adataGraph = LazyLoader('adataGraph')
adataMerge = LazyLoader('adataMerge')
adataCluster = LazyLoader('adataCluster')
adataBuild = LazyLoader('adataBuild')
