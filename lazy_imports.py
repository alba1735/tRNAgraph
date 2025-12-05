import importlib

class LazyLoader:
    """Lazy loader for internal tRNAgraph submodules."""
    def __init__(self, module_name):
        self.module_name = module_name
        self.module = None

    def __getattr__(self, name):
        if self.module is None:
            self.module = importlib.import_module(f'.{self.module_name}', package='tRNAgraph')
        return getattr(self.module, name)


class ExternalLazyLoader:
    """Lazy loader for external packages (not part of tRNAgraph)."""
    def __init__(self, module_name):
        self.module_name = module_name
        self.module = None

    def __getattr__(self, name):
        if self.module is None:
            self.module = importlib.import_module(self.module_name)
        return getattr(self.module, name)


# Plot modules
plotsBar = LazyLoader('plotsBar')
plotsCount = LazyLoader('plotsCount')
plotsCluster = LazyLoader('plotsCluster')
plotsCompare = LazyLoader('plotsCompare')
plotsCorrelation = LazyLoader('plotsCorrelation')
plotsCoverage = LazyLoader('plotsCoverage')
plotsHeatmap = LazyLoader('plotsHeatmap')
plotsSeqlogo = LazyLoader('plotsSeqlogo')
plotsPca = LazyLoader('plotsPca')
plotsRadar = LazyLoader('plotsRadar')
plotsVolcano = LazyLoader('plotsVolcano')

# Legacy plot modules
plotsTrimmingStats = LazyLoader('plotsTrimmingStats')

# External cluster modules (umap-learn and hdbscan packages)
umap = ExternalLazyLoader('umap')
hdbscan = ExternalLazyLoader('hdbscan')

# Tools modules
toolsMap = LazyLoader('toolsMap')
toolsTDatabase = LazyLoader('toolsTDatabase')
toolsTrim = LazyLoader('toolsTrim')
toolsTG = LazyLoader('toolsTG')
toolsCountReads = LazyLoader('toolsCountReads')
toolsGetCoverage = LazyLoader('toolsGetCoverage')
toolsTrackHub = LazyLoader('toolsTrackHub')
