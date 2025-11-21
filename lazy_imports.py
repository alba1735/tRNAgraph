import importlib

class LazyLoader:
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
plotsLegacyReadLength = LazyLoader('plotsLegacyReadLength')
plotsLegacyMismatch = LazyLoader('plotsLegacyMismatch')
plotsTrimmingStats = LazyLoader('plotsTrimmingStats')
plotsLegacyLocusCoverage = LazyLoader('plotsLegacyLocusCoverage')
plotsLegacyMismatchBoxplot = LazyLoader('plotsLegacyMismatchBoxplot')
plotsLegacyGeneFeatures = LazyLoader('plotsLegacyGeneFeatures')
plotsLegacyFeatureTypes = LazyLoader('plotsLegacyFeatureTypes')
plotsLegacyPCA = LazyLoader('plotsLegacyPCA')
plotsLegacyScatter = LazyLoader('plotsLegacyScatter')
plotsLegacyCCA = LazyLoader('plotsLegacyCCA')
plotsLegacyVolcano = LazyLoader('plotsLegacyVolcano')
plotsLegacyCoverage = LazyLoader('plotsLegacyCoverage')

# Cluster modules
umap = LazyLoader('umap')
hdbscan = LazyLoader('hdbscan')

# Tools modules
toolsMap = LazyLoader('toolsMap')
toolsTDatabase = LazyLoader('toolsTDatabase')
toolsTrim = LazyLoader('toolsTrim')
toolsTG = LazyLoader('toolsTG')
