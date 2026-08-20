import anndata as ad
import os
import multiprocessing
import json
import warnings
from pydantic import ValidationError
warnings.filterwarnings("ignore", message="Attempting to set identical low and high ylims")
from . import toolsTG
from .toolsSchemas import GraphFilterConfig, ColormapFile
from .lazy_imports import (
    plotsCount, plotsCluster, plotsCompare, plotsCorrelation,
    plotsCoverage, plotsHeatmap, plotsSeqlogo, plotsPca, plotsRadar, plotsVolcano
)

class anndataGrapher:
    '''
    Class to generate graphs from an AnnData object by calling the appropriate graphing functions
    '''
    def __init__(self, args):
        self.args = args
        self.adata_original = ad.read_h5ad(self.args.anndata)
        # Resolve the requested normalization:split-tag ONCE, into a working copy, so every
        # downstream plot module reads .X/.obs[...]/.uns[...] exactly as it does for the
        # full/default variant -- see toolsTG.build_variant_view() for why this resolved copy
        # must never be written back to self.args.anndata directly.
        self.variant_spec = toolsTG.parse_variant(self.adata_original, getattr(self.args, 'variant', 'norm:full'))
        self.adata = toolsTG.build_variant_view(self.adata_original, self.variant_spec)
        self.config_name = 'default'
        # Resolve grouping-column args before anything below reads them -- the log2FC
        # precompute a few lines down reads self.args.heatgrp/volgrp directly, ahead of the
        # plots*.py modules that separately validate their own grp argument.
        self._resolve_grp_args()
        # Load cmap dict for each graph type
        self.cmap_dict = {'cluster':self.args.clustergrp, 'compare':self.args.comparegrp1, \
                          'coverage':self.args.covgrp, 'pca':self.args.pcacolors, 'radar':self.args.radargrp, \
                          'volcano':self.args.volgrp}
        # Load all graph types if specified
        if self.args.graphtypes == 'all' or 'all' in self.args.graphtypes:
            self.args.graphtypes = ['cluster', 'correlation', 'count', 'coverage', 'heatmap', 'logo', 'pca', 'radar', 'volcano']
            self.args.clusteroverview = True
        # Load max threads available unless specified
        if self.args.threads == 0:
            try:
                # This is a linux only function but is less likely to cause problems than multiprocessing.cpu_count()
                self.args.threads = len(os.sched_getaffinity(0))
            except:
                self.args.threads = multiprocessing.cpu_count()
        # Load config file if specified
        if self.args.config:
            print('Loading config file: ' + self.args.config)
            with open(self.args.config, 'r') as f:
                raw_config = json.load(f)
            try:
                config = GraphFilterConfig.model_validate(raw_config)
            except ValidationError as e:
                raise ValueError(f'Invalid config file {self.args.config}:\n{e}') from e
            self.args.output += '/' + config.name
            self.config_name = config.name
            print(toolsTG.builder(self.args.output))
            if config.obs is not None or config.obs_r is not None:
                # Dictionary of uns columns and values to filter by as groups and samples since the coulmns are different from the main adata obs
                obs_dict = {i:True for i in self.adata.uns['amino_counts'].columns.values}
                obs_dict.update({i:True for i in self.adata.uns['type_counts'].columns.values})
                filter_dict = dict(config.obs or {})
                if config.obs_r is not None:
                    # Add the inverse of the obs_r filter to the filter_dict
                    for k,v in config.obs_r.items():
                        filter_dict[k] = [i for i in self.adata.obs[k].unique() if i not in v]
                for k,v in filter_dict.items():
                    print('Filtering AnnData object by observation: ' + k + ' , ' + str(v))
                    # Filter all uns columns by the observation and update the obs_dict
                    sub_obs_dict = dict(zip(self.adata.obs['sample'], self.adata.obs[k]))
                    sub_obs_dict.update(dict(zip(self.adata.obs['group'], self.adata.obs[k])))
                    obs_dict = {i:False if sub_obs_dict.get(i,False) not in v else True for i,j in obs_dict.items()}
                    # Filter the adata object by the observation
                    self.adata = self.adata[self.adata.obs[k].isin(v), :]
                # Filter uns columns by the obs_dict
                uns_dict = self.adata.uns.copy() # This was the only way I could get the uns columns to update without an implicit copy warning
                for uns_key in ['amino_counts', 'anticodon_counts', 'nontRNA_counts', 'type_counts']:
                    uns_value = self.adata.uns[uns_key].loc[:, [i for i in self.adata.uns[uns_key].columns.values if obs_dict[i]]].copy()
                    uns_dict[uns_key] = uns_value
                self.adata.uns = uns_dict
            if config.var is not None or config.var_r is not None:
                filter_dict = dict(config.var or {})
                if config.var_r is not None:
                    # Add the inverse of the var_r filter to the filter_dict
                    for k,v in config.var_r.items():
                        filter_dict[k] = [i for i in self.adata.var[k].unique() if i not in v]
                for k,v in filter_dict.items():
                    print('Filtering AnnData object by variable: ' + k + ' , ' + str(v))
                    self.adata = self.adata[:, self.adata.var[k].isin(v)]
            print('Config file loaded.\n')

        else:
            self.args.config = {}
        # Load the colormap if specified
        if self.args.colormap:
            print('Loading colormap file: ' + self.args.colormap)
            with open(self.args.colormap, 'r') as f:
                raw_colormap = json.load(f)
            try:
                self.args.colormap = ColormapFile.model_validate(raw_colormap).root
            except ValidationError as e:
                raise ValueError(f'Invalid colormap file {self.args.colormap}:\n{e}') from e
            print('Colormap loaded.\n')
        # Check for heatmap or volcano in graph types and if present check for readcount_cutoff in log2FC - Will precompute this and save it back to uns
        # This is done now to prevent saving issues later if multiprocessing is used
        log2FC_dict = self.adata.uns['log2FC'].copy()
        if 'heatmap' in self.args.graphtypes or 'volcano' in self.args.graphtypes:
            for grp in list(set([self.args.heatgrp, self.args.volgrp])):
                for readtype in [f'nreads_{i}_norm' for i in self.args.diffrts]: #list(set(self.args.heatrts+self.args.volrts))]:
                    for cutoff in list(set([self.args.heatcutoff, self.args.volcutoff])):
                        toolsTG.adataLog2FC(self.adata, grp, readtype, readcount_cutoff=cutoff, config_name=self.config_name, overwrite=self.args.regen_uns).main()
        if 'volcano' in self.args.graphtypes:
            # The volcano combined overview page always uses these two read types (mirroring
            # PCA's default --pcareadtypes), regardless of what --diffrts requests.
            for readtype in ['nreads_total_unique_norm', 'nreads_total_norm']:
                toolsTG.adataLog2FC(self.adata, self.args.volgrp, readtype, readcount_cutoff=self.args.volcutoff, config_name=self.config_name, overwrite=self.args.regen_uns).main()
        if log2FC_dict != self.adata.uns['log2FC'] or self.args.regen_uns:
            print('The log2FC uns dictionary has been updated.\n')
            # Persist onto the ORIGINAL (unresolved) adata, into the correct namespaced
            # location -- never write self.adata (the resolved view) back to disk, since for a
            # split variant it would overwrite the real full/default variant's data.
            if self.variant_spec.tag == 'full':
                self.adata_original.uns['log2FC'] = self.adata.uns['log2FC']
            else:
                self.adata_original.uns.setdefault('size_splits', {}).setdefault(self.variant_spec.tag, {})['log2FC'] = self.adata.uns['log2FC']
            self.adata_original.write(self.args.anndata)

    def _resolve_grp_args(self):
        '''
        Validate the CLI's grouping-column parameters against the resolved AnnData object
        once, up front, before any downstream code reads them -- self.args.covgrp/comparegrp1/
        comparegrp2/heatgrp/volgrp are each read directly (and more than once: the log2FC
        precompute above, cmap_dict, and the per-plot-type calls in plot()), so resolving them
        here keeps every call site consistent instead of only patching the plots*.py modules
        that happen to validate their own grp argument.
        '''
        for param_name in ('covgrp', 'comparegrp1', 'comparegrp2', 'heatgrp', 'volgrp'):
            setattr(self.args, param_name, toolsTG.resolve_grp_column(self.adata, getattr(self.args, param_name), param_name))

    def main(self):
        # Generate graphs
        if self.args.verbose:
            print('Generating graphs with the following parameters:\n')
            for i in self.args.__dict__: print(f'{i}: {self.args.__dict__[i]}')
            print('')
        # Remove coverage from self.args.graphtypes and add it to non_pooled_graphs
        non_pooled_graphs = []
        # if 'bar' in self.args.graphtypes:
        #     self.args.graphtypes.remove('bar')
        #     if 'count' not in self.args.graphtypes:
        #         self.args.graphtypes.append('count')
        if 'coverage' in self.args.graphtypes:
            self.args.graphtypes.remove('coverage')
            non_pooled_graphs.append('coverage')
        # Pool if applicable
        if self.args.threads > 1 and len(self.args.graphtypes) > 1:
            print(f'Multithreading enabled with {self.args.threads} threads to generate plots\n')
            # Create a multiprocessing pool
            pool = multiprocessing.Pool(self.args.threads)
            # Generate graphs
            pool_output = pool.starmap(self.plot, [(gt, True) for gt in self.args.graphtypes])
            # Close the pool and wait for the work to finish
            pool.close()
            pool.join()
            for gt in non_pooled_graphs:
                self.plot(gt)
            for po in pool_output:
                if po:
                    print(po + '\n')
        else:
            # Combine non_pooled_graphs with self.args.graphtypes
            self.args.graphtypes += non_pooled_graphs
            for gt in self.args.graphtypes:
                self.plot(gt)

    def plot(self, gt, threaded=None):
        if threaded:
            threaded = f'Generating {gt} plots...\n'
        else:
            print(f'Generating {gt} plots...')
        adata_c = self.adata.copy()
        # Define the colormap to use for the graph type
        colormap = None
        cmapgrp = gt
        if gt == 'coverage':
            if not self.args.covgrp:
                cmapgrp = 'group'
        cmappar = self.cmap_dict.get(cmapgrp, None)
        if self.args.colormap:
            if cmappar in self.args.colormap:
                colormap = self.args.colormap[cmappar]
        # Create the output directory (namespaced by --variant when non-default, so different
        # --variant runs into the same --output don't overwrite each other's files)
        output = self.args.output + '/' + gt + '/'
        if self.variant_spec.raw != 'norm:full':
            output = self.args.output + '/' + gt + '/' + self.variant_spec.raw.replace(':', '_') + '/'
        if threaded:
            threaded += toolsTG.builder(output) + '\n'
        else:
            print(toolsTG.builder(output))
        # Plot specific parameters
        if gt == 'cluster':
            if not 'cluster_runinfo' in self.adata.uns:
                if threaded:
                    threaded += 'No cluster run information found in AnnData object. Please run the cluster command first.\n'
                else:
                    print('No cluster run information found in AnnData object. Please run the cluster command first.\n')
            else:
                threaded = plotsCluster.visualizer(adata_c, self.args.clustergrp, self.args.clusteroverview, self.args.clusternumeric, self.args.clusterlabels, self.args.clustermask, colormap, output, threaded=threaded).generate_plots()
        if gt == 'compare':
            threaded = plotsCompare.visualizer(adata_c, self.args.comparegrp1, self.args.comparegrp2, colormap, output, threaded=threaded)
        if gt == 'correlation':
            threaded = plotsCorrelation.visualizer(adata_c, self.args.corrmethod, self.args.corrgroup, output, threaded=threaded)
        if gt == 'count':
            threaded = plotsCount.visualizer(adata_c, self.args.colormap if self.args.colormap else {}, output, threaded=threaded)
        if gt == 'coverage':
            pcV = plotsCoverage.visualizer(adata_c, self.args.threads, self.args.covgrp, self.args.covobs, self.args.covtype, self.args.covgap, self.args.covmethod, colormap, output)
            # Generate folders/subfolders if coveragecombine is specified
            print(toolsTG.builder(f'{output}{self.args.covobs}/'))
            print(toolsTG.builder(f'{output}{self.args.covobs}/low_coverage/'))
            # Generate coverage plots with combine or split pdfs
            if not self.args.combinedpdfonly:
                print('Generating individual coverage plots pdfs...')
                pcV.generate_split()
            print('Generating combined coverage plots pdf...')
            pcV.generate_combine()
        if gt == 'heatmap':
            threaded = plotsHeatmap.visualizer(adata_c, self.args.heatgrp, self.args.diffrts, self.args.heatcutoff, self.args.heatbound, self.args.heatsubplots, output, threaded=threaded, config_name=self.config_name, overwrite=self.args.regen_uns)
        if gt == 'logo':
            plotsSeqlogo.visualizer(adata_c, self.args.logogrp, self.args.logomanualgrp, self.args.logomanualname, self.args.logopseudocount, self.args.logosize, self.args.ccatail, self.args.pseudogenes, self.args.logornamode, output).generate_plots()
        if gt == 'pca':
            threaded = plotsPca.visualizer(adata_c, self.args.pcamarkers, self.args.pcacolors, self.args.pcareadtypes, colormap, output, threaded=threaded)
        if gt == 'radar':
            if 'all' in self.args.radarmethod:
                self.args.radarmethod = ['mean', 'median', 'max', 'sum']
            for radarmethod in self.args.radarmethod:
                pRd = plotsRadar.visualizer(adata_c, self.args.radargrp, radarmethod, self.args.radarscaled, colormap, output, threaded=threaded)
                threaded = pRd.isotype_plots()
        if gt == 'volcano':
            threaded = plotsVolcano.visualizer(adata_c, self.args.volgrp, self.args.diffrts, self.args.volcutoff, output, colormap=colormap, toplabels=self.args.vollabels, threaded=threaded, config_name=self.config_name, overwrite=self.args.regen_uns)
        # Return threaded output  
        if threaded:
            threaded += f'{gt.capitalize()} plots generated!\n'
            return threaded
        else:
            print(f'{gt.capitalize()} plots generated!\n')
