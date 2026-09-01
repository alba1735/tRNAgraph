import anndata as ad
import itertools
import os
import multiprocessing
import json
import logging
from pydantic import ValidationError
from . import toolsTG
from .toolsSchemas import GraphFilterConfig, OUTPUT_FORMATS
from .lazy_imports import (
    plotsCount, plotsCluster, plotsCompare, plotsCorrelation,
    plotsCoverage, plotsHeatmap, plotsMismatch, plotsSeqlogo, plotsPca, plotsRadar, plotsVolcano
)

#: Every `graph` option whose value names something inside the AnnData object, with the
#: vocabulary it is drawn from. Kept as one table so a newly added option is validated by
#: being listed here rather than by remembering to call a validator at its use site.
GRAPH_LABEL_PARAMS = (
    ('clustergrp', 'obs'),
    ('clusterlabels', 'obs'),
    ('comparegrp1', 'obs'),
    ('comparegrp2', 'obs'),
    ('corrgroup', 'obs'),
    ('covgrp', 'obs'),
    ('covobs', 'obs'),
    ('heatgrp', 'obs'),
    ('logogrp', 'obs'),
    ('pcamarkers', 'obs'),
    ('pcacolors', 'obs'),
    ('radargrp', 'obs'),
    ('volgrp', 'obs'),
    ('covtype', 'coverage'),
    ('diffrts', 'readtype'),
    ('pcareadtypes', 'readtype'),
)


#: Graph types whose FILENAMES already record the read basis, so a basis directory segment
#: would only duplicate identical output into a second place. Coverage's filenames carry
#: --covtype's category; PCA's carry the readtype label, which differs by basis.
BASIS_IN_FILENAME_TYPES = frozenset({'coverage', 'pca'})


class anndataGrapher:
    '''
    Class to generate graphs from an AnnData object by calling the appropriate graphing functions
    '''
    def __init__(self, args):
        self.logger = logging.getLogger(__name__)
        self.args = args
        self.quiet = getattr(args, 'quiet', False)
        # Read and validate --config FIRST, and apply its `flags` block, because nearly
        # everything below is derived from options that block can set: the variant view, the
        # read basis, covtype, the resolved grouping columns, and the graphtypes expansion.
        # Only the obs/var filtering is deferred to its original position further down, since
        # that needs the variant view this parse precedes.
        self.config = self._load_config()
        self._apply_config_flags()
        self.adata_original = ad.read_h5ad(self.args.anndata)
        # Resolve the requested normalization:split-tag ONCE, into a working copy, so every
        # downstream plot module reads .X/.obs[...]/.uns[...] exactly as it does for the
        # full/default variant -- see toolsTG.build_variant_view() for why this resolved copy
        # must never be written back to self.args.anndata directly.
        self.variant_spec = toolsTG.parse_variant(self.adata_original, getattr(self.args, 'variant', 'norm:full'))
        self.adata = toolsTG.build_variant_view(self.adata_original, self.variant_spec)
        # Resolve the read basis ONCE, here, for the whole command -- graphs plot unique
        # (transcript-specific) counts unless --allreads says otherwise. This must happen
        # before _precompute_and_persist_log2fc() below, which reads self.args.diffrts to
        # decide which log2FC entries to compute and cache: resolving later would populate
        # the cache against the wrong basis and every heatmap/volcano would silently read it.
        self.read_basis = toolsTG.read_basis(getattr(self.args, 'allreads', False))
        self.args.covtype = toolsTG.resolve_covtype(getattr(self.args, 'covtype', None), self.read_basis)
        self.config_name = 'default'
        # Resolve grouping-column args before anything below reads them -- the log2FC
        # precompute a few lines down reads self.args.heatgrp/volgrp directly, ahead of the
        # plots*.py modules that separately validate their own grp argument.
        self._validate_label_args()
        # Load cmap dict for each graph type
        self.cmap_dict = {'cluster':self.args.clustergrp, 'compare':self.args.comparegrp1, \
                          'coverage':self.args.covgrp, 'pca':self.args.pcacolors, 'radar':self.args.radargrp, \
                          'volcano':self.args.volgrp}
        # Load all graph types if specified
        if self.args.graphtypes == 'all' or 'all' in self.args.graphtypes:
            self.args.graphtypes = ['cluster', 'correlation', 'count', 'coverage', 'heatmap', 'logo', 'mismatch', 'pca', 'radar', 'volcano']
            self.args.clusteroverview = True
        # Load max threads available unless specified
        if self.args.threads == 0:
            try:
                # This is a linux only function but is less likely to cause problems than multiprocessing.cpu_count()
                self.args.threads = len(os.sched_getaffinity(0))
            except:
                self.args.threads = multiprocessing.cpu_count()
        # Apply the config's obs/var filters. The file was already read and its `flags` block
        # applied at the top of __init__; only the filtering waits until here, because it
        # needs the variant view built above.
        if self.config is not None:
            config = self.config
            self.args.output += '/' + config.name
            self.config_name = config.name
            self.logger.info(toolsTG.builder(self.args.output))
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
                    self.logger.info('Filtering AnnData object by observation: ' + k + ' , ' + str(v))
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
                    self.logger.info('Filtering AnnData object by variable: ' + k + ' , ' + str(v))
                    self.adata = self.adata[:, self.adata.var[k].isin(v)]
            self.logger.info('Config file loaded.\n')

        else:
            self.args.config = {}
        # Load the style file if specified. It carries both the palette and the presentation
        # settings, so `self.style` is consulted for per-graph-type figure settings while
        # `self.args.colormap` keeps its existing meaning (the colors block) for every plot
        # module that already takes a colormap argument.
        fmt = getattr(self.args, 'format', None)
        if fmt is not None and fmt not in OUTPUT_FORMATS:
            raise ValueError(f"--format must be one of {list(OUTPUT_FORMATS)}, got '{fmt}'.")
        self.style = None
        self.args.colormap = None
        if getattr(self.args, 'style', None):
            self.style = toolsTG.load_style_file(self.args.style, self.logger)
            self.args.colormap = self.style.colors or None
            self.logger.info('Style loaded.\n')
        # Check for heatmap or volcano in graph types and if present check for readcount_cutoff in log2FC - Will precompute this and save it back to uns
        # This is done now to prevent saving issues later if multiprocessing is used
        self._precompute_and_persist_log2fc()

    def _load_config(self):
        '''Read and validate --config, or None when none was given.'''
        if not getattr(self.args, 'config', None):
            return None
        self.logger.info('Loading config file: ' + self.args.config)
        with open(self.args.config, 'r') as f:
            raw_config = json.load(f)
        try:
            return GraphFilterConfig.model_validate(raw_config)
        except ValidationError as e:
            raise ValueError(f'Invalid config file {self.args.config}:\n{e}') from e

    def _apply_config_flags(self):
        '''
        Lay the config's `flags` block over the args, skipping anything typed on the command
        line so a CLI flag always wins -- the same precedence --style already promises.

        `cli_specified` is a frozenset of parameter names built in cli.py from click's
        ParameterSource. It is absent when a namespace is constructed directly (the Python
        API, and tests), in which case nothing is treated as explicitly typed and the file
        applies in full.
        '''
        if self.config is None or self.config.flags is None:
            return
        typed = getattr(self.args, 'cli_specified', None) or frozenset()
        for key, value in self.config.flags.model_dump(exclude_none=True).items():
            if key in typed:
                self.logger.info(
                    f'Config sets {key}, but --{key} was given on the command line; '
                    f'keeping the command-line value.')
                continue
            setattr(self.args, key, value)
            self.logger.info(f'Config sets {key} = {value!r}')

    def _precompute_and_persist_log2fc(self):
        '''
        Precompute log2FC for every {heatgrp,volgrp} x diffrts x {heatcutoff,volcutoff} combo
        (plus volcano's fixed overview readtypes) the CURRENT invocation's flags actually need,
        and persist anything newly computed back to the original h5ad -- so a later `graph` run
        with matching flags (e.g. regenerating figures) hits the cache instead of recomputing.
        Runs once, before the graph-type multiprocessing.Pool is created, so it's safe to use
        real parallelism here (n_cpus=self.args.threads) unlike the same-shaped on-demand calls
        inside plotsVolcano.py/plotsHeatmap.py, which run inside pooled workers and must stay at
        adataLog2FC's safe default of n_cpus=1 to avoid a nested-process-pool deadlock.

        "Was anything actually new" is tracked via adataLog2FC.main()'s own computed_fresh flag,
        not by diffing the uns['log2FC'] dict before/after -- a before/after comparison is
        fundamentally the wrong tool here: a shallow dict.copy() shares the same nested dict
        objects between the "before"/"after" snapshots, so adding a new readtype/cutoff under an
        already-existing config_name/compare path (the common case, 'default'/'group' by default)
        never registers as a difference (equality short-circuits on object identity) -- and using
        a deep copy instead just trades that silent bug for a crash (`ValueError: truth value of
        a DataFrame is ambiguous`), since dict equality then actually has to compare the nested
        DataFrame values directly.
        '''
        threads = getattr(self.args, 'threads', None) or None
        any_computed = False
        if 'heatmap' in self.args.graphtypes or 'volcano' in self.args.graphtypes:
            for grp in list(set([self.args.heatgrp, self.args.volgrp])):
                for readtype in self.resolved_diffrts():
                    for cutoff in list(set([self.args.heatcutoff, self.args.volcutoff])):
                        log2fc = toolsTG.adataLog2FC(self.adata, grp, readtype, readcount_cutoff=cutoff, config_name=self.config_name, overwrite=self.args.regen_uns, n_cpus=threads, lfc_shrink=getattr(self.args, 'lfcshrink', True))
                        log2fc.main()
                        any_computed = any_computed or log2fc.computed_fresh
        if 'volcano' in self.args.graphtypes:
            # The volcano combined overview page always uses these two read types (mirroring
            # PCA's default --pcareadtypes), regardless of what --diffrts requests.
            for readtype in plotsVolcano.OVERVIEW_TRNA_READTYPES:
                log2fc = toolsTG.adataLog2FC(self.adata, self.args.volgrp, readtype, readcount_cutoff=self.args.volcutoff, config_name=self.config_name, overwrite=self.args.regen_uns, n_cpus=threads, lfc_shrink=getattr(self.args, 'lfcshrink', True))
                log2fc.main()
                any_computed = any_computed or log2fc.computed_fresh
        if any_computed or self.args.regen_uns:
            self.logger.info('The log2FC uns dictionary has been updated.\n')
            # Persist onto the ORIGINAL (unresolved) adata, into the correct namespaced
            # location -- never write self.adata (the resolved view) back to disk, since for a
            # split variant it would overwrite the real full/default variant's data.
            if self.variant_spec.tag == 'full':
                self.adata_original.uns['log2FC'] = self.adata.uns['log2FC']
            else:
                self.adata_original.uns.setdefault('size_splits', {}).setdefault(self.variant_spec.tag, {})['log2FC'] = self.adata.uns['log2FC']
            self.adata_original.write(self.args.anndata)

    def resolved_diffrts(self):
        '''
        --diffrts as obs column names for the resolved read basis. Kept as a method rather
        than resolved in place so the bare readtype survives in self.args for the output
        path and for --verbose, and so heatmap/volcano and the log2FC precompute cannot
        drift apart over which column a given readtype means.
        '''
        return [toolsTG.resolve_readtype(rt, self.read_basis, self.adata) for rt in self.args.diffrts]

    def _validate_label_args(self):
        '''
        Check every label-valued option against the resolved AnnData object once, up front,
        before any downstream code reads them.

        This replaces a fallback that substituted 'sample' for an unrecognised column and
        carried on: a typo'd --covgrp then produced a complete, plausible set of figures
        grouped by the wrong thing, which is worse than an error because nothing about the
        output says so. It also covers every option rather than the five that had a
        fallback -- --clustergrp, --clusterlabels, --corrgroup, --logogrp, --pcamarkers,
        --pcacolors and --radargrp were validated nowhere and failed at first use inside
        pandas, with a message naming neither the flag nor the alternatives.

        Batched deliberately: a graph command carries a dozen of these, so reporting them one
        at a time turns fixing a command line into as many round trips as there are typos.
        '''
        requests = []
        for param_name, domain in GRAPH_LABEL_PARAMS:
            value = getattr(self.args, param_name, None)
            # An unset optional (e.g. --clusterlabels) means "no such dimension", not a
            # column named None.
            if value is None:
                continue
            values = value if isinstance(value, (list, tuple)) else [value]
            requests.extend((param_name, v, domain) for v in values if v is not None)
        toolsTG.validate_labels(self.adata, requests, extra_problems=self._parameter_problems())

    def _parameter_problems(self):
        '''
        Combinations that cannot work even though every label in them exists.

        Reported alongside the unknown labels rather than raised separately, so one pass finds
        everything wrong with the command.
        '''
        problems = []
        # Gated on compare actually being requested, because --comparegrp1 and --comparegrp2
        # BOTH DEFAULT to 'group': an unconditional check would abort every ordinary run of a
        # command that never asked for a compare plot. compare is excluded from `-g all` by
        # design, and this runs before the 'all' expansion, so it appears in graphtypes only
        # when named explicitly -- which is exactly when the collision matters.
        graphtypes = getattr(self.args, 'graphtypes', []) or []
        grp1 = getattr(self.args, 'comparegrp1', None)
        if 'compare' in graphtypes and grp1 is not None and grp1 == getattr(self.args, 'comparegrp2', None):
            # The compare plot takes a fold change BETWEEN --comparegrp2 values within each
            # --comparegrp1 value, so one column for both is not a comparison. Left
            # unchecked, log2fc_compare_df pivots on a duplicated column and pandas raises
            # "Grouper for 'group' not 1-dimensional", which names neither flag.
            problems.append(
                f"--comparegrp1 and --comparegrp2 both name '{grp1}'; they must name different "
                f"columns, since the fold change is taken between --comparegrp2 values within "
                f"each --comparegrp1 value"
            )
        corrmask = getattr(self.args, 'corrmask', None)
        if corrmask is not None and corrmask not in plotsCorrelation.CORR_MASK_CHOICES:
            # Checked here rather than where the matrix is drawn: correlation is one of ten
            # graph types, so a run can spend most of its time before reaching it, and a typo
            # should cost nothing.
            problems.append(
                f"--corrmask '{corrmask}' is not one of "
                f"{', '.join(plotsCorrelation.CORR_MASK_CHOICES)}"
            )
        return problems

    def _compute_graph_weight(self, gt):
        '''
        Best-effort upfront estimate of how many plots/pages graph type `gt` will produce, used
        as its toolsTG.PhaseTracker weight so the outer graphing bar reflects real proportions --
        coverage can produce hundreds/thousands of plots (one per --covobs value), most other
        types only a handful -- rather than equal-weighting every type the same. A rough
        estimator, not a source of truth: wrapped in try/except so a weight-computation quirk (an
        unexpected column, an edge-case config) never blocks the actual graphing run, falling
        back to a modest constant instead.
        '''
        try:
            if gt == 'coverage':
                return max(1, int(self.adata.obs[self.args.covobs].nunique()))
            if gt == 'pca':
                # One call per --pcareadtypes entry, plus the non-tRNA and combined-RNA overview
                # plots (each producing its own evr/scatter/pairplot trio).
                return max(1, len(self.args.pcareadtypes) + 2)
            if gt == 'radar':
                methods = self.args.radarmethod if isinstance(self.args.radarmethod, list) else [self.args.radarmethod]
                aminos = int(self.adata.obs['amino'].nunique()) if 'amino' in self.adata.obs.columns else 20
                return max(1, len(methods) * 2 * aminos)
            if gt == 'volcano':
                pairs = max(1, len(list(itertools.combinations(self.adata.obs[self.args.volgrp].dropna().unique(), 2))))
                return max(1, len(self.args.diffrts) * pairs)
            if gt == 'heatmap':
                pairs = max(1, len(list(itertools.combinations(self.adata.obs[self.args.heatgrp].dropna().unique(), 2))))
                return max(1, len(self.args.diffrts) * pairs)
            if gt == 'cluster':
                if self.args.clustergrp in self.adata.obs.columns:
                    return max(1, int(self.adata.obs[self.args.clustergrp].nunique()))
                return 3
            if gt == 'compare':
                if self.args.comparegrp2 in self.adata.obs.columns:
                    return max(1, 2 * int(self.adata.obs[self.args.comparegrp2].nunique()))
                return 2
            if gt == 'correlation':
                return max(1, len([c for c in self.adata.obs.columns if '_norm' in c]))
            if gt == 'count':
                return 2
            if gt == 'mismatch':
                # Two overview pages, the read-level histogram, and one page per amino acid.
                aminos = int(self.adata.obs['amino'].nunique()) if 'amino' in self.adata.obs.columns else 20
                return max(1, 3 + aminos)
            if gt == 'logo':
                if self.args.logogrp in self.adata.obs.columns:
                    return max(1, int(self.adata.obs[self.args.logogrp].nunique()))
                return 4
        except Exception:
            pass
        return 1

    def _output_dir_for(self, gt):
        '''
        Where one graph type's files go.

        The path names CONTENT, not the flags that were typed: a segment is added only when
        nothing already in the path distinguishes that output. `--variant` always qualifies,
        since the same plot of a different normalization is a different picture and no
        filename records which one it was.

        The read basis qualifies every graph type in BASIS_IN_FILENAME_TYPES, which are the
        ones whose filenames already name it. Coverage is there because `--covtype` names the
        category. PCA is there because every filename carries plotsPca._readtype_label()
        ('total' against 'total_unique') -- and because plotsPca pins the both-bases
        comparison in OVERVIEW_TRNA_READTYPES so it survives whatever --pcareadtypes asks
        for, which meant a default run and an --allreads run emitted exactly the same two
        files into two directories.
        '''
        output = self.args.output + '/' + gt + '/'
        if self.variant_spec.raw != 'norm:full':
            output += self.variant_spec.raw.replace(':', '_') + '/'
        if gt not in BASIS_IN_FILENAME_TYPES and self.read_basis != toolsTG.READ_BASIS_UNIQUE:
            output += self.read_basis + '/'
        return output

    def _plot_wrapper(self, args):
        gt, threaded = args
        return self.plot(gt, threaded)

    def main(self):
        # Generate graphs
        if self.args.verbose:
            self.logger.info('Generating graphs with the following parameters:\n')
            for i in self.args.__dict__: self.logger.info(f'{i}: {self.args.__dict__[i]}')
            self.logger.info('')
        # Remove coverage from self.args.graphtypes and add it to non_pooled_graphs
        non_pooled_graphs = []
        # if 'bar' in self.args.graphtypes:
        #     self.args.graphtypes.remove('bar')
        #     if 'count' not in self.args.graphtypes:
        #         self.args.graphtypes.append('count')
        if 'coverage' in self.args.graphtypes:
            self.args.graphtypes.remove('coverage')
            non_pooled_graphs.append('coverage')
        # One shared outer tracker spans every graph type -- weighted by _compute_graph_weight()
        # so coverage's hundreds/thousands of plots are proportionally represented next to the
        # other types' handful each.
        all_graphtypes = self.args.graphtypes + non_pooled_graphs
        self.phase_tracker = toolsTG.PhaseTracker(
            phases=all_graphtypes, logger=self.logger, desc="Graphing",
            weights=[self._compute_graph_weight(gt) for gt in all_graphtypes],
        )
        # Pool if applicable
        if self.args.threads > 1 and len(self.args.graphtypes) > 1:
            self.logger.info(f'Multithreading enabled with {self.args.threads} threads to generate plots\n')
            # Create a multiprocessing pool
            pool = multiprocessing.Pool(self.args.threads)
            # imap() (ORDERED, not imap_unordered) -- each completed type's `with phase(): pass`
            # must register against the SAME position in `all_graphtypes` that the tracker was
            # built from, so the phase-level message names the right graph type.
            pool_output = []
            for po in pool.imap(self._plot_wrapper, [(gt, True) for gt in self.args.graphtypes]):
                with self.phase_tracker.phase():
                    pass
                pool_output.append(po)
            # Close the pool and wait for the work to finish
            pool.close()
            pool.join()
            for gt in non_pooled_graphs:
                with self.phase_tracker.phase():
                    self.plot(gt)
            for po in pool_output:
                if po:
                    self.logger.info(po + '\n')
        else:
            # Combine non_pooled_graphs with self.args.graphtypes
            self.args.graphtypes += non_pooled_graphs
            for gt in self.args.graphtypes:
                with self.phase_tracker.phase():
                    self.plot(gt)

    def style_for(self, gt):
        '''Resolved presentation settings for one graph type, CLI flags already applied.'''
        return toolsTG.resolve_plot_style(getattr(self, 'style', None), gt,
                                          format=getattr(self.args, 'format', None))

    def plot(self, gt, threaded=None):
        # rcParam-expressible style settings (font size, default figure size) have to be in
        # force while the figures are CREATED, not when they are saved, so the whole dispatch
        # runs inside the context rather than each module handling it.
        with toolsTG.style_context(self.style_for(gt)):
            return self.dispatch_plot(gt, threaded)

    def dispatch_plot(self, gt, threaded=None):
        if threaded:
            threaded = f'Generating {gt} plots...\n'
        else:
            self.logger.info(f'Generating {gt} plots...')
        adata_c = self.adata.copy()
        settings = self.style_for(gt)
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
        output = self._output_dir_for(gt)
        if threaded:
            threaded += toolsTG.builder(output) + '\n'
        else:
            self.logger.info(toolsTG.builder(output))
        # Plot specific parameters
        if gt == 'cluster':
            if not 'cluster_runinfo' in self.adata.uns:
                if threaded:
                    threaded += 'No cluster run information found in AnnData object. Please run the cluster command first.\n'
                else:
                    self.logger.warning('No cluster run information found in AnnData object. Please run the cluster command first.\n')
            else:
                threaded = plotsCluster.visualizer(adata_c, self.args.clustergrp, self.args.clusteroverview, self.args.clusternumeric, self.args.clusterlabels, self.args.clustermask, colormap, output, threaded=threaded, read_basis=self.read_basis, settings=settings).generate_plots()
        if gt == 'compare':
            threaded = plotsCompare.visualizer(adata_c, self.args.comparegrp1, self.args.comparegrp2, colormap, output, threaded=threaded, read_basis=self.read_basis, settings=settings)
        if gt == 'correlation':
            threaded = plotsCorrelation.visualizer(adata_c, self.args.corrmethod, self.args.corrgroup, output, threaded=threaded, is_full_variant=self.variant_spec.tag == 'full', read_basis=self.read_basis, settings=settings, corr_mask=getattr(self.args, 'corrmask', None))
        if gt == 'count':
            threaded = plotsCount.visualizer(adata_c, self.args.colormap if self.args.colormap else {}, output, threaded=threaded, settings=settings)
        if gt == 'coverage':
            pcV = plotsCoverage.visualizer(adata_c, self.args.threads, self.args.covgrp, self.args.covobs, self.args.covtype, self.args.covgap, self.args.covmethod, colormap, output, phase_tracker=self.phase_tracker, quiet=self.quiet, settings=settings)
            # plotsCoverage owns its own layout below `output`: per-category plots under
            # <covtype-alias>/, the specificity overview at the base.
            pcV.build_output_dirs()
            # Generate coverage plots with combine or split pdfs
            if not self.args.combinedpdfonly:
                self.logger.info('Generating individual coverage plots pdfs...')
                pcV.generate_split()
            self.logger.info('Generating combined coverage plots pdf...')
            pcV.generate_combine()
            self.logger.info('Generating coverage specificity overview pdf...')
            pcV.generate_partition_overview()
        if gt == 'heatmap':
            threaded = plotsHeatmap.visualizer(adata_c, self.args.heatgrp, self.resolved_diffrts(), self.args.heatcutoff, self.args.heatbound, self.args.heatsubplots, output, threaded=threaded, config_name=self.config_name, overwrite=self.args.regen_uns, settings=settings)
        if gt == 'logo':
            plotsSeqlogo.visualizer(adata_c, self.args.logogrp, self.args.logomanualgrp, self.args.logomanualname, self.args.logopseudocount, self.args.logosize, self.args.ccatail, self.args.pseudogenes, self.args.logornamode, output, read_basis=self.read_basis, settings=settings).generate_plots()
        if gt == 'mismatch':
            threaded = plotsMismatch.visualizer(adata_c, self.args.colormap if self.args.colormap else {}, output, self.args.mismatchpseudocount, threaded=threaded, settings=settings).generate_plots()
        if gt == 'pca':
            threaded = plotsPca.visualizer(adata_c, self.args.pcamarkers, self.args.pcacolors, self.args.pcareadtypes, colormap, output, threaded=threaded, is_full_variant=self.variant_spec.tag == 'full', read_basis=self.read_basis, settings=settings)
        if gt == 'radar':
            if 'all' in self.args.radarmethod:
                self.args.radarmethod = ['mean', 'median', 'max', 'sum']
            for radarmethod in self.args.radarmethod:
                pRd = plotsRadar.visualizer(adata_c, self.args.radargrp, radarmethod, self.args.radarscaled, colormap, output, threaded=threaded, read_basis=self.read_basis, settings=settings)
                threaded = pRd.isotype_plots()
        if gt == 'volcano':
            threaded = plotsVolcano.visualizer(adata_c, self.args.volgrp, self.resolved_diffrts(), self.args.volcutoff, output, colormap=colormap, toplabels=self.args.vollabels, threaded=threaded, config_name=self.config_name, overwrite=self.args.regen_uns, is_full_variant=self.variant_spec.tag == 'full', xlim=self.args.volxlim, settings=settings)
        # Return threaded output  
        if threaded:
            threaded += f'{gt.capitalize()} plots generated!\n'
            return threaded
        else:
            self.logger.info(f'{gt.capitalize()} plots generated!\n')
