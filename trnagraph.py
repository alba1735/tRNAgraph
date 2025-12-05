#!/usr/bin/env python3

import anndata as ad
import os
import sys
import json
import contextlib
import typer
from typing import Optional, List
from types import SimpleNamespace
# Custom functions
try:
    from . import toolsTestSuite
    from .lazy_imports import (
        toolsMap, toolsTDatabase, toolsTrim, toolsTG
    )
    from .adataGraph import anndataGrapher
    from .adataMerge import anndataMerger
    from .adataCluster import anndataCluster
    from .adataBuild import trax2anndata
except ImportError:
    # Fallback for script execution: add parent directory to path and import as package
    # This ensures that relative imports inside sub-modules (like adataGraph) work correctly
    import sys
    import os
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    
    from tRNAgraph import toolsTestSuite
    from tRNAgraph.lazy_imports import (
        toolsMap, toolsTDatabase, toolsTrim, toolsTG
    )
    from tRNAgraph.adataGraph import anndataGrapher
    from tRNAgraph.adataMerge import anndataMerger
    from tRNAgraph.adataCluster import anndataCluster
    from tRNAgraph.adataBuild import trax2anndata

def _main_logic(args):
    '''
    Main function for running the logic and calling the appropriate class
    '''
    if args.mode == 'makedb':
        if not os.path.isfile(args.genome):
            raise Exception('Error: genome fasta file does not exist.')
        print('Building tRNA database...')
        toolsTDatabase.tRNADatabaseBuilder(args).main()
        print('Done!\n')
    elif args.mode == 'trim':
        # Check for fastp
        import shutil
        if shutil.which('fastp') is None:
            raise Exception("Error: 'fastp' is not installed or not in PATH. Please install it (e.g., 'conda install -c bioconda fastp').")
        # Validate manifest existence
        if not os.path.isfile(args.manifest):
            raise Exception(f'Error: Manifest file does not exist: {args.manifest}')
        print('Starting fastp trimming pipeline...')
        toolsTrim.FastpTrimmer(args).process()
        print('Done!\n')
    elif args.mode == 'map':
        print('Mapping samples...')
        toolsMap.MapSamples(args).main()
        print('Done!\n')
    # Read database object or create one from trax run if none provided
    elif args.mode == 'build':
        # If experiment is provided, we are running analysis
        if hasattr(args, 'experiment') and args.experiment:
            if not args.database or not args.samples:
                raise Exception('Error: database and samples are required for analysis.')
            args.traxdir = args.experiment # Output of analysis is in experiment directory
        elif not args.traxdir:
             raise Exception('Error: traxdir is required if not running analysis.')

        # Clean the path to the trax directory
        args.traxdir = os.path.abspath(args.traxdir)
        
        # Raise exception if metadata file is empty or doesn't exist
        if args.metadata:
            if not os.path.isfile(args.metadata):
                raise Exception('Error: metadata file does not exist.')
        print('Building AnnData object...')
        # Create path to output directory if it doesn't exist
        print(toolsTG.builder(args.output))
        # Create AnnData object
        analysis_args = args if hasattr(args, 'experiment') and args.experiment else None
        trax2anndata(args.traxdir, args.metadata, args.output, analysis_args).create()
        print('Done!\n')
    elif args.mode == 'merge':
        # Raise exception if h5ad file is empty or doesn't exist
        if not os.path.isfile(args.anndata1):
            raise Exception('Error: first h5ad file does not exist.')
        if not os.path.isfile(args.anndata2):
            raise Exception('Error: second h5ad file does not exist.')
        # Create path to output directory if it doesn't exist
        print(toolsTG.builder(args.output))
        # Merge AnnData objects
        print('Merging database objects...\n')
        anndataMerger(args).merge()
        print('Done!\n')
    elif args.mode == 'cluster':
        # Raise exception if h5ad file is empty or doesn't exist
        if not os.path.isfile(args.anndata):
            raise Exception('Error: h5ad file does not exist.')
        # Create path to output directory if it doesn't exist
        print(toolsTG.builder(args.output))
        print('Clustering data from database object...\n')
        anndataCluster(args).main()
        print('Done!\n')
    elif args.mode == 'graph':
        # Create output directory if it doesn't exist
        args.output = os.path.abspath(args.output)
        print(toolsTG.builder(args.output))
        # Raise exception if h5ad file is empty or doesn't exist
        if not os.path.isfile(args.anndata):
            raise Exception('Error: h5ad file does not exist.')
        print('Graphing data from database object...\n')
        anndataGrapher(args).main()
        print('Done!\n')
    elif args.mode == 'log2fc':
        # Raise exception if h5ad file is empty or doesn't exist
        if not os.path.isfile(args.anndata):
            raise Exception('Error: h5ad file does not exist.')
        # Load the AnnData object
        adata = ad.read_h5ad(args.anndata)
        # Load config file for name if specified
        config_name = 'default'
        if args.config:
            with open(args.config, 'r') as f:
                args.config = json.load(f)
            if 'name' in args.config:
                # self.args.output += '/' + self.args.config['name']
                config_name = args.config['name']
                # print(toolsTG.builder(self.args.output))
            else:
                raise ValueError('Config file must contain a "name" field')
        print('Calculating log2FC for database object...\n')
        adata_copy = adata.copy()
        log2FC_dict = adata.uns['log2FC']
        for readtype in [f'nreads_{i}_norm' for i in args.readtypes]:
            for cutoff in args.cutoff:
                toolsTG.adataLog2FC(adata_copy, args.group, readtype, readcount_cutoff=cutoff, config_name=config_name, overwrite=True).main()
        # if log2FC_dict.items() != adata_copy.uns['log2FC'].items(): # Can fix this to be more efficient later
        print('The log2FC uns dictionary has been updated.\nWriting h5ad database object to: ' + args.anndata)
        adata_copy.write(args.anndata)
        # else:
        # print('The log2FC uns dictionary has not been updated.\n')
        print('Done!\n')
    elif args.mode == 'csv':
        args.output = os.path.abspath(args.output)
        # Add the name of the h5ad file to the output directory minus the extension account for periods in the name removing just .h5ad
        args.output += '/' + '.'.join(os.path.basename(args.anndata).split('.')[:-1]) + '/'
        print(toolsTG.builder(args.output))
        adata = ad.read_h5ad(args.anndata)
        print('Writing csv files to: ' + args.output)
        adata.write_csvs(args.output, skip_data=False)
        print('Done!\n')
    elif args.mode == 'test':
        toolsTestSuite.demoPipeline(args).main()
        print('Done!\n')



app = typer.Typer(
    help="tRNAgraph is a tool for analyzing tRNA-seq data generated from tRAX. It can be used to create an AnnData database object from \
            a trax output folder, or to analyze an existing database object and generate expanded visulizations. The database object can also be used to \
            perform further multivariate analysis such as clustering and classification of readcoverages.",
    add_completion=False,
    no_args_is_help=True
)

preprocess_app = typer.Typer(help="Preprocess raw fastq/fasta files for tRNA analysis", no_args_is_help=True)
app.add_typer(preprocess_app, name="preprocess")

tools_app = typer.Typer(help="Extra utilities for working with tRNAgraph objects", no_args_is_help=True)
app.add_typer(tools_app, name="tools")

def run_logic(args):
    # Set log file if specified
    if args.log:
        sys.stdout = open(args.log, 'w')
    # Run main function
    if args.quiet:
        with open(os.devnull, 'w') as f, contextlib.redirect_stdout(f):
            _main_logic(args)
    else:
        _main_logic(args)

@preprocess_app.command("makedb", help="Build bowtie2 index from gtRNAdb/tRNAScan-SE output and reference genome")
def makedb(
    genome: str = typer.Option(..., "-g", "--genome", help="Specify location of the reference genome fasta file"),
    trnaout: str = typer.Option(..., "-t", "--trnaout", help="Specify location of the tRNAScan-SE out file"),
    trnafa: str = typer.Option(..., "-r", "--trnafa", help="Specify location of the tRNA reference fasta file"),
    namemap: str = typer.Option(..., "-m", "--namemap", help="Specify location of the tRNA name mapping file"),
    addtrna: Optional[str] = typer.Option(None, "--addtrna", help="Specify location of additional tRNA sequences file"),
    addseqs: Optional[str] = typer.Option(None, "--addseqs", help="Specify location of additional sequences file"),
    orgmode: str = typer.Option("euk", "-s", "--orgmode", help="Specify organism mode used for tRNAScan-SE"),
    forcecca: bool = typer.Option(False, "--forcecca", help="Force addition of CCA tail"),
    threads: int = typer.Option(0, "-n", "--threads", help="Specify number of threads to use (default: cpu_max)"),
    output: str = typer.Option("db", "-o", "--output", help="Specify output directory/name for bowtie2 index files"),
    log: Optional[str] = typer.Option(None, "--log", help="Log output to file"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
):
    args = SimpleNamespace(
        mode='makedb', genome=genome, trnaout=trnaout, trnafa=trnafa, namemap=namemap,
        addtrna=addtrna, addseqs=addseqs, orgmode=orgmode, forcecca=forcecca,
        threads=threads, output=output, log=log, quiet=quiet
    )
    run_logic(args)

@preprocess_app.command("trim", help="Trim, merge, and extract UMIs from fastq files using fastp")
def trim(
    runname: str = typer.Option(..., "-r", "--runname", help="Name of the run (used for output filenames)"),
    manifest: str = typer.Option(..., "-i", "--manifest", help="Tab-delimited file: SampleName <tab> R1_Path [<tab> R2_Path]"),
    adapter1: Optional[str] = typer.Option(None, "-a1", "--adapter1", help="Adapter sequence for R1 (optional, fastp auto-detects)"),
    adapter2: Optional[str] = typer.Option(None, "-a2", "--adapter2", help="Adapter sequence for R2 (optional, fastp auto-detects)"),
    length: int = typer.Option(15, "-l", "--length", help="Minimum length of sequence after trimming"),
    umilength: int = typer.Option(0, "-u", "--umilength", help="Length of UMI (0 to disable)"),
    umi3: bool = typer.Option(False, "--umi3", help="UMI is at the 3-prime end (Default is 5-prime)"),
    threads: int = typer.Option(0, "-n", "--threads", help="Total number of threads to use (0 = all available)"),
    log: Optional[str] = typer.Option(None, "--log", help="Log output to file"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
    verbose: bool = typer.Option(False, "-v", "--verbose", help="Print detailed command execution"),
):
    args = SimpleNamespace(
        mode='trim', runname=runname, manifest=manifest, adapter1=adapter1, adapter2=adapter2,
        length=length, umilength=umilength, umi3=umi3, threads=threads, log=log, quiet=quiet, verbose=verbose
    )
    run_logic(args)

@preprocess_app.command("map", help="Map reads to tRNA database")
def map_cmd(
    experiment: str = typer.Option(..., "-e", "--experiment", help="Experiment name to be used"),
    database: str = typer.Option(..., "-d", "--database", help="Name of the tRNA database"),
    samples: str = typer.Option(..., "-s", "--samples", help="Sample file"),
    lazy: bool = typer.Option(False, "--lazy", help="Skip mapping reads if bam files exist"),
    minnontrnasize: int = typer.Option(20, "--minnontrnasize", help="Minimum read length for non-tRNAs"),
    local: bool = typer.Option(False, "--local", help="Use local bam mapping"),
    threads: int = typer.Option(8, "-n", "--threads", help="Number of threads to use with Bowtie2 (default: 8)"),
    skipcheck: bool = typer.Option(False, "--skipcheck", help="Skips the check that the fq files match bam files"),
    bamdir: Optional[str] = typer.Option(None, "--bamdir", help="Directory for placing bam files (default: bam/<experimentname>)"),
    log: Optional[str] = typer.Option(None, "--log", help="Log output to file"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
):
    args = SimpleNamespace(
        mode='map', experiment=experiment, database=database, samples=samples,
        lazy=lazy, minnontrnasize=minnontrnasize, local=local, threads=threads, skipcheck=skipcheck,
        bamdir=bamdir, log=log, quiet=quiet
    )
    run_logic(args)

@app.command("build", help="Build a h5ad AnnData object from a tRAX run")
def build(
    traxdir: Optional[str] = typer.Option(None, "-i", "--traxdir", help="Specify location of trax directory (required if not running analysis)"),
    metadata: str = typer.Option(..., "-m", "--metadata", help="Specify a metadata file to create annotations, you can also use the sample file used to generate tRAX DB"),
    output: str = typer.Option("h5ad/trnagraph.h5ad", "-o", "--output", help="Specify output h5ad file"),
    # Analysis arguments
    experiment: Optional[str] = typer.Option(None, "-e", "--experiment", help="Experiment name (required for analysis)"),
    database: Optional[str] = typer.Option(None, "-d", "--database", help="Name of the tRNA database (required for analysis)"),
    samples: Optional[str] = typer.Option(None, "-s", "--samples", help="Sample file (required for analysis)"),
    gtf: Optional[str] = typer.Option(None, "--gtf", help="The ensembl gene list for that species"),
    pairs: Optional[str] = typer.Option(None, "--pairs", help="List of sample pairs to compare"),
    bed: Optional[List[str]] = typer.Option(None, "--bed", help="Additional bed files for feature list"),
    nofrag: bool = typer.Option(False, "--nofrag", help="Omit fragment determination (Used for TGIRT mapping)"),
    nosizefactors: bool = typer.Option(False, "--nosizefactors", help="Don't use Deseq size factors in plotting"),
    maxmismatches: Optional[str] = typer.Option(None, "--maxmismatches", help="Maximum allowed mismatches"),
    mincoverage: Optional[str] = typer.Option(None, "--mincoverage", help="Minimum read count for coverage plots"),
    minnontrnasize: int = typer.Option(20, "--minnontrnasize", help="Minimum read length for non-tRNAs"),
    hub: bool = typer.Option(False, "--hub", help="Make a track hub"),
    hubonly: bool = typer.Option(False, "--hubonly", help="Only make the track hub"),
    dumpother: bool = typer.Option(False, "--dumpother", help="Dump 'other' features when counting gene types"),
    bamdir: Optional[str] = typer.Option(None, "--bamdir", help="Directory for placing bam files (default: bam/<experimentname>)"),
    uniqueonly: bool = typer.Option(False, "--uniqueonly", help="Show only unique coverage"),
    dispfittype: str = typer.Option("mean", "--dispfittype", help="DESeq2 dispersion fit type: 'mean' (default, robust for small samples) or 'parametric' (for 5+ replicates per condition)"),
    threads: int = typer.Option(8, "-n", "--threads", help="Number of threads to use (default: 8)"),
    
    log: Optional[str] = typer.Option(None, "--log", help="Log output to file"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
):
    args = SimpleNamespace(
        mode='build', traxdir=traxdir, metadata=metadata, output=output,
        experiment=experiment, database=database, samples=samples, gtf=gtf, pairs=pairs,
        bed=bed, nofrag=nofrag, nosizefactors=nosizefactors, maxmismatches=maxmismatches,
        mincoverage=mincoverage, minnontrnasize=minnontrnasize, hub=hub, hubonly=hubonly,
        dumpother=dumpother, bamdir=bamdir, uniqueonly=uniqueonly, dispfittype=dispfittype, threads=threads,
        log=log, quiet=quiet
    )
    run_logic(args)

@app.command("merge", help="Merge data from two existing h5ad AnnData objects")
def merge(
    anndata1: str = typer.Option(..., "-i1", "--anndata1", help="Specify location of first h5ad object"),
    anndata2: str = typer.Option(..., "-i2", "--anndata2", help="Specify location of second h5ad object"),
    dropno: bool = typer.Option(False, "--dropno", help="Drop non tRNAs genes that are not present in both AnnData objects"),
    droprna: bool = typer.Option(False, "--droprna", help="Drop RNA categories that are not present in both AnnData objects"),
    output: str = typer.Option("h5ad/trnagraph.merge.h5ad", "-o", "--output", help="Specify output h5ad file"),
    log: Optional[str] = typer.Option(None, "--log", help="Log output to file"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
):
    args = SimpleNamespace(
        mode='merge', anndata1=anndata1, anndata2=anndata2, dropno=dropno, droprna=droprna,
        output=output, log=log, quiet=quiet
    )
    run_logic(args)

@app.command("cluster", help="Cluster data from an existing h5ad AnnData object")
def cluster(
    anndata: str = typer.Option(..., "-i", "--anndata", help="Specify location of h5ad object"),
    randomstate: Optional[int] = typer.Option(None, "-r", "--randomstate", help="Specify random state for UMAP if you want to have a static seed"),
    readcutoff: int = typer.Option(20, "-t", "--readcutoff", help="Specify readcount cutoff to use for clustering"),
    coveragetype: List[str] = typer.Option(['uniquecoverage', 'readstarts', 'readends', 'mismatchedbases', 'deletions'], "-v", "--coveragetype", help="Specify coverage types for umap clustering treated as features"),
    ncomponentsmp: int = typer.Option(2, "-c1", "--ncomponentsmp", help="Specify number of components to use for UMAP clustering of samples"),
    ncomponentgrp: int = typer.Option(2, "-c2", "--ncomponentgrp", help="Specify number of components to use for UMAP clustering of groups"),
    neighborclusmp: int = typer.Option(150, "-l1", "--neighborclusmp", help="Specify number of neighbors to use for UMAP clustering of samples"),
    neighborclusgrp: int = typer.Option(40, "-l2", "--neighborclusgrp", help="Specify number of neighbors to use for UMAP clustering of groups"),
    neighborstdsmp: int = typer.Option(75, "-n1", "--neighborstdsmp", help="Specify number of neighbors to use for UMAP projection plotting of samples"),
    neighborstdgrp: int = typer.Option(20, "-n2", "--neighborstdgrp", help="Specify number of neighbors to use for UMAP projection plotting of groups"),
    hdbscanminsampsmp: int = typer.Option(6, "-d1", "--hdbscanminsampsmp", help="Specify minsamples size to use for HDBSCAN clustering of samples"),
    hdbscanminsampgrp: int = typer.Option(3, "-d2", "--hdbscanminsampgrp", help="Specify minsamples size to use for HDBSCAN clustering of groups"),
    hdbscanminclusmp: int = typer.Option(30, "-b1", "--hdbscanminclusmp", help="Specify min cluster size to use for HDBSCAN clustering of samples"),
    hdbscanminclugrp: int = typer.Option(10, "-b2", "--hdbscanminclugrp", help="Specify min cluster size to use for HDBSCAN clustering of groups"),
    mindist: float = typer.Option(0.1, "-m", "--mindist", help="Specify minimum distance to use for UMAP clustering"),
    variancethreshold: float = typer.Option(0.1, "-e", "--variancethreshold", help="Specify variance threshold to use for feature selection"),
    umapstatsmetrics: str = typer.Option("euclidean", "-us", "--umapstatsmetrics", help="Specify UMAP statistics metrics to use for feature selection"),
    hdbstatsmetrics: str = typer.Option("euclidean", "-uh", "--hdbstatsmetrics", help="Specify hdbscan statistics metrics to use for feature selection with UMAP"),
    clusterobsexperimental: List[str] = typer.Option([], "--clusterobsexperimental", help="This is an experimental feature to add columns from adata.obs to the adata.var and adata.X to be used for clustering"),
    overwrite: bool = typer.Option(False, "-w", "--overwrite", help="Overwrite existing cluster information in AnnData object"),
    output: str = typer.Option("h5ad/trnagraph.cluster.h5ad", "-o", "--output", help="Specify output directory"),
    log: Optional[str] = typer.Option(None, "--log", help="Log output to file"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
):
    args = SimpleNamespace(
        mode='cluster', anndata=anndata, randomstate=randomstate, readcutoff=readcutoff, coveragetype=coveragetype,
        ncomponentsmp=ncomponentsmp, ncomponentgrp=ncomponentgrp, neighborclusmp=neighborclusmp, neighborclusgrp=neighborclusgrp,
        neighborstdsmp=neighborstdsmp, neighborstdgrp=neighborstdgrp, hdbscanminsampsmp=hdbscanminsampsmp, hdbscanminsampgrp=hdbscanminsampgrp,
        hdbscanminclusmp=hdbscanminclusmp, hdbscanminclugrp=hdbscanminclugrp, mindist=mindist, variancethreshold=variancethreshold,
        umapstatsmetrics=umapstatsmetrics, hdbstatsmetrics=hdbstatsmetrics, clusterobsexperimental=clusterobsexperimental,
        overwrite=overwrite, output=output, log=log, quiet=quiet
    )
    run_logic(args)

@app.command("graph", help="Graph data from an existing h5ad AnnData object")
def graph(
    anndata: str = typer.Option(..., "-i", "--anndata", help="Specify location of h5ad object"),
    output: str = typer.Option("figures", "-o", "--output", help="Specify output directory"),
    graphtypes: List[str] = typer.Option(["all"], "-g", "--graphtypes", help="Specify graphs to create, if not specified it will default to 'all'"),
    config: Optional[str] = typer.Option(None, "--config", help="Specify a json file containing observations/variables to filter out and other config options"),
    colormap: Optional[str] = typer.Option(None, "--colormap", help="Specify a json file containing colormaps for the graphs"),
    regen_uns: bool = typer.Option(False, "--regen_uns", help="Force regenerate uns log2fc data if it would be generated again"),
    threads: int = typer.Option(0, "-n", "--threads", help="Specify number of threads to use (default: cpu_max)"),
    log: Optional[str] = typer.Option(None, "--log", help="Log output to file"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
    verbose: bool = typer.Option(False, "-v", "--verbose", help="Print verbose output to stdout"),
    barcol: str = typer.Option("group", "--barcol", help="Specify AnnData column to of what the individal stacks of bars will be"),
    bargrp: str = typer.Option("amino", "--bargrp", help="Specify AnnData column to of what will stack in bar columns"),
    barsubgrp: Optional[str] = typer.Option(None, "--barsubgrp", help="Specify AnnData column for secondary spliting of bars into subplots"),
    barsort: Optional[str] = typer.Option(None, "--barsort", help="Specify AnnData column to sort the bars by"),
    barlabel: Optional[str] = typer.Option(None, "--barlabel", help="Specify wether to label the bars using a different AnnData column"),
    clustergrp: str = typer.Option("amino", "--clustergrp", help="Specify AnnData column to group by"),
    clusterlabels: Optional[str] = typer.Option(None, "--clusterlabels", help="Specify a AnnData column of names to use for the clusters instead of the default and will place them on the plot"),
    clusteroverview: bool = typer.Option(False, "--clusteroverview", help="Specify wether to generate an overview of the clusters"),
    clusternumeric: bool = typer.Option(False, "--clusternumeric", help="Specify wether to the cluster category is numeric"),
    clustermask: bool = typer.Option(False, "--clustermask", help="Specify wether to mask the cluster plots to annotated HDBSCAN clusters"),
    comparegrp1: str = typer.Option("group", "--comparegrp1", help="Specify AnnData column as main comparative group"),
    comparegrp2: str = typer.Option("group", "--comparegrp2", help="Specify AnnData column to group by"),
    corrmethod: str = typer.Option("pearson", "--corrmethod", help="Specify correlation method"),
    corrgroup: str = typer.Option("sample", "--corrgroup", help="Specify a grouping variable to generate correlation matrices for"),
    covgrp: str = typer.Option("group", "--covgrp", help="Specify a grouping variable to generate coverage plots for"),
    covobs: str = typer.Option("trna", "--covobs", help="Specify the basis for each individual coverage plot"),
    covtype: str = typer.Option("uniquecoverage", "--covtype", help="Specify a coverage type for coverage plots corresponding to trax coverage file outputs"),
    covgap: bool = typer.Option(False, "--covgap", help="Specify wether to include gaps in coverage plots"),
    covmethod: str = typer.Option("mean", "--covmethod", help="Specify method to use for coverage plots when combining multiple groups"),
    combinedpdfonly: bool = typer.Option(False, "--combinedpdfonly", help="Do not generate single tRNA coverage plot PDFs for every tRNA, only keep the combined output"),
    heatgrp: str = typer.Option("group", "--heatgrp", help="Specify group to use for heatmap"),
    diffrts: List[str] = typer.Option(['wholecounts_unique', 'fiveprime_unique', 'threeprime_unique', 'other_unique', 'total_unique'], "--diffrts", help="Specify readtypes to use for heatmap/volcano"),
    heatcutoff: int = typer.Option(80, "--heatcutoff", help="Specify readcount cutoff to use for heatmap"),
    heatbound: int = typer.Option(25, "--heatbound", help="Specify range to use for bounding the heatmap to top and bottom counts"),
    heatsubplots: bool = typer.Option(False, "--heatsubplots", help="Specify wether to generate subplots for each comparasion in addition to the sum"),
    pcamarkers: str = typer.Option("sample", "--pcamarkers", help="Specify AnnData column to use for PCA markers"),
    pcacolors: str = typer.Option("group", "--pcacolors", help="Specify AnnData column to color PCA markers by"),
    pcareadtypes: List[str] = typer.Option(['total_unique', 'total'], "--pcareadtypes", help="Specify read types to use for PCA markers"),
    radargrp: str = typer.Option("group", "--radargrp", help="Specify AnnData column to group by"),
    radarmethod: List[str] = typer.Option(['mean'], "--radarmethod", help="Specify method to use for radar plots"),
    radarscaled: bool = typer.Option(False, "--radarscaled", help="Specify wether to scale the radar plots to 100%% (optional)"),
    logogrp: str = typer.Option("amino", "--logogrp", help="Specify AnnData column to group sequences by"),
    logomanualgrp: Optional[List[str]] = typer.Option(None, "--logomanualgrp", help="Specify a manual group of tRNAs to use for seqlogo plots instead of using the AnnData column"),
    logomanualname: Optional[str] = typer.Option(None, "--logomanualname", help="Specify a name for the manual group of tRNAs output file, will be ignored and timestamped if not specified"),
    logopseudocount: int = typer.Option(20, "--logopseudocount", help="Specify the number of pseudocounts to add to each position when calculating as ratio of the bases in the pool of RNAs"),
    logosize: str = typer.Option("noloop", "--logosize", help="Specify the sequence size to use for the logo plots from presets"),
    ccatail: bool = typer.Option(True, "--ccatail", flag_value=False, help="Specify wether to keep the CCA tail from the sequences"),
    pseudogenes: bool = typer.Option(True, "--pseudogenes", flag_value=False, help="Specify wether to keep the pseudo-tRNAs (tRX)"),
    logornamode: bool = typer.Option(False, "--logornamode", help="Specify wether to print the output as RNA rather than DNA"),
    volgrp: str = typer.Option("group", "--volgrp", help="Specify group to use for volcano plot"),
    volcutoff: int = typer.Option(80, "--volcutoff", help="Specify readcount cutoff to use for volcano plot"),
):
    args = SimpleNamespace(
        mode='graph', anndata=anndata, output=output, graphtypes=graphtypes, config=config, colormap=colormap,
        regen_uns=regen_uns, threads=threads, log=log, quiet=quiet, verbose=verbose, barcol=barcol, bargrp=bargrp,
        barsubgrp=barsubgrp, barsort=barsort, barlabel=barlabel, clustergrp=clustergrp, clusterlabels=clusterlabels,
        clusteroverview=clusteroverview, clusternumeric=clusternumeric, clustermask=clustermask, comparegrp1=comparegrp1,
        comparegrp2=comparegrp2, corrmethod=corrmethod, corrgroup=corrgroup, covgrp=covgrp, covobs=covobs, covtype=covtype,
        covgap=covgap, covmethod=covmethod, combinedpdfonly=combinedpdfonly, heatgrp=heatgrp, diffrts=diffrts,
        heatcutoff=heatcutoff, heatbound=heatbound, heatsubplots=heatsubplots, pcamarkers=pcamarkers, pcacolors=pcacolors,
        pcareadtypes=pcareadtypes, radargrp=radargrp, radarmethod=radarmethod, radarscaled=radarscaled, logogrp=logogrp,
        logomanualgrp=logomanualgrp, logomanualname=logomanualname, logopseudocount=logopseudocount, logosize=logosize,
        ccatail=ccatail, pseudogenes=pseudogenes, logornamode=logornamode, volgrp=volgrp, volcutoff=volcutoff
    )
    run_logic(args)

@tools_app.command("log2fc", help="Compute log2fc data from an existing h5ad AnnData object")
def log2fc(
    anndata: str = typer.Option(..., "-i", "--anndata", help="Specify location of h5ad object"),
    group: str = typer.Option("group", "-g", "--group", help="Specify group to use for log2fc from obs"),
    readtypes: List[str] = typer.Option(['wholecounts_unique', 'fiveprime_unique', 'threeprime_unique', 'other_unique', 'total_unique'], "-r", "--readtypes", help="Specify readtypes to generate log2fc for"),
    cutoff: List[int] = typer.Option([80], "-x", "--cutoff", help="Specify readcounts cutoff to use for log2fc"),
    config: Optional[str] = typer.Option(None, "-c", "--config", help="Specify a json file containing observations/variables to filter out and other config options"),
    log: Optional[str] = typer.Option(None, "--log", help="Log output to file"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
):
    args = SimpleNamespace(
        mode='log2fc', anndata=anndata, group=group, readtypes=readtypes, cutoff=cutoff, config=config, log=log, quiet=quiet
    )
    run_logic(args)

@tools_app.command("csv", help="Output .h5ad to CSV")
def csv_cmd(
    anndata: str = typer.Option(..., "-i", "--anndata", help="Specify location of h5ad object"),
    output: str = typer.Option("csv", "-o", "--output", help="Specify output directory"),
    log: Optional[str] = typer.Option(None, "--log", help="Log output to file"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
):
    args = SimpleNamespace(
        mode='csv', anndata=anndata, output=output, log=log, quiet=quiet
    )
    run_logic(args)

@tools_app.command("test", help="Run pipeline demo tests")
def test(
    metadata: bool = typer.Option(False, "--metadata", help="Run metadata download test"),
    fastq: bool = typer.Option(False, "--fastq", help="Run fastq download test"),
    trna: bool = typer.Option(False, "--trna", help="Run tRNA download test"),
    genome: bool = typer.Option(False, "--genome", help="Run genome download test"),
    trim: bool = typer.Option(False, "--trim", help="Run trim test"),
    makedb: bool = typer.Option(False, "--makedb", help="Run makedb test"),
    map: bool = typer.Option(False, "--map", help="Run map test"),
    hubonly: bool = typer.Option(False, "--hubonly", help="Run map test with hubonly flag"),
    maponly: bool = typer.Option(False, "--maponly", help="Run map test with maponly flag"),
    build: bool = typer.Option(False, "--build", help="Run build test"),
    cluster: bool = typer.Option(False, "--cluster", help="Run cluster test"),
    merge: bool = typer.Option(False, "--merge", help="Run merge test"),
    graph: bool = typer.Option(False, "--graph", help="Run graph test"),
    all: bool = typer.Option(False, "--all", help="Run all tests (default)"),
    cleanrun: bool = typer.Option(False, "--cleanrun", help="Clean up test files after running tests"),
    directory: Optional[str] = typer.Option(None, "-d", "--directory", help="Specify directory to run tests in"),
    log: Optional[str] = typer.Option(None, "--log", help="Log output to file"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
):
    args = SimpleNamespace(
        mode='test', metadata=metadata, fastq=fastq, trna=trna, genome=genome, trim=trim,
        makedb=makedb, map=map, hubonly=hubonly, maponly=maponly, build=build, cluster=cluster, merge=merge, graph=graph, all=all, cleanrun=cleanrun, directory=directory, log=log, quiet=quiet
    )
    run_logic(args)

if __name__ == '__main__':
    app()
