#!/usr/bin/env python3

import os
import sys
import json
import shutil
import logging
import contextlib
import typer
from datetime import datetime
from typing import Optional, List
from types import SimpleNamespace

# Lazy imports are handled via the lazy_imports module to improve startup time
# These objects are proxies that only import the actual module when an attribute is accessed
try:
    from .modules.lazy_imports import (
        toolsMap, toolsTDatabase, toolsTrim, toolsTG,
        toolsTestSuite, adataGraph, adataMerge, adataCluster, adataBuild,
        anndata, matplotlib
    )
    from .modules import env_check
    from . import __version__
except ImportError:
    # Fallback for script execution: add parent directory to path and import as package
    import sys
    import os
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    from tRNAgraph.modules.lazy_imports import (
        toolsMap, toolsTDatabase, toolsTrim, toolsTG,
        toolsTestSuite, adataGraph, adataMerge, adataCluster, adataBuild,
        anndata, matplotlib
    )
    from tRNAgraph.modules import env_check
    from tRNAgraph import __version__


def _adata_basename(anndata_path: str) -> str:
    """Basename of an .h5ad path, extension stripped -- used to disambiguate log filenames for
    commands centered on a single anndata object, since multiple .h5ad files commonly live in
    the same directory."""
    return os.path.splitext(os.path.basename(anndata_path))[0]


class _Tee:
    """
    Duplicates writes to multiple underlying streams. Used to make print() calls (not yet
    converted to logging -- most of the codebase, per the roadmap's "Logging" item) land in
    both the persisted .log/ file and the console (unless --quiet), without needing every call
    site converted first.
    """
    def __init__(self, *streams):
        self._streams = [s for s in streams if s is not None]

    def write(self, data):
        for s in self._streams:
            s.write(data)

    def flush(self):
        for s in self._streams:
            s.flush()

    def isatty(self):
        return any(getattr(s, 'isatty', lambda: False)() for s in self._streams)


def configure_logging(log_path: str, quiet: bool) -> logging.Logger:
    """
    Configure the shared 'trnagraph' logger's handlers for one CLI invocation: a FileHandler on
    `log_path` (always attached -- file logging is unconditional, independent of --quiet, which
    only ever suppresses the console) plus a StreamHandler on the real console (only if not
    quiet, tagged `_is_console_handler` so toolsTG.progress_iterator() can find and temporarily
    detach exactly this one handler for a live rich display, without touching the FileHandler).
    Per Python's own logging documentation, library/module code should never configure its own
    handlers, only call `logging.getLogger(__name__)` and log -- messages then propagate up
    from e.g. 'trnagraph.modules.toolsTrim' to this 'trnagraph' logger for free.
    """
    logger = logging.getLogger('trnagraph')
    logger.setLevel(logging.INFO)
    logger.propagate = False
    for handler in logger.handlers[:]:
        handler.close()
        logger.removeHandler(handler)

    file_handler = logging.FileHandler(log_path, mode='w')
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(file_handler)

    if not quiet:
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setFormatter(logging.Formatter('%(message)s'))
        stream_handler._is_console_handler = True
        logger.addHandler(stream_handler)

    return logger

@contextlib.contextmanager
def handle_output(quiet: bool, tool: str, destination: Optional[str] = None, name_suffix: Optional[str] = None):
    """
    Context manager wrapping one CLI command's whole run. Always persists a timestamped log to
    ./.log/ (untracked by git; see .gitignore's existing `.log` entry), independent of --quiet,
    which only suppresses the console, never the file. On success, moves that log into
    `destination` (the command's real output location) if given. On any exception, prints a
    warning to stderr (visible even under --quiet, since a failed run should never be silent)
    pointing at the log still sitting in .log/, and deliberately skips the move -- the
    destination may not have been fully or correctly produced by a run that crashed.

    `name_suffix` (e.g. an input anndata's basename) is appended to the filename for commands
    where multiple objects could plausibly share one directory.
    """
    os.makedirs('.log', exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    filename = f'{timestamp}_{tool}' + (f'_{name_suffix}' if name_suffix else '') + '.log'
    log_path = os.path.join('.log', filename)

    real_stdout = sys.stdout
    logger = configure_logging(log_path, quiet)
    failed = False
    try:
        with open(log_path, 'a') as log_fileobj:
            tee = _Tee(log_fileobj, None if quiet else real_stdout)
            with contextlib.redirect_stdout(tee):
                yield
    except Exception:
        failed = True
        sys.stderr.write(f"WARNING: {tool} failed -- see {log_path} for details.\n")
        raise
    finally:
        for handler in logger.handlers[:]:
            handler.close()
            logger.removeHandler(handler)
        if not failed and destination:
            os.makedirs(destination, exist_ok=True)
            shutil.move(log_path, os.path.join(destination, filename))

app = typer.Typer(
    help="tRNAgraph is a tool for for advanced analysis of tRNA-seq data.",
    add_completion=False,
    no_args_is_help=True
)

def validate_environment():
    """
    Validates that the environment is set up correctly.
    """
    env_check.validate_environment()

def version_callback(value: bool):
    if value:
        typer.echo(f"tRNAgraph {__version__}")
        raise typer.Exit()

@app.callback()
def main_callback(
    skip_env_check: bool = typer.Option(False, "--skip-env-check", help="Skip environment validation checks"),
    version: Optional[bool] = typer.Option(
        None, "--version", callback=version_callback, is_eager=True, help="Show the tRNAgraph version and exit"
    )
):
    """
    tRNAgraph is a tool for for advanced analysis of tRNA-seq data.
    """
    if not skip_env_check:
        validate_environment()

preprocess_app = typer.Typer(help="Preprocess raw fastq/fasta files for tRNA analysis", no_args_is_help=True)
app.add_typer(preprocess_app, name="preprocess")

analyze_app = typer.Typer(help="Analyze tRNA-seq data (Build, Addsplit, Cluster)", no_args_is_help=True)
app.add_typer(analyze_app, name="analyze")

tools_app = typer.Typer(help="Extra utilities for working with tRNAgraph objects", no_args_is_help=True)
app.add_typer(tools_app, name="tools")

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
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
):
    # -o is a name prefix (e.g. "references/vibrChol1/trnadb/vibrChol1_db"), not itself a
    # directory -- the index files land in its dirname.
    destination = os.path.dirname(output) or "."
    with handle_output(quiet, tool="makedb", destination=destination):
        if not os.path.isfile(genome):
            raise Exception('Error: genome fasta file does not exist.')

        args = SimpleNamespace(
            mode='makedb', genome=genome, trnaout=trnaout, trnafa=trnafa, namemap=namemap,
            addtrna=addtrna, addseqs=addseqs, orgmode=orgmode, forcecca=forcecca,
            threads=threads, output=output, quiet=quiet
        )

        print('Building tRNA database...')
        toolsTDatabase.tRNADatabaseBuilder(args).main()
        print('Done!\n')

@preprocess_app.command("trim", help="Trim, merge, and extract UMIs from fastq files using fastp")
def trim(
    input: str = typer.Option(..., "-i", "--input", help="Tab-delimited manifest file: SampleName <tab> R1_Path [<tab> R2_Path]"),
    adapter1: Optional[str] = typer.Option(None, "-a1", "--adapter1", help="Adapter sequence for R1 (optional, fastp auto-detects)"),
    adapter2: Optional[str] = typer.Option(None, "-a2", "--adapter2", help="Adapter sequence for R2 (optional, fastp auto-detects)"),
    length: int = typer.Option(15, "-l", "--length", help="Minimum length of sequence after trimming"),
    umilength: int = typer.Option(0, "-u", "--umilength", help="Length of UMI (0 to disable)"),
    umi3: bool = typer.Option(False, "--umi3", help="UMI is at the 3-prime end (Default is 5-prime)"),
    threads: int = typer.Option(0, "-n", "--threads", help="Total number of threads to use (0 = all available)"),
    colormap: Optional[str] = typer.Option(None, "--colormap", help="Specify a json file containing colormaps for the trim stats plot (top-level 'trimtype' key)"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
    verbose: bool = typer.Option(False, "-v", "--verbose", help="Print detailed command execution"),
):
    # Same directory _generate_summary() writes trim_stats.csv/trim_feature_types.pdf to: the
    # first manifest entry's output-prefix directory, falling back to processed/trimmed exactly
    # like toolsTrim.py's own FastpTrimmer._construct_command() does.
    destination = "processed/trimmed"
    if os.path.isfile(input):
        with open(input, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    first_prefix = line.split()[0]
                    destination = os.path.dirname(first_prefix) or "processed/trimmed"
                    break
    with handle_output(quiet, tool="trim", destination=destination):
        if shutil.which('fastp') is None:
            raise Exception("Error: 'fastp' is not installed or not in PATH. Please install it (e.g., 'conda install -c bioconda fastp').")
        if not os.path.isfile(input):
            raise Exception(f'Error: Manifest file does not exist: {input}')

        args = SimpleNamespace(
            mode='trim', input=input, adapter1=adapter1, adapter2=adapter2,
            length=length, umilength=umilength, umi3=umi3, threads=threads, colormap=colormap,
            quiet=quiet, verbose=verbose
        )

        print('Starting fastp trimming pipeline...')
        toolsTrim.FastpTrimmer(args).process()
        print('Done!\n')

@preprocess_app.command("map", help="Map reads to tRNA database")
def map_cmd(
    output: str = typer.Option(..., "-o", "--output", help="Experiment name to be used"),
    database: str = typer.Option(..., "-d", "--database", help="Name of the tRNA database"),
    input: str = typer.Option(..., "-i", "--input", help="Specify a metadata file to create annotations"),
    lazy: bool = typer.Option(False, "--lazy", help="Skip mapping reads if bam files exist"),
    minnontrnasize: int = typer.Option(20, "--minnontrnasize", help="Minimum read length for non-tRNAs"),
    local: bool = typer.Option(False, "--local", help="Use local bam mapping"),
    threads: int = typer.Option(8, "-n", "--threads", help="Number of threads to use with Bowtie2 (default: 8)"),
    skipcheck: bool = typer.Option(False, "--skipcheck", help="Skips the check that the fq files match bam files"),
    bamdir: Optional[str] = typer.Option(None, "--bamdir", help="Directory for placing bam files (default: processed/<output>_bam)"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
):
    # -o is an "experiment name", not itself a directory map writes to -- the actual mapped BAM
    # output goes to --bamdir (or its default).
    destination = bamdir or f"processed/{output}_bam"
    with handle_output(quiet, tool="map", destination=destination):
        args = SimpleNamespace(
            mode='map', output=output, database=database, input=input,
            lazy=lazy, minnontrnasize=minnontrnasize, local=local, threads=threads, skipcheck=skipcheck,
            bamdir=bamdir, quiet=quiet
        )

        print('Mapping samples...')
        toolsMap.MapSamples(args).main()
        print('Done!\n')

@analyze_app.command("build", help="Build a h5ad AnnData object from a tRNAgraph preprocess run")
def build(
    input: str = typer.Option(..., "-i", "--input", help="Specify a metadata file to create annotations"),
    output: str = typer.Option("h5ad", "-o", "--output", help="Specify output directory (h5ad file named <dirname>.h5ad will be created inside)"),
    database: str = typer.Option(..., "-d", "--database", help="Name of the tRNA database"),
    # Analysis arguments
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
    filterother: bool = typer.Option(False, "--filterother", help="Dump reads counted in the 'other' type category (the 'other' row in typecounts.txt/uns['type_counts']) to a separate BAM file for inspection -- i.e. reads matching no tRNA, bed, or GTF-annotated feature"),
    bamdir: Optional[str] = typer.Option(None, "--bamdir", help="Directory for placing bam files (default: processed/<output>_bam)"),
    uniqueonly: bool = typer.Option(False, "--uniqueonly", help="Show only unique coverage"),
    dispfittype: str = typer.Option("parametric", "--dispfittype", help="DESeq2 dispersion fit type: 'parametric' (default) or 'mean' (robust for small samples)"),
    threads: int = typer.Option(8, "-n", "--threads", help="Number of threads to use (default: 8)"),
    readlengthsplit: Optional[int] = typer.Option(None, "-c", "--readlengthsplit", help="Read length cutoff for splitting (generates additional under/over analyses)"),
    overwritebams: bool = typer.Option(False, "--overwritebams", help="Force overwrite of existing BAM files during map/split"),
    savesplitbams: bool = typer.Option(False, "--savesplitbams", help="Keep the split BAM files (under --bamdir/u<N>,o<N>) created for --readlengthsplit instead of deleting them once merged into the AnnData object"),
    vst: str = typer.Option("log1p", "--vst", help="Variance Stabilizing Transformation method [vst, log1p, none]"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
):
    # Output is a directory path - h5ad filename is based on the directory basename
    output_dir = os.path.abspath(output)
    with handle_output(quiet, tool="build", destination=output_dir):
        if not os.path.isfile(input):
            raise Exception('Error: metadata file does not exist.')

        print('Building AnnData object...')
        h5ad_filename = os.path.basename(output_dir) + '.h5ad'
        print(toolsTG.builder(output_dir))

        # Update output to be the full h5ad path inside the directory
        full_output_path = os.path.join(output_dir, h5ad_filename)

        args = SimpleNamespace(
            mode='build', input=input, output=full_output_path,
            database=database, gtf=gtf, pairs=pairs,
            bed=bed, nofrag=nofrag, nosizefactors=nosizefactors, maxmismatches=maxmismatches,
            mincoverage=mincoverage, minnontrnasize=minnontrnasize, hub=hub, hubonly=hubonly,
            filterother=filterother, bamdir=bamdir, uniqueonly=uniqueonly, dispfittype=dispfittype, threads=threads,
            readlengthsplit=readlengthsplit, overwritebams=overwritebams, savesplitbams=savesplitbams,
            vst=vst, quiet=quiet
        )

        # Create AnnData object
        adataBuild.AnnDataBuilder(output_dir, input, full_output_path, args).create()
        print('Done!\n')

@analyze_app.command("addsplit", help="Add an additional read-length split variant (under/over cutoff pair) to an existing h5ad AnnData object")
def addsplit(
    anndata_path: str = typer.Option(..., "-i", "--anndata", help="Existing h5ad object to add a split variant to"),
    readlengthsplit: int = typer.Option(..., "-c", "--readlengthsplit", help="New read length cutoff to add (generates u<N>/o<N> variants)"),
    metadata: Optional[str] = typer.Option(None, "--metadata", help="Metadata file (default: recovered from the object's own build provenance)"),
    bamdir: Optional[str] = typer.Option(None, "--bamdir", help="Original (unsplit) BAM directory (default: recovered from provenance)"),
    database: Optional[str] = typer.Option(None, "-d", "--database", help="Override tRNA database (default: recovered from provenance)"),
    gtf: Optional[str] = typer.Option(None, "--gtf", help="Override GTF path (default: recovered from provenance)"),
    dispfittype: Optional[str] = typer.Option(None, "--dispfittype", help="Override DESeq2 dispersion fit type (default: recovered from provenance)"),
    vst: Optional[str] = typer.Option(None, "--vst", help="VST strategy for this split's vst layer (default: recovered from provenance)"),
    overwritebams: bool = typer.Option(False, "--overwritebams", help="Force overwrite of existing split BAM files"),
    savesplitbams: bool = typer.Option(False, "--savesplitbams", help="Keep the split BAM files (under --bamdir/u<N>,o<N>) instead of deleting them once merged into the AnnData object"),
    threads: int = typer.Option(8, "-n", "--threads", help="Number of threads to use (default: 8)"),
    output: Optional[str] = typer.Option(None, "-o", "--output", help="Output h5ad path (default: overwrite the input file in place)"),
    overwrite: bool = typer.Option(False, "-w", "--overwrite", help="Overwrite this cutoff's u<N>/o<N> data if already present in the object"),
    force: bool = typer.Option(False, "--force", help="Proceed even if explicitly-overridden parameters conflict with the object's original build provenance"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
):
    # add_split() itself resolves output_path = args.output or args.anndata (and writes there);
    # replicate that resolution here to know the destination before any work starts.
    destination = os.path.dirname(os.path.abspath(output or anndata_path))
    with handle_output(quiet, tool="addsplit", destination=destination, name_suffix=_adata_basename(anndata_path)):
        if not os.path.isfile(anndata_path):
            raise Exception('Error: h5ad file does not exist.')

        args = SimpleNamespace(
            mode='addsplit', anndata=anndata_path, readlengthsplit=readlengthsplit, metadata=metadata,
            bamdir=bamdir, database=database, gtf=gtf, dispfittype=dispfittype, vst=vst,
            overwritebams=overwritebams, savesplitbams=savesplitbams, threads=threads, output=output, overwrite=overwrite, force=force,
            quiet=quiet
        )

        print('Adding split variant to database object...\n')
        adataBuild.add_split(args)
        print('Done!\n')

@analyze_app.command("cluster", help="Cluster data from an existing h5ad AnnData object")
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
    variant: str = typer.Option("norm:full", "--variant", help="Select which normalization:split-tag to cluster, e.g. 'norm:u60'. norm is one of norm/raw/allfeatures/vst; tag is 'full' or an added split tag. Default 'norm:full' is today's default behavior"),
    overwrite: bool = typer.Option(False, "-w", "--overwrite", help="Overwrite existing cluster information in AnnData object"),
    output: str = typer.Option("trnagraph.cluster.h5ad", "-o", "--output", help="Specify output h5ad file path"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
):
    output_path = os.path.abspath(output)
    output_dir = os.path.dirname(output_path)
    with handle_output(quiet, tool="cluster", destination=output_dir or ".", name_suffix=_adata_basename(anndata)):
        if not os.path.isfile(anndata):
            raise Exception('Error: h5ad file does not exist.')

        if output_dir:
            print(toolsTG.builder(output_dir))

        args = SimpleNamespace(
            mode='cluster', anndata=anndata, randomstate=randomstate, readcutoff=readcutoff, coveragetype=coveragetype,
            ncomponentsmp=ncomponentsmp, ncomponentgrp=ncomponentgrp, neighborclusmp=neighborclusmp, neighborclusgrp=neighborclusgrp,
            neighborstdsmp=neighborstdsmp, neighborstdgrp=neighborstdgrp, hdbscanminsampsmp=hdbscanminsampsmp, hdbscanminsampgrp=hdbscanminsampgrp,
            hdbscanminclusmp=hdbscanminclusmp, hdbscanminclugrp=hdbscanminclugrp, mindist=mindist, variancethreshold=variancethreshold,
            umapstatsmetrics=umapstatsmetrics, hdbstatsmetrics=hdbstatsmetrics, clusterobsexperimental=clusterobsexperimental,
            variant=variant, overwrite=overwrite, output=output_path, quiet=quiet
        )
        
        print('Clustering data from database object...\n')
        adataCluster.anndataCluster(args).main()
        print('Done!\n')

@app.command("graph", help="Graph data from an existing h5ad AnnData object")
def graph(
    anndata: str = typer.Option(..., "-i", "--input", help="Specify location of h5ad object"),
    output: str = typer.Option("figures", "-o", "--output", help="Specify output directory"),
    graphtypes: List[str] = typer.Option(['all', 'cluster', 'correlation', 'count', 'coverage', 'heatmap', 'logo', 'pca', 'radar', 'volcano'], "-g", "--graphtypes", help="Specify graphs to create, if not specified it will default to 'all'"),
    config: Optional[str] = typer.Option(None, "--config", help="Specify a json file containing observations/variables to filter out and other config options"),
    colormap: Optional[str] = typer.Option(None, "--colormap", help="Specify a json file containing colormaps for the graphs"),
    regen_uns: bool = typer.Option(False, "--regen_uns", help="Force regenerate uns log2fc data if it would be generated again"),
    variant: str = typer.Option("norm:full", "--variant", help="Select which normalization:split-tag to plot, e.g. 'raw:full', 'norm:u60', 'allfeatures:o60'. norm is one of norm/raw/allfeatures/vst; tag is 'full' or an added split tag (e.g. 'u60'). Default 'norm:full' is today's default behavior"),
    threads: int = typer.Option(0, "-n", "--threads", help="Specify number of threads to use (default: cpu_max)"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
    verbose: bool = typer.Option(False, "-v", "--verbose", help="Print verbose output to stdout"),
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
    covtype: str = typer.Option("uniquecoverage", "--covtype", help="Specify a coverage type for coverage plots corresponding to coverage file outputs"),
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
    vollabels: Optional[int] = typer.Option(None, "--vollabels", help="Specify number of top significant markers to label on each volcano plot; omit to label all significant markers, or pass 0 to disable labels"),
):
    output_path = os.path.abspath(output)
    with handle_output(quiet, tool="graph", destination=output_path, name_suffix=_adata_basename(anndata)):
        # Set matplotlib backend to Agg to avoid display issues
        matplotlib.use('Agg')

        if not os.path.isfile(anndata):
            raise Exception('Error: h5ad file does not exist.')

        print(toolsTG.builder(output_path))

        args = SimpleNamespace(
            mode='graph', anndata=anndata, output=output_path, graphtypes=graphtypes, config=config, colormap=colormap,
            regen_uns=regen_uns, variant=variant, threads=threads, quiet=quiet, verbose=verbose, clustergrp=clustergrp, clusterlabels=clusterlabels,
            clusteroverview=clusteroverview, clusternumeric=clusternumeric, clustermask=clustermask, comparegrp1=comparegrp1,
            comparegrp2=comparegrp2, corrmethod=corrmethod, corrgroup=corrgroup, covgrp=covgrp, covobs=covobs, covtype=covtype,
            covgap=covgap, covmethod=covmethod, combinedpdfonly=combinedpdfonly, heatgrp=heatgrp, diffrts=diffrts,
            heatcutoff=heatcutoff, heatbound=heatbound, heatsubplots=heatsubplots, pcamarkers=pcamarkers, pcacolors=pcacolors,
            pcareadtypes=pcareadtypes, radargrp=radargrp, radarmethod=radarmethod, radarscaled=radarscaled, logogrp=logogrp,
            logomanualgrp=logomanualgrp, logomanualname=logomanualname, logopseudocount=logopseudocount, logosize=logosize,
            ccatail=ccatail, pseudogenes=pseudogenes, logornamode=logornamode, volgrp=volgrp, volcutoff=volcutoff,
            vollabels=vollabels
        )
        
        print('Graphing data from database object...\n')
        adataGraph.anndataGrapher(args).main()
        print('Done!\n')

@tools_app.command("log2fc", help="Compute log2fc data from an existing h5ad AnnData object")
def log2fc(
    anndata_path: str = typer.Option(..., "-i", "--anndata", help="Specify location of h5ad object"),
    group: str = typer.Option("group", "-g", "--group", help="Specify group to use for log2fc from obs"),
    readtypes: List[str] = typer.Option(['wholecounts_unique', 'fiveprime_unique', 'threeprime_unique', 'other_unique', 'total_unique'], "-r", "--readtypes", help="Specify readtypes to generate log2fc for"),
    cutoff: List[int] = typer.Option([80], "-x", "--cutoff", help="Specify readcounts cutoff to use for log2fc"),
    config: Optional[str] = typer.Option(None, "-c", "--config", help="Specify a json file containing observations/variables to filter out and other config options"),
    variant: str = typer.Option("norm:full", "--variant", help="Select which normalization:split-tag to compute log2fc for, e.g. 'norm:u60'. norm is one of norm/raw/allfeatures/vst; tag is 'full' or an added split tag. Default 'norm:full' is today's default behavior"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
):
    # log2fc always writes back into its own input file in place -- that file's directory is
    # the destination.
    destination = os.path.dirname(os.path.abspath(anndata_path))
    with handle_output(quiet, tool="log2fc", destination=destination, name_suffix=_adata_basename(anndata_path)):
        if not os.path.isfile(anndata_path):
            raise Exception('Error: h5ad file does not exist.')

        # Load the AnnData object
        # Note: using anndata.read_h5ad from lazy_imports
        adata = anndata.read_h5ad(anndata_path)

        # Load config file for name if specified
        config_name = 'default'
        config_data = None
        if config:
            with open(config, 'r') as f:
                config_data = json.load(f)
            if 'name' in config_data:
                config_name = config_data['name']
            else:
                raise ValueError('Config file must contain a "name" field')

        print('Calculating log2FC for database object...\n')
        # Resolve the requested normalization:split-tag into a working view so the readtype
        # obs-column lookups below transparently read the right variant's data -- see
        # toolsTG.build_variant_view() for why this view must never be written back directly.
        variant_spec = toolsTG.parse_variant(adata, variant)
        adata_view = toolsTG.build_variant_view(adata, variant_spec)

        for readtype in [f'nreads_{i}_norm' for i in readtypes]:
            for c in cutoff:
                toolsTG.adataLog2FC(adata_view, group, readtype, readcount_cutoff=c, config_name=config_name, overwrite=True).main()

        # Persist the newly-computed log2FC cache back onto the ORIGINAL (unresolved) adata,
        # into the correct namespaced location -- never write the resolved view itself, which
        # would overwrite the real full/default variant's data with the split variant's values.
        if variant_spec.tag == 'full':
            adata.uns['log2FC'] = adata_view.uns.get('log2FC', {})
        else:
            adata.uns.setdefault('size_splits', {}).setdefault(variant_spec.tag, {})['log2FC'] = adata_view.uns.get('log2FC', {})

        print('The log2FC uns dictionary has been updated.\nWriting h5ad database object to: ' + anndata_path)
        adata.write(anndata_path)
        print('Done!\n')

@tools_app.command("csv", help="Output .h5ad to CSV")
def csv_cmd(
    anndata_path: str = typer.Option(..., "-i", "--anndata", help="Specify location of h5ad object"),
    output: str = typer.Option("csv", "-o", "--output", help="Specify output directory"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
):
    # Add the name of the h5ad file to the output directory minus the extension
    output_path = os.path.abspath(output) + '/' + '.'.join(os.path.basename(anndata_path).split('.')[:-1]) + '/'
    with handle_output(quiet, tool="csv", destination=output_path, name_suffix=_adata_basename(anndata_path)):
        print(toolsTG.builder(output_path))

        adata = anndata.read_h5ad(anndata_path)
        print('Writing csv files to: ' + output_path)
        adata.write_csvs(output_path, skip_data=False)
        print('Done!\n')

@tools_app.command("merge", help="Merge data from two existing h5ad AnnData objects")
def merge(
    anndata1: str = typer.Option(..., "-i1", "--anndata1", help="Specify location of first h5ad object"),
    anndata2: str = typer.Option(..., "-i2", "--anndata2", help="Specify location of second h5ad object"),
    dropno: bool = typer.Option(False, "--dropno", help="Drop non tRNAs genes that are not present in both AnnData objects"),
    droprna: bool = typer.Option(False, "--droprna", help="Drop RNA categories that are not present in both AnnData objects"),
    output: str = typer.Option("trnagraph.merge.h5ad", "-o", "--output", help="Specify output h5ad file path"),
    force: bool = typer.Option(False, "--force", help="Proceed even if the two AnnData objects' build provenance (database/gtf) conflicts"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
):
    # Output is a h5ad file path - its parent directory is the destination. Named after the
    # newly-merged output file itself (not either input), since there's no single input identity.
    output_path = os.path.abspath(output)
    output_dir = os.path.dirname(output_path)
    with handle_output(quiet, tool="merge", destination=output_dir or ".", name_suffix=_adata_basename(output_path)):
        if not os.path.isfile(anndata1):
            raise Exception('Error: first h5ad file does not exist.')
        if not os.path.isfile(anndata2):
            raise Exception('Error: second h5ad file does not exist.')

        if output_dir:
            print(toolsTG.builder(output_dir))

        args = SimpleNamespace(
            mode='merge', anndata1=anndata1, anndata2=anndata2, dropno=dropno, droprna=droprna,
            output=output_path, force=force, quiet=quiet
        )

        print('Merging database objects...\n')
        adataMerge.anndataMerger(args).merge()
        print('Done!\n')

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
    build: bool = typer.Option(False, "--build", help="Run build test (no split)"),
    split_build: bool = typer.Option(False, "--split-build", help="Run build test with read length split"),
    cluster: bool = typer.Option(False, "--cluster", help="Run cluster test"),
    merge: bool = typer.Option(False, "--merge", help="Run merge test"),
    graph: bool = typer.Option(False, "--graph", help="Run graph test (no split)"),
    split_graph: bool = typer.Option(False, "--split-graph", help="Run graph test with read length split"),
    all: bool = typer.Option(False, "--all", help="Run all tests, forcing a clean workspace and full redownload"),
    skip_download: bool = typer.Option(False, "--skip-download", help="Skip metadata/fastq/tRNA/genome download steps and run everything else (downloads are already skipped by default when the target files are present; this forces it regardless)"),
    cleanrun: bool = typer.Option(False, "--cleanrun", help="Clean up test files after running tests"),
    directory: Optional[str] = typer.Option(None, "-d", "--directory", help="Specify directory to run tests in"),
    log: Optional[str] = typer.Option(None, "--log", help="Log output to file"),
    quiet: bool = typer.Option(False, "-q", "--quiet", help="Suppress output to stdout"),
):
    # Deliberately NOT wrapped in handle_output(): tools test keeps its own independent,
    # always-overwritten toolsTestSuite.log (unchanged, out of scope for the .log/ redesign --
    # see docs/roadmap.md). demoPipeline configures its own handlers entirely and has no bare
    # print() calls left for handle_output's stdout-tee to catch, so nothing here needs it.
    args = SimpleNamespace(
        mode='test', metadata=metadata, fastq=fastq, trna=trna, genome=genome, trim=trim,
        makedb=makedb, map=map, hubonly=hubonly, build=build,
        split_build=split_build, cluster=cluster, merge=merge, graph=graph, split_graph=split_graph,
        all=all, skip_download=skip_download, cleanrun=cleanrun, directory=directory, log=log, quiet=quiet
    )
    toolsTestSuite.demoPipeline(args).main()
    print('Done!\n')

if __name__ == '__main__':
    app()
