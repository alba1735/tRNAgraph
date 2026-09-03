import pandas as pd
import anndata as ad
import os
import shutil
import datetime
import subprocess
import logging
import numpy as np
from types import SimpleNamespace
from typing import Dict, Optional, Tuple
from multiprocessing import cpu_count
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import sklearn.preprocessing
from pydantic import ValidationError
from . import toolsTG
from .toolsSchemas import VariantContribution, MetadataFile, PairsFile
from .lazy_imports import toolsMap, toolsCountReads, toolsGetCoverage, toolsTrackHub

logger = logging.getLogger(__name__)


def _analysis_pipeline_phase_names():
    '''
    The fixed, named phase sequence AnalysisPipeline.run() executes -- excluding the `hubonly`
    early-return short-circuit, which isn't phase-tracked at all (a rare, effectively instant
    path with nothing worth reporting progress on). DESeq2 size factors are always computed
    (the removed `--nosizefactors` flag was confirmed broken -- see docs/roadmap.md), so
    "Analyzing counts"/"Analyzing unique counts" are unconditionally part of the sequence.
    '''
    return [
        "Counting Reads", "Analyzing counts", "Counting Read Types",
        "Analyzing unique counts", "Generating Read Coverage plots",
    ]


def _full_build_phase_names(analysis_args):
    '''
    The complete phase sequence for one `analyze build` run, spanning both AnalysisPipeline's
    per-variant counting/DESeq2/coverage phases and AnnDataBuilder.create()'s own X/layers/VST/
    write phases -- built upfront (everything here is knowable from `analysis_args` before any
    work starts) so a single shared PhaseTracker can report accurate cumulative progress across
    the whole command, not just its own class's slice of it. When --readlengthsplit is set, the
    AnalysisPipeline phase block repeats once per split variant (see AnnDataBuilder's settled
    design: a variant restarts/repeats the same phase sequence, labeled, rather than adding a
    third rendered progress level).
    '''
    vst_strategy = str(getattr(analysis_args, 'vst', 'vst')).lower()
    readlengthsplit = getattr(analysis_args, 'readlengthsplit', None)

    analysis_phases = _analysis_pipeline_phase_names()
    phases = list(analysis_phases)
    phases.append("Building AnnData object")
    if vst_strategy != 'none':
        phases.append("Computing VST")
    if readlengthsplit:
        phases.extend(analysis_phases)  # Under<N>
        phases.extend(analysis_phases)  # Over<N>
    phases.append("Writing h5ad")
    return phases


# The type labels in typecounts.txt/typerealcounts.txt that count as tRNA. toolsCountReads
# writes the first four as literals; every other row in those files is a GTF biotype out of
# emblbiotypes, a --bed feature set, or an extra-sequence class.
#
# Mt_tRNA is included even though it is GTF-driven rather than database-derived -- gtRNAdb and
# tRNAscan-SE exclude mitochondrial tRNAs, so a makedb-built database contains none and this row
# is structurally 0.0 today. Mitochondrial tRNAs are tRNAs, they are planned to become
# first-class (see docs/roadmap.md), and at roughly 60-75nt a read-length cutoff partitions them
# meaningfully, unlike the non-tRNA feature classes this filter targets. Keeping the data
# flowing now is cheaper than dropping it and re-adding it later.
TRNA_TYPE_LABELS = ('tRNA', 'tRNA_antisense', 'pretRNA', 'pretRNA_antisense', 'Mt_tRNA')


def _load_split_nontrna_read_counts(normalizedcounts_allfeatures_path, fallback_columns, fix_index):
    '''
    The non-tRNA feature counts for uns['nontRNA_counts'], read from the all-feature-normalized
    counts file.

    Non-tRNA features must not be normalized against tRNA/tRX-controlled size factors -- those
    represent the tRNA population, not the whole library -- so this reads the all-feature-
    controlled file rather than the primary one.

    That file does not exist for a read-length split variant, which excludes non-tRNA features
    and therefore runs no all-feature DESeq2 pass at all; the answer there is an empty frame over
    the sample axis. Reading it unconditionally is what crashed the first real split build after
    all-feature normalization became complete-variant only.
    '''
    if not os.path.exists(normalizedcounts_allfeatures_path):
        return pd.DataFrame(columns=list(fallback_columns))
    counts = fix_index(pd.read_csv(normalizedcounts_allfeatures_path, sep='\t', header=0))
    non_trna = counts[~counts.index.str.contains('tRNA')]
    return non_trna[~non_trna.index.str.contains('tRX')]


def _filter_rows_by_label(path, keep_label, has_header):
    '''
    Drop rows whose first tab-separated field fails `keep_label`, in place, returning how many
    were removed. A no-op when the file is absent or nothing matches.

    Filters line-wise rather than round-tripping through pandas, so every kept line is preserved
    byte for byte. That is not fastidiousness: these files use tRAX's convention where the header
    names only the data columns, one field shorter than the data rows, and readers depend on it
    -- `pd.read_csv(sep='\t')` with no index_col auto-detects the first column as the index only
    because of that raggedness. pandas' to_csv always writes a leading separator, which turned
    the labels into an ordinary string column and made plotsCount's `df*100/df.sum()` fail with
    "unsupported operand type(s) for /: 'str' and 'str'". `has_header` is explicit because
    genetypes.txt has no header line at all, and reading it with index_col=0 silently consumed
    its first real data row.
    '''
    if not os.path.exists(path):
        return 0
    with open(path) as handle:
        lines = handle.readlines()
    if not lines:
        return 0

    header = lines[:1] if has_header else []
    kept, removed = [], 0
    for line in lines[len(header):]:
        if not line.strip():
            continue
        if keep_label(line.split('\t', 1)[0].strip()):
            kept.append(line)
        else:
            removed += 1
    if removed:
        with open(path, 'w') as handle:
            handle.writelines(header + kept)
    return removed


def _load_mismatch_counts(path, size_factors):
    '''
    Load `<exp>-mismatches.txt` back into raw read counts for uns['mismatch_counts'].

    toolsCountReads.printmismatchcounts() writes this histogram -- one row per
    (mismatch count, tRNA/non-tRNA), one column per sample -- already divided by each
    sample's size factor. The AnnData object keeps the raw counts instead, so that
    division is undone here. The mismatch plot draws per-sample proportions, which the
    size factor cannot move either way; raw is stored because it is what was measured.

    Returns None when the file is absent, which is what a pre-existing object built
    before this key existed looks like to `graph`.
    '''
    if not os.path.isfile(path):
        return None
    df = pd.read_csv(path, sep='\t')
    for column in df.columns:
        if column in ('count', 'type'):
            continue
        factor = size_factors.get(column, 1.0) if hasattr(size_factors, 'get') else 1.0
        df[column] = (df[column] * factor).round().astype('int64')
    return df


def _trna_only_mismatch_counts(mismatch_counts):
    '''
    Drop the non-tRNA rows from a read-level mismatch histogram, for split variants.

    Read-length splits exclude non-tRNA features entirely (their length range is far wider
    than a tRNA's, so a shared cutoff divides the two pools on different terms), so a split's
    mismatch histogram must not carry a 'nontrna' series either.
    '''
    if mismatch_counts is None:
        return None
    return mismatch_counts[mismatch_counts['type'] != 'nontrna'].reset_index(drop=True)


def _filter_nontrna_rows_from_type_counts_file(path):
    '''
    Drop non-tRNA rows from a type-indexed count file (typecounts.txt/typerealcounts.txt).
    Same split-variant rationale as _filter_nontrna_rows_from_counts_file(), but keyed on the
    type-label vocabulary those files use instead of on feature names.
    '''
    return _filter_rows_by_label(path, lambda label: label in TRNA_TYPE_LABELS, has_header=True)


def _filter_nontrna_rows_from_counts_file(path, has_header=True):
    '''
    Drop non-tRNA feature rows from a feature-indexed count/annotation file.

    Used only for read-length split variants. A split cutoff partitions tRNA reads by design,
    but non-tRNA features are not classified by that criterion at all and span a far wider size
    range, so their per-split numbers record where the cutoff happened to fall rather than
    anything about the data -- on a real hg38 build the `other` biotype came out 12.41M under
    60nt against 0.99M over it. Filtering here, after the counting pass has written its files
    and before DESeq2 reads them, keeps toolsCountReads' classification path byte-identical for
    every variant; that module's ordering behaviour is validated against tRAX and is not worth
    perturbing for a reason unrelated to classification.

    Pass has_header=False for genetypes.txt, which has no header line.
    '''
    return _filter_rows_by_label(path, toolsTG.is_trna_feature, has_header)


class AnalysisPipeline:
    def __init__(self, args, expname=None, phase_tracker=None, variant_label=None, split_tag=None):
        self.logger = logging.getLogger(__name__)
        self.args = args
        self.dbname = args.database
        if expname:
            self.expname = expname
        else:
            self.expname = os.path.dirname(args.output)
        self.samplefilename = args.input
        self.ensgtf = args.gtf
        self.bedfiles = args.bed
        self.bamdir = args.bamdir if args.bamdir else os.path.join("processed", "bam")
        if args.threads:
            self.cores = args.threads
        else:
            try:
                self.cores = min(8, cpu_count())
            except Exception:
                self.cores = 8
        self.minnontrnasize = args.minnontrnasize
        self.maxmismatches = args.maxmismatches
        self.pairfile = args.pairs
        self.hubonly = args.hubonly
        self.makehubs = args.hub
        self.filterother = args.filterother
        self.dispfittype = getattr(args, 'dispfittype', 'mean')  # Default to 'mean' for robustness
        self.quiet = getattr(args, 'quiet', False)
        self.variant_label = variant_label
        # The split variant this pipeline is building ('u60'/'o60'), or None for the full/
        # complete variant. Distinct from variant_label, which is a human-readable phase-tracker
        # string ("Under 60") and must not be given control-flow meaning -- a later cosmetic
        # change to a log label would otherwise silently alter what gets analysed.
        self.split_tag = split_tag
        # A caller that spans a larger sequence (e.g. AnnDataBuilder, whose own phases -- and
        # potentially multiple split-variant repeats of this class's phases -- surround this
        # instance's slice of the work) passes its own shared tracker in; otherwise this instance
        # tracks just its own phases, so AnalysisPipeline stays usable standalone.
        self.phase_tracker = phase_tracker if phase_tracker is not None else toolsTG.PhaseTracker(
            phases=_analysis_pipeline_phase_names(), logger=self.logger, desc="Build",
        )

        self.trnainfo = toolsMap.trnadatabase(self.dbname)
        
        # Determine results and graphs directory names
        default_results_dir_name, default_graphs_dir_name = toolsTG.variant_dir_names(args)
        results_dir_name = getattr(args, 'results_dir_name', None) or default_results_dir_name
        graphs_dir_name = getattr(args, 'graphs_dir_name', None) or default_graphs_dir_name

        self.expinfo = toolsMap.expdatabase(self.expname, results_dir_name, graphs_dir_name)

    def run(self):
        # Create directories
        if not os.path.exists(self.expname):
            os.makedirs(self.expname)
        
        # graphsdir is never written to by this command -- nothing to create it for
        if not os.path.exists(self.expinfo.resultsdir):
            os.makedirs(self.expinfo.resultsdir)

        # Subdirectories in results. All-feature normalization is complete-variant only, so a
        # split variant gets no allfeature/ directory to leave empty.
        subdirs = ["mismatch", "pretRNAs", "unique", "trna"]
        if self.split_tag is None:
            subdirs.append("allfeature")
        for subdir in subdirs:
            path = os.path.join(self.expinfo.resultsdir, subdir)
            if not os.path.exists(path):
                os.makedirs(path)


        if self.hubonly:
            self.logger.info("Generating Track Hub...")
            self.createtrackhub()
            return

        self.makefeaturebed()

        # Counting Reads
        self.logger.info("Counting Reads")
        with self.phase_tracker.phase(variant=self.variant_label):
            self.countfeatures()
        # Before run_deseq2() reads these counts, so a split's size factors are never estimated
        # over non-tRNA features. Idempotent: the type-indexed files don't exist yet and are
        # skipped here, then filtered by the second call after counttypes() writes them.
        self.apply_split_variant_filters()

        # DESeq2 analysis using PyDESeq2
        self.logger.info("Analyzing counts")
        with self.phase_tracker.phase(variant=self.variant_label):
            self.run_deseq2()

        # Counting Read Types
        self.logger.info("Counting Read Types")
        with self.phase_tracker.phase(variant=self.variant_label):
            self.counttypes()
        self.apply_split_variant_filters()

        # Run DESeq2 on unique count files (they're generated by counttypes())
        self.logger.info("Analyzing unique counts")
        with self.phase_tracker.phase(variant=self.variant_label):
            self.run_unique_deseq2()

        # Coverage plots
        self.logger.info("Generating Read Coverage plots")
        orgtype = self.trnainfo.getorgtype()
        with self.phase_tracker.phase(variant=self.variant_label):
            self.gettrnacoverage(orgtype)

        if self.makehubs:
            self.createtrackhub()
            
        # Write run info
        self.write_runinfo()

    def write_runinfo(self):
        runinfo_file = os.path.join(self.expinfo.resultsdir, f"{os.path.basename(self.expinfo.expname)}-runinfo.txt")
        import time
        import subprocess
        
        try:
            # Assuming we are in a git repo
            git_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            git_version = subprocess.check_output(["git", "describe", "--always"], cwd=git_dir).decode().strip()
            git_hash = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=git_dir).decode().strip()
        except Exception:
            git_version = "unknown"
            git_hash = "unknown"

        with open(runinfo_file, "w") as f:
            f.write("Starting\n")
            f.write(f"expname\t{self.expname}\n")
            f.write(f"time\t{time.time()} ({time.strftime('%m/%d/%Y')})\n")
            f.write(f"samplefile\t{os.path.abspath(self.samplefilename)}\n")
            f.write(f"dbname\t{os.path.abspath(self.dbname)}\n")
            f.write(f"git version\t{git_version}\n")
            f.write(f"git version hash\t{git_hash}\n")
            # Reconstruct command line roughly (can be improved if we pass args string)
            cmd = f"trnagraph analyze build --experiment={self.expname} --database={self.dbname} --samples={self.samplefilename}"
            if self.pairfile:
                cmd += f" --pairs={self.pairfile}"
            f.write(f"command\t{cmd}\n")
            f.write("\n")

    def makefeaturebed(self):
        allfeatfile = open(self.expinfo.allfeats, "w")
        for currfeature in toolsTG.readbed(self.trnainfo.maturetrnas):
            print(currfeature.bedstring(), file=allfeatfile)
        for currfeature in toolsTG.readbed(self.trnainfo.locifile):
            print(currfeature.bedstring(), file=allfeatfile)
        if self.ensgtf:
            for currfeature in toolsTG.readgtf(self.ensgtf):
                print(currfeature.bedstring(name = currfeature.data["genename"]), file=allfeatfile)
        if self.bedfiles:
            for currbed in self.bedfiles:
                for currfeature in toolsTG.readbed(currbed):
                    print(currfeature.bedstring(), file=allfeatfile)
        allfeatfile.close()

    def countfeatures(self):
        toolsCountReads.countreads_main(samplefile=self.samplefilename, ensemblgtf=self.ensgtf,
                            maturetrnas=[self.trnainfo.maturetrnas], bamdir=self.bamdir,
                            otherseqs=self.trnainfo.otherseqs, trnaloci=[self.trnainfo.locifile],
                            removepseudo=True, genetypefile=self.expinfo.genetypes,
                            trnatable=self.trnainfo.trnatable, countfile=self.expinfo.genecounts,
                            bedfile=self.bedfiles if self.bedfiles else [], trnacounts=self.expinfo.trnacounts,
                            trnaends=self.expinfo.trnaendfile, trnauniquecounts=self.expinfo.trnauniquefile,
                            cores=self.cores, maxmismatches=self.maxmismatches, quiet=self.quiet)

    def apply_split_variant_filters(self):
        '''
        Strip non-tRNA rows from this variant's count outputs, for split variants only.

        Runs after the counting passes have written their files and before DESeq2 reads them, so
        toolsCountReads' classification path stays byte-identical for every variant. A no-op for
        the complete variant, and a no-op with no --gtf, where no non-tRNA feature was ever
        counted in the first place.
        '''
        if self.split_tag is None:
            return
        removed = 0
        # genecounts (readcounts.txt) carries a header naming the samples; genetypes.txt has no
        # header line at all and starts straight into data.
        removed += _filter_nontrna_rows_from_counts_file(self.expinfo.genecounts, has_header=True)
        removed += _filter_nontrna_rows_from_counts_file(self.expinfo.genetypes, has_header=False)
        for path in (self.expinfo.genetypecounts, self.expinfo.genetyperealcounts):
            removed += _filter_nontrna_rows_from_type_counts_file(path)
        if removed:
            self.logger.info(
                f"Split variant '{self.split_tag}': removed {removed} non-tRNA rows from the "
                f"count outputs. A read-length cutoff cannot meaningfully partition non-tRNA "
                f"features, so they are carried only by the complete variant."
            )

    def run_deseq2(self):
        # 1. Main Counts - tRNA/tRX-controlled size factors (default; drives adata.X/obs/raw)
        self.run_deseq2_on_file(
            counts_file=self.expinfo.genecounts,
            norm_counts_file=self.expinfo.normalizedcounts,
            size_factors_file=self.expinfo.sizefactors,
            output_dir=self.expinfo.resultsdir,
            prefix="",
            use_trna_control=True
        )

        # 1b. Main Counts - all-feature size factors (secondary, kept for comparison).
        # Complete-variant only. A split variant has had its non-tRNA features removed, so
        # "all features" and "tRNA only" would be the same set and this pass would merely
        # duplicate the primary one; computing it from the unfiltered counts instead -- which is
        # what used to happen -- estimated size factors over a non-tRNA pool that a tRNA-length
        # cutoff had arbitrarily truncated, and then normalized real tRNA data with them via the
        # user-reachable `--variant allfeatures:<tag>`.
        if self.split_tag is None:
            self.run_deseq2_on_file(
                counts_file=self.expinfo.genecounts,
                norm_counts_file=self.expinfo.normalizedcounts_allfeatures,
                size_factors_file=self.expinfo.allfeaturesizefactors,
                output_dir=os.path.join(self.expinfo.resultsdir, "allfeature"),
                prefix="allfeature_",
                use_trna_control=False
            )

        # 2. tRNA Counts
        self.run_deseq2_on_file(
            counts_file=self.expinfo.trnacounts,
            norm_counts_file=self.expinfo.normalizedtrnacounts,
            size_factors_file=self.expinfo.trnasizefactors,
            output_dir=os.path.join(self.expinfo.resultsdir, "trna"),
            prefix="trna_"
        )
    
    def run_unique_deseq2(self):
        """Run DESeq2 on unique count files (must be called after counttypes())"""
        # 3. Unique tRNA Counts
        self.run_deseq2_on_file(
            counts_file=self.expinfo.trnauniqcountsfile,
            norm_counts_file=os.path.join(self.expinfo.resultsdir, "unique", f"{os.path.basename(self.expinfo.expname)}-uniquetrnas_normalizedreadcounts.txt"),
            size_factors_file=os.path.join(self.expinfo.resultsdir, "unique", f"{os.path.basename(self.expinfo.expname)}-uniquetrnas_SizeFactors.txt"),
            output_dir=os.path.join(self.expinfo.resultsdir, "unique"),
            prefix="uniquetrnas_"
        )
        
        # 4. Unique Amino Counts
        self.run_deseq2_on_file(
            counts_file=self.expinfo.trnauniqaminofile,
            norm_counts_file=os.path.join(self.expinfo.resultsdir, "unique", f"{os.path.basename(self.expinfo.expname)}-uniqueaminos_normalizedreadcounts.txt"),
            size_factors_file=os.path.join(self.expinfo.resultsdir, "unique", f"{os.path.basename(self.expinfo.expname)}-uniqueaminos_SizeFactors.txt"),
            output_dir=os.path.join(self.expinfo.resultsdir, "unique"),
            prefix="uniqueaminos_"
        )
        
        # 5. Unique Anticodon Counts
        self.run_deseq2_on_file(
            counts_file=self.expinfo.trnauniqanticodonfile,
            norm_counts_file=os.path.join(self.expinfo.resultsdir, "unique", f"{os.path.basename(self.expinfo.expname)}-uniqueanticodons_normalizedreadcounts.txt"),
            size_factors_file=os.path.join(self.expinfo.resultsdir, "unique", f"{os.path.basename(self.expinfo.expname)}-uniqueanticodons_SizeFactors.txt"),
            output_dir=os.path.join(self.expinfo.resultsdir, "unique"),
            prefix="uniqueanticodons_"
        )

    def run_deseq2_on_file(self, counts_file, norm_counts_file, size_factors_file, output_dir, prefix, use_trna_control=False):
        # Load counts

        if not os.path.exists(counts_file):
            self.logger.warning(f"Warning: Counts file {counts_file} not found. Skipping DESeq2 for {prefix}.")
            return

        counts_df = pd.read_csv(counts_file, sep='\t', index_col=0)
        # Transpose because PyDESeq2 expects samples as rows
        counts_df = counts_df.T
        
        # Load sample info
        try:
            # Use sep=None to auto-detect delimiter (handles tabs or spaces)
            sample_df = pd.read_csv(self.samplefilename, sep=None, engine='python', header=None)
            
            # Drop header if present
            if str(sample_df.iloc[0, 0]).lower() == 'fastq':
                sample_df = sample_df.iloc[1:]
                
            # Select columns 1 and 2
            if len(sample_df.columns) >= 3:
                sample_df = sample_df.iloc[:, [1, 2]]
                sample_df.columns = ['sample', 'condition']
                sample_df['replicate'] = sample_df['condition'] # Assign replicate same as condition
            else:
                 self.logger.error("Error: Metadata file must have at least 3 columns (fastq, sample, group)")
                 return

            sample_df.set_index('sample', inplace=True)
        except Exception as e:
            self.logger.error(f"Error reading sample file {self.samplefilename}: {e}")
            return

        # Filter samples that are in counts
        # Note: counts_df.index might be sample names, ensure they match sample_df.index
        common_samples = counts_df.index.intersection(sample_df.index)
        if len(common_samples) < len(counts_df.index):
             self.logger.warning(f"Warning: Some samples in counts file {counts_file} are missing from sample file.")
        
        counts_df = counts_df.loc[common_samples]
        sample_df = sample_df.loc[common_samples]
        
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Run DESeq2
        if counts_df.empty or (counts_df.sum().sum() == 0):
            self.logger.warning(f"Warning: Counts matrix {counts_file} is empty or all zeros. Skipping DESeq2.")
            # Create dummy size factors
            pd.DataFrame(1.0, index=sample_df.index, columns=['sizeFactor']).T.to_csv(size_factors_file, sep=' ', index=False)
            # Create empty normalized counts
            counts_df.T.to_csv(norm_counts_file, sep='\t')
            return

        try:
            control_genes = None
            if use_trna_control:
                # Select only features containing 'tRNA' or 'tRX' for computing size factors
                trna_features = [f for f in counts_df.columns if toolsTG.is_trna_feature(f)]
                if trna_features:
                    control_genes = trna_features
                    self.logger.info(f"Using {len(control_genes)} tRNA features for size factor calculation.")
                else:
                    self.logger.warning("Warning: No tRNA features found for size factor computation. Defaulting to all features.")

            dds = DeseqDataSet(counts=counts_df, metadata=sample_df, design_factors="condition", fit_type=self.dispfittype, control_genes=control_genes)
            dds.deseq2()


        except Exception as e:
            self.logger.warning(f"Warning: DESeq2 failed for {prefix}: {e}")
            # Create dummy size factors
            pd.DataFrame(1.0, index=sample_df.index, columns=['sizeFactor']).T.to_csv(size_factors_file, sep=' ', index=False)
            # Save raw counts as normalized counts (fallback)
            counts_df.T.to_csv(norm_counts_file, sep='\t')
            return
        
        # Save normalized counts
        if 'normed_counts' in dds.layers:
            norm_counts = dds.layers['normed_counts']
            norm_counts = pd.DataFrame(norm_counts, index=dds.obs_names, columns=dds.var_names)
            norm_counts = norm_counts.T # Transpose back
            norm_counts.to_csv(norm_counts_file, sep='\t')
        else:
            self.logger.warning("Warning: 'normed_counts' not found in dds.layers. Using raw counts.")
            counts_df.T.to_csv(norm_counts_file, sep='\t')
        
        # Save size factors
        if 'size_factors' in dds.obs:
            size_factors = dds.obs['size_factors']
            pd.DataFrame(size_factors.values, index=dds.obs_names, columns=['sizeFactor']).T.to_csv(size_factors_file, sep=' ', index=False)
        elif 'size_factors' in dds.obsm:
            size_factors = dds.obsm['size_factors']
            pd.DataFrame(size_factors.values, index=dds.obs_names, columns=['sizeFactor']).T.to_csv(size_factors_file, sep=' ', index=False)
        else:
             self.logger.warning("Warning: 'size_factors' not found in dds.obs or dds.obsm. Using 1.0.")
             pd.DataFrame(1.0, index=dds.obs_names, columns=['sizeFactor']).T.to_csv(size_factors_file, sep=' ', index=False)

        # Run pairwise comparisons if pairs file is provided
        de_results = {}
        if self.pairfile:
            de_results = self.run_pairwise_de(dds, sample_df, output_dir)
        
        # Generate aggregated files (padjs, logvals, avgs, etc.)
        self.write_aggregated_de_results(dds, sample_df, de_results, output_dir, prefix)

    def run_pairwise_de(self, dds, sample_df, output_dir):
        de_results = {}
        try:
            # Use sep='\s+' to handle any whitespace
            pairs_df = pd.read_csv(self.pairfile, sep=r'\s+', engine='python', header=None, names=['Sample1', 'Sample2'])
            PairsFile(path=self.pairfile, pairs=list(pairs_df.itertuples(index=False, name=None)))
        except Exception as e:
            self.logger.error(f"Error reading pairs file: {e}")
            return de_results

        # Ensure output directory exists
        de_out_dir = os.path.join(output_dir, "de_results")
        os.makedirs(de_out_dir, exist_ok=True)

        for index, row in pairs_df.iterrows():
            sample1 = row['Sample1']
            sample2 = row['Sample2']
            
            # Check if samples exist in metadata or are valid conditions
            is_sample1 = sample1 in sample_df.index
            is_sample2 = sample2 in sample_df.index
            
            is_cond1 = sample1 in sample_df['condition'].values
            is_cond2 = sample2 in sample_df['condition'].values
            
            if not (is_sample1 or is_cond1) or not (is_sample2 or is_cond2):
                self.logger.warning(f"{sample1} or {sample2} not found in metadata (as sample or condition). Skipping pair.")
                continue

            if is_sample1:
                cond1 = sample_df.loc[sample1, 'condition']
            else:
                cond1 = sample1
                
            if is_sample2:
                cond2 = sample_df.loc[sample2, 'condition']
            else:
                cond2 = sample2
            
            if cond1 == cond2:
                self.logger.warning(f"Comparison between same condition {cond1}. Skipping DE.")
                continue

            self.logger.info(f"Running DE for {cond1} vs {cond2}...")
            
            try:
                stat_res = DeseqStats(dds, contrast=["condition", cond1, cond2])
                stat_res.summary()
                res_df = stat_res.results_df
                
                # Save results
                out_file = os.path.join(de_out_dir, f"{cond1}_vs_{cond2}.txt")
                res_df.to_csv(out_file, sep='\t')
                self.logger.info(f"Saved DE results to {out_file}")
                
                # Store for aggregation
                comp_name = f"{cond1}_{cond2}" # Match tRAX naming convention if possible, or just use cond1_cond2
                de_results[comp_name] = res_df
                
            except Exception as e:
                self.logger.error(f"Error running DE for {cond1} vs {cond2}: {e}")
                
        return de_results

    def write_aggregated_de_results(self, dds, sample_df, de_results, output_dir, prefix):
        """
        Generates aggregated output files:
        - padjs.txt
        - logvals.txt (log2FoldChange)
        - avgs.txt (mean normalized counts per condition)
        - medians.txt (median normalized counts per condition)
        - dispersions.txt
        - combine.txt
        """
        self.logger.info(f"Generating aggregated DE results for {prefix}...")
        
        # 1. Dispersions
        # PyDESeq2 stores dispersions in dds.var DataFrame
        dispersions = None
        if hasattr(dds, 'var') and 'dispersions' in dds.var.columns:
            dispersions = dds.var['dispersions'].values
        elif hasattr(dds, 'varm') and 'dispersions' in dds.varm:
            dispersions = dds.varm['dispersions']
        elif hasattr(dds, 'varm') and 'genewise_dispersions' in dds.varm:
            dispersions = dds.varm['genewise_dispersions']
        elif hasattr(dds, 'dispersions'):
            dispersions = dds.dispersions
        
        if dispersions is not None:
            disp_file = os.path.join(output_dir, f"{os.path.basename(self.expinfo.expname)}-{prefix}dispersions.txt")
            # Write without header, space-separated to match tRAX format
            with open(disp_file, 'w') as f:
                for gene_name, disp_val in zip(dds.var_names, dispersions):
                    f.write(f"{gene_name} {disp_val}\n")
        else:
            self.logger.warning(f"Could not find dispersions in dds object for {prefix}")

        
        # 2. Normalized Counts Stats (Avgs, Medians)
        if 'normed_counts' in dds.layers:
            norm_counts = pd.DataFrame(dds.layers['normed_counts'], index=dds.obs_names, columns=dds.var_names).T
        else:
            # Fallback to raw counts if normalization failed (shouldn't happen if dds succeeded)
            norm_counts = pd.DataFrame(dds.X, index=dds.obs_names, columns=dds.var_names).T

        conditions = sample_df['condition'].unique()
        avgs_df = pd.DataFrame(index=norm_counts.index)
        medians_df = pd.DataFrame(index=norm_counts.index)
        
        for cond in conditions:
            samples = sample_df[sample_df['condition'] == cond].index
            # Filter samples that are in norm_counts columns
            valid_samples = [s for s in samples if s in norm_counts.columns]
            if valid_samples:
                avgs_df[cond] = norm_counts[valid_samples].mean(axis=1)
                medians_df[cond] = norm_counts[valid_samples].median(axis=1)
                
        # tRAX's -avgs.txt is DESeq2's baseMean -- the mean of normalized counts across ALL
        # samples, invariant to the contrast requested. analyzecounts.R's colgetavgname() takes
        # column 1 of each results() object and writes it once per pairwise comparison, so tRAX's
        # columns are labelled per comparison while every one of them holds the same numbers.
        # Emit the quantity once under its real name rather than porting that duplication.
        # Computed from the normalized counts directly (DESeq2's own definition) rather than
        # lifted out of a results() object, so it is still written when no comparison ran.
        basemean_df = pd.DataFrame({'baseMean': norm_counts.mean(axis=1)})

        avgs_file = os.path.join(output_dir, f"{os.path.basename(self.expinfo.expname)}-{prefix}avgs.txt")
        # The per-group means are the more useful quantity of the two and are kept -- under a
        # name that does not collide with tRAX's, which holds baseMean above. No tRAX file
        # corresponds to this one.
        groupavgs_file = os.path.join(output_dir, f"{os.path.basename(self.expinfo.expname)}-{prefix}groupavgs.txt")
        medians_file = os.path.join(output_dir, f"{os.path.basename(self.expinfo.expname)}-{prefix}medians.txt")
        basemean_df.to_csv(avgs_file, sep='\t')
        avgs_df.to_csv(groupavgs_file, sep='\t')
        medians_df.to_csv(medians_file, sep='\t')
        
        # 3. Aggregated DE Results (padjs, logvals)
        if not de_results:
            return

        padjs_df = pd.DataFrame(index=norm_counts.index)
        logvals_df = pd.DataFrame(index=norm_counts.index)
        
        for comp_name, res_df in de_results.items():
            # res_df index should match norm_counts index
            # Rename columns to comparison name
            if 'padj' in res_df.columns:
                padjs_df[comp_name] = res_df['padj']
            if 'log2FoldChange' in res_df.columns:
                logvals_df[comp_name] = res_df['log2FoldChange']
            
        padjs_file = os.path.join(output_dir, f"{os.path.basename(self.expinfo.expname)}-{prefix}padjs.txt")
        logvals_file = os.path.join(output_dir, f"{os.path.basename(self.expinfo.expname)}-{prefix}logvals.txt")
        padjs_df.to_csv(padjs_file, sep='\t')
        logvals_df.to_csv(logvals_file, sep='\t')
        
        # 4. Combine File
        # Format: log2_Comp1, log2_Comp2, ..., pval_Comp1, ..., Cond1_avg, Cond2_avg...
        combine_df = pd.DataFrame(index=norm_counts.index)
        
        # Add log2FC
        for col in logvals_df.columns:
            combine_df[f"log2_{col}"] = logvals_df[col]
            
        # Add pvals
        for col in padjs_df.columns:
            combine_df[f"pval_{col}"] = padjs_df[col]
            
        # Add per-group medians. tRAX's -combine.txt is cbind(alllogvals, allprobs, medcountmat)
        # -- its trailing per-group columns are medians, matching its own -medians.txt exactly.
        # These were means here, under the same group labels, so a column-aligned diff against
        # tRAX compared a mean against a median without either side looking wrong.
        for col in medians_df.columns:
            combine_df[col] = medians_df[col]
            
        combine_file = os.path.join(output_dir, f"{os.path.basename(self.expinfo.expname)}-{prefix}combine.txt")
        combine_df.to_csv(combine_file, sep='\t')

    def counttypes(self):
        toolsCountReads.main(sizefactors=self.expinfo.sizefactors, combinereps=True,
                            bamdir=self.bamdir, otherseqs=self.trnainfo.otherseqs,
                            samplefile=self.samplefilename, maturetrnas=[self.trnainfo.maturetrnas],
                            trnatable=self.trnainfo.trnatable, trnaaminofile=self.expinfo.trnaaminofile,
                            trnaanticodonfile=self.expinfo.trnaanticodonfile, ensemblgtf=self.ensgtf,
                            trnaloci=[self.trnainfo.locifile], countfile=self.expinfo.genetypecounts,
                            realcountfile=self.expinfo.genetyperealcounts, mismatchfile=self.expinfo.mismatchcountfile,
                            bedfile=self.bedfiles, readlengthfile=self.expinfo.trnalengthfile,
                            countfrags=False, bamnofeature=self.filterother,
                            uniquename=self.expinfo.uniquename, cores=self.cores,
                            quiet=self.quiet)

    def gettrnacoverage(self, orgtype):
        toolsGetCoverage.main(samplefile=self.samplefilename, bedfile=[self.trnainfo.maturetrnas],
                             locibed=[self.trnainfo.locifile], locistk=self.trnainfo.locialign,
                             bamdir=self.bamdir, lociedgemargin=30, sizefactors=self.expinfo.sizefactors,
                             orgtype=orgtype, locicoverage=self.expinfo.locicoveragefile,
                             stkfile=self.trnainfo.trnaalign, numfile=self.trnainfo.trnanums,
                             locinums=self.trnainfo.locinums, allcoverage=self.expinfo.trnacoveragefile,
                             trnafasta=self.trnainfo.trnafasta, cores=self.cores,
                             uniqcoverage=self.expinfo.trnauniqcoveragefile,
                             # Hardcoded, not a parameter: tRAX's getcoverage.py ignores its own
                             # --uniqueonly flag and always passes True, so coverage is unique-reads-only
                             # there. Turning this back into a configurable value is what commit dca2673
                             # did by accident, silently changing every coverage file (~1.28x tRAX on
                             # hg38). The unique/multi breakdown survives regardless, in the coverage
                             # table's own columns and in obs's nreads_*_unique_* vs nreads_* pairs.
                             uniqueonly=True, sigmismatchbed=self.expinfo.sigmismatchbed,
                             sigmismatchtable=self.expinfo.sigmismatchfile)

    def createtrackhub(self):
        hub_builder = toolsTrackHub.TrackHubBuilder(
            genomedatabase=self.trnainfo.dbname,
            samplefilename=self.samplefilename,
            expname=self.expinfo.expname,
            threads=self.cores
        )
        hub_builder.run()


def _precompute_default_log2fc(view, threads=None):
    '''
    Precompute log2FC for the 'group' comparison's two overview readtypes
    (nreads_total_unique_norm, nreads_total_norm) across the standard tRNA-seq read-count
    cutoffs, giving graph-time volcano/heatmap runs at these common combos an immediate cache
    hit instead of needing an on-demand fit. Shared by _adata_build_() (the full/default
    variant) and merge_variant_into_adata() (each split variant), so both stay in sync.

    Called once per variant, outside any multiprocessing.Pool, so -- unlike the identically-
    shaped on-demand calls inside plotsVolcano.py/plotsHeatmap.py, which run inside pooled
    worker processes and must stay at adataLog2FC's safe default of n_cpus=1 to avoid a
    nested-process-pool deadlock -- it's safe to use real parallelism here. `threads` should be
    the same --threads budget already governing the rest of this build/addsplit command, not an
    independent "use everything" default.
    '''
    for cutoff in [20, 40, 80, 100, 200]:  # common read cutoffs for tRNAseq
        toolsTG.adataLog2FC(view, 'group', 'nreads_total_unique_norm', readcount_cutoff=cutoff, config_name='default', overwrite=True, n_cpus=threads).main()
        toolsTG.adataLog2FC(view, 'group', 'nreads_total_norm', readcount_cutoff=cutoff, config_name='default', overwrite=True, n_cpus=threads).main()


# The acceptor-stem regions the obs['fragment'] heuristic averages coverage over. Named
# explicitly, rather than tracking whatever `location` happens to call them, because splitting
# 73-76 out into 'threeprime_end' would otherwise drop the discriminator base and CCA tail from
# the 3' mean -- and the CCA tail is exactly where 3' fragments pile up, so fragment calls would
# shift as a side effect of a renaming. These tuples preserve the pre-rename membership; whether
# the heuristic *should* exclude the tail is a separate question, deliberately not settled here.
FRAGMENT_FIVEPRIME_REGIONS = ('fiveprime_extra', 'fiveprime_acceptorstem')
FRAGMENT_THREEPRIME_REGIONS = ('threeprime_acceptorstem', 'threeprime_end')


def _sprinzl_location_maps() -> Tuple[Dict[str, str], Dict[str, str], Dict[str, Optional[str]]]:
    '''
    Map each Sprinzl position to its structural region, the half of the tRNA it falls in, and
    the region's short code -- written into `adata.var['location']`, `['half']` and
    `['location_code']` respectively.

    The region definitions are ported from tRNAscan-SE's SprinzlPos.pm `ss_pos` table (Patricia
    Chan, UCSC), the same source the Lowe Lab's own tooling uses, so tRNAgraph's annotation
    agrees with tRNAscan-SE, tRAX and tDRnamer rather than with a parallel invention. Short
    codes are that table's own (5P1/3P1/L1/... ); the descriptive names it pairs them with are
    "5p Acceptor Stem", "Acceptor-D-arm-linker", "Variable Stem" and so on.

    Notable consequences relative to what tRNAgraph used to emit:
      - The acceptor stem is 1-7 paired with 66-72, symmetric. It previously ran -1..7 against
        66..76, folding the discriminator base and the CCA tail into a seven-base-pair stem.
      - 73-76 become their own region ('threeprime_end', L7): 73 is the discriminator base and
        74-76 the CCA tail, neither of which is part of the stem.
      - The 5' extra base (-1, e.g. the G-1 of tRNA-His) gets its own region rather than being
        counted as a stem position. It is absent from bacterial position tables entirely, so
        this region simply has no rows under `--orgmode bact`.
      - 44-48 are the variable LOOP and the e-series positions the variable STEM -- the reverse
        of the old 'anticodon_to_t_internal'/'extensionloop' naming.
      - The D- and T-stems are split into 5'/3' halves, matching how the acceptor and anticodon
        stems were already handled.

    Every key must be a `str`: `var['positions']` holds strings, so an int key silently matches
    nothing and the position comes out NaN.
    '''
    loc_dict: Dict[str, str] = {}
    loc_half_dict: Dict[str, str] = {}
    loc_code_dict: Dict[str, Optional[str]] = {}

    def assign(positions, region, half, code):
        for position in positions:
            loc_dict[position] = region
            loc_half_dict[position] = half
            loc_code_dict[position] = code

    def rng(start, stop):
        return [str(i) for i in range(start, stop)]

    # 5' extra base. SprinzlPos.pm's position list starts at 1, so there is no canonical code
    # for it; inventing one would misrepresent the scheme, so it is left unset.
    assign(['-1'], 'fiveprime_extra', 'fiveprime', None)
    # Acceptor stem: 1-7 pairs with 66-72 (1:72, 2:71, ... 7:66).
    assign(rng(1, 8), 'fiveprime_acceptorstem', 'fiveprime', '5P1')
    assign(rng(66, 73), 'threeprime_acceptorstem', 'threeprime', '3P1')
    # Discriminator base (73) and CCA tail (74-76).
    assign(rng(73, 77), 'threeprime_end', 'threeprime', 'L7')
    # Acceptor stem to D-arm linker.
    assign(['8', '9'], 'a_to_d_internal', 'fiveprime', 'L1')
    # D-arm: stem 10-13 pairs with 22-25, loop between them.
    assign(rng(10, 14), 'fiveprime_dstem', 'fiveprime', '5P2')
    assign(rng(22, 26), 'threeprime_dstem', 'fiveprime', '3P2')
    assign(rng(14, 22) + ['17a', '20a', '20b'], 'dloop', 'fiveprime', 'L2')
    # D-arm to anticodon-arm linker.
    assign(['26'], 'd_to_anticodon_internal', 'fiveprime', 'L3')
    # Anticodon arm: stem 27-31 pairs with 39-43, loop 32-38 between them.
    assign(rng(27, 32), 'fiveprime_anticodonstem', 'center', '5P3')
    assign(rng(39, 44), 'threeprime_anticodonstem', 'center', '3P3')
    assign(rng(32, 39), 'anticodonloop', 'center', 'L4')
    # Variable arm: 44-48 are the loop, the e-series positions the stem.
    assign(rng(44, 49), 'variableloop', 'threeprime', 'L5')
    assign(['e' + str(i) for i in range(1, 20)], 'variablestem', 'threeprime', 'P4')
    # T-arm: stem 49-53 pairs with 61-65, T-Psi-C loop 54-60 between them.
    assign(rng(49, 54), 'fiveprime_tstem', 'threeprime', '5P5')
    assign(rng(61, 66), 'threeprime_tstem', 'threeprime', '3P5')
    assign(rng(54, 61), 'tloop', 'threeprime', 'L6')

    return loc_dict, loc_half_dict, loc_code_dict


class AnnDataBuilder():
    '''
    Create h5ad AnnData object
    '''
    def __init__(self, resultsdir, metadata, output, analysis_args=None, results_dir_name=None, graphs_dir_name=None):
        '''
        Initialize AnnDataBuilder object.

        `results_dir_name`/`graphs_dir_name` let this be used as a "loader-only" instance
        (analysis_args=None) that reads an already-generated results/graphs directory --
        e.g. a size-split variant's `results/u60`/`graphs/u60` -- without re-running the
        analysis pipeline. This is how split-variant contributions get computed for merging
        into an existing AnnData object, both at initial build time (_apply_readlength_split_)
        and via the incremental `analyze addsplit` command.
        '''
        self.logger = logging.getLogger(__name__)
        self.analysis_args = analysis_args
        self.metadata_path = metadata
        # Populated by _handle_preprocessing() when splitting runs -- tracks which u{N}/o{N}
        # split-BAM directories already had every sample's file present *before* this run's own
        # splitting step, so _apply_readlength_split_()'s cleanup only deletes what this run
        # actually created.
        self._split_dirs_preexisted = {}

        # Parse and validate the metadata/samples file FIRST, before any expensive/risky
        # pipeline work below (AnalysisPipeline's counting/coverage generation) runs -- a
        # malformed metadata file (e.g. a duplicate sample row) previously wasn't caught until
        # deep inside that pipeline run, as an opaque pandas reshape/pivot error.
        metadata_type = '\t' if metadata.endswith('.tsv') else ',' if metadata.endswith('.csv') else None
        try:
            self.metadata = pd.read_csv(metadata, sep=metadata_type, header=None, index_col=None, engine='python' if metadata_type is None else None)
        except Exception:
            raise ValueError(f'Could not read metadata file, check to make sure it is formated correctly: {metadata}')

        # Validate and drop fastq column
        if len(self.metadata.columns) < 3:
             raise ValueError(f'Metadata file must have at least 3 columns (fastq, sample, group): {metadata}')
        # Drop the first column (fastq)
        self.metadata = self.metadata.iloc[:, 1:]
        # turn first row into a list then remove it from the dataframe
        first_row = self.metadata.iloc[0].dropna().values.tolist()
        if len(first_row) >= 2 and first_row[0] == 'sample' and first_row[1] == 'group':
            self.observations = first_row
            self.metadata = self.metadata.iloc[1:]
        else:
            self.observations = ['sample', 'group'] + [f'metadata_{i}' for i in range(2, len(first_row))]

        # Make sure that observations are not going to be generated automatically
        auto_obs = ['trna', 'iso', 'amino', 'deseq2_sizefactor', 'refseq', 'dataset', 'pseudogene']
        if any(x in self.observations for x in auto_obs):
            raise ValueError(f'The following observation categories will automatically be generated please remove these if you included them: {auto_obs}')
        # Make sure that observations are unique
        if len(self.observations) != len(set(self.observations)):
            raise ValueError(f'Observation categories must be unique, please remove duplicates from the observation catgories you wish to generate: {self.observations}')
        # Make sure that sample and group are the first two observations
        if self.observations[0] != 'sample' or self.observations[1] != 'group':
            raise ValueError(f'The first two observation categories must be "sample" and "group" please reorder your observation categories to match the following: ["sample", "group", ...]: {self.observations}')
        # Add manual observations to obs list if they are not provided or if the length of the observations list does not match the number of parameters in the coverage file
        if len(self.observations) != len(self.metadata.columns):
            diff_obs_count = len(self.metadata.columns)-len(self.observations)
            self.logger.warning(f'Number of observations does not match number of parameters in coverage file by {diff_obs_count}. To create a more specific database object, please provide the correct number of observations.')
            if diff_obs_count > 0:
                self.logger.info(f'Adding {diff_obs_count} observations to the end of the list')
                self.observations += ['obs_' + str(x) for x in range(diff_obs_count)]
            if diff_obs_count < 0:
                self.logger.info(f'Removing {abs(diff_obs_count)} observations from the end of the list')
                self.observations = self.observations[:diff_obs_count]
        # Validate the parsed rows (duplicate/empty sample names, no data rows) at read time,
        # before they're indexed into self.metadata below -- a duplicate sample name previously
        # went undetected here and silently produced duplicate adata.obs rows downstream.
        try:
            MetadataFile(
                path=metadata,
                header=self.observations,
                rows=[[None if pd.isna(x) else str(x) for x in row] for row in self.metadata.itertuples(index=False)],
            )
        except ValidationError as e:
            raise ValueError(f'Invalid metadata file {metadata}:\n{e}') from e
        # Add align observation categories to metadata
        self.metadata.columns = self.observations
        self.metadata.set_index('sample', inplace=True)
        # Sample names as given in the metadata file, kept before the to_dict() conversion below
        # so create() can check them against the samples actually present in the coverage file.
        self.metadata_samples = set(self.metadata.index)
        self.metadata = self.metadata.to_dict()

        # Run analysis pipeline if args are provided
        if analysis_args:
            # Handle auto-mapping and splitting if needed
            self._handle_preprocessing(analysis_args)

            # One shared tracker spans the WHOLE build -- this class's own later phases in
            # create() (Building AnnData object / Computing VST / Writing h5ad), the main
            # AnalysisPipeline run below, and (via _apply_readlength_split_) each split-variant's
            # own repeat of AnalysisPipeline's phases -- so the reported percentage reflects
            # total build progress, not just whichever phase happens to be running.
            self.phase_tracker = toolsTG.PhaseTracker(
                phases=_full_build_phase_names(analysis_args), logger=self.logger, desc="Build",
            )

            self.logger.info("Running analysis pipeline (Full Dataset)...")
            pipeline = AnalysisPipeline(analysis_args, phase_tracker=self.phase_tracker)
            pipeline.run()
        else:
            self.phase_tracker = None

        # Initialize expdatabase to get file paths (toolsTG.variant_dir_names tolerates
        # analysis_args being None -- getattr(None, 'readlengthsplit', None) is falsy, so it
        # just falls through to the plain 'results'/'graphs' default, same as no split requested)
        default_results_dir_name, default_graphs_dir_name = toolsTG.variant_dir_names(analysis_args)
        if results_dir_name is None:
            results_dir_name = getattr(analysis_args, 'results_dir_name', None) or default_results_dir_name
        if graphs_dir_name is None:
            graphs_dir_name = getattr(analysis_args, 'graphs_dir_name', None) or default_graphs_dir_name

        self.expinfo = toolsMap.expdatabase(resultsdir, results_dir_name, graphs_dir_name)

        # Add unique feature column to coverage file for alignment and sorting
        coverageresults = self.expinfo.trnacoveragefile
        self.coverage = pd.read_csv(coverageresults, sep='\t', header=0)
        self.coverage['uniquefeat'] = self.coverage['Feature'] + '_' + self.coverage['Sample']
        self.positions = pd.unique(self.coverage['position'])
        # Add size factors to coverage file
        sizefactors = self.expinfo.sizefactors
        self.size_factors = pd.read_csv(sizefactors, sep=" ", header=0).to_dict('index')[0]
        self.size_factors_list = None
        # Secondary, all-feature-controlled size factors kept for comparison against the
        # tRNA-controlled default. Absent for a read-length split variant, which excludes
        # non-tRNA features and so runs no all-feature DESeq2 pass at all.
        if os.path.exists(self.expinfo.allfeaturesizefactors):
            self.size_factors_allfeatures = pd.read_csv(self.expinfo.allfeaturesizefactors, sep=" ", header=0).to_dict('index')[0]
        else:
            self.size_factors_allfeatures = None
        # For adding unique counts to coverage file
        trnauniquecounts = self.expinfo.trnauniqcountsfile #'-trnauniquecounts.txt'
        self.unique_counts = pd.read_csv(trnauniquecounts, sep='\t', header=0).to_dict('index')
        # For adding normalized read counts to coverage file split by read type
        def _fix_index(df):
            # Fallback for when the feature/gene index isn't auto-detected on read: PyDESeq2's
            # var_names has no .name set, so to_csv writes it as an unnamed column, which
            # pd.read_csv then reads back as a numeric-dtype "Unnamed: 0" column instead of the index.
            if df.index.dtype != 'object':
                if 'Unnamed: 0' in df.columns:
                    df = df.set_index('Unnamed: 0')
                    df.index.name = None
                elif df.shape[1] > 0 and df.iloc[:, 0].astype(str).str.contains('tRNA').any():
                    df = df.set_index(df.columns[0])
                    df.index.name = None
            return df

        normalizedreadcounts = self.expinfo.normalizedcounts
        normalized_read_counts = _fix_index(pd.read_csv(normalizedreadcounts, sep='\t', header=0))

        # Non-tRNA features must not be normalized against tRNA/tRX-controlled size
        # factors -- those are only representative of the tRNA population, not the whole library.
        # Use the all-feature-controlled normalized counts (same combined counts file, but
        # size factors estimated over all features) so adata.uns['nontRNA_counts'] is on a
        # statistically appropriate scale for non-tRNA analysis.
        self.non_trna_read_counts = _load_split_nontrna_read_counts(
            self.expinfo.normalizedcounts_allfeatures, normalized_read_counts.columns, _fix_index)

        # Clean all non tRNAs from normalized read counts by removing all rows that do not have a tRNA in the feature column
        normalized_read_counts = normalized_read_counts[(normalized_read_counts.index.str.contains('tRNA')) | (normalized_read_counts.index.str.contains('tRX'))]
        self.read_types = pd.unique(normalized_read_counts.index.str.split('_').str[1])
        self.normalized_read_counts = normalized_read_counts.to_dict('index')
        # For adding anticoodon counts to coverage file
        anticodoncounts = self.expinfo.trnaanticodonfile
        self.anticodon_counts = pd.read_csv(anticodoncounts, sep='\t')
        # For adding type counts to coverage file
        typecounts = self.expinfo.genetypecounts
        self.type_counts = pd.read_csv(typecounts, sep='\t')
        typerealcounts = self.expinfo.genetyperealcounts
        self.type_real_counts = pd.read_csv(typerealcounts, sep='\t')
        # For adding amino acid counts to coverage file
        aminoacidcounts = self.expinfo.trnaaminofile
        self.amino_counts = pd.read_csv(aminoacidcounts, sep='\t')
        # Read-level mismatch histogram, stored raw -- see _load_mismatch_counts()
        self.mismatch_counts = _load_mismatch_counts(self.expinfo.mismatchcountfile, self.size_factors)
        # Create list of reference sequences from actualbase column of coverage file - skips gap positions
        self.seqs = self._seq_build_()
        self.seqs_full = self._seq_build_(gap=True)
        # Add trnagraph version to run info based on github hash - Will be changed to git describe once the package is deployed
        git_version, git_hash = _get_git_version_()

        # Capture the CLI flags used to build this object, sanitizing None (anndata's h5ad uns writer
        # doesn't reliably round-trip None values, so substitute a string sentinel like adataCluster.py does)
        # `cli_specified` is bookkeeping for the --config merge (a frozenset of the option
        # names the user actually typed), not a setting this object was built with -- and
        # h5ad cannot write a frozenset, so recording it fails the whole build at the point
        # the object is saved.
        run_flags = {}
        if self.analysis_args:
            run_flags = toolsTG.sanitize_flags(self.analysis_args)
        self.trnagraph_run_info = {'expname': resultsdir.split('/')[-1],
                                   'time': os.popen('date').read().rstrip(),
                                   'trnagraph_directory': resultsdir,
                                   'git version': git_version,
                                   'git version hash': git_hash,
                                   'flags': run_flags}
        # Output file name
        self.output = output
        # Names of coverage types to add to adata object from coverage file
        self.cov_types = ['coverage', 'readstarts', 'readends', 'uniquecoverage', 'multitrnacoverage',
                          'multianticodoncoverage', 'multiaminocoverage','tRNAreadstotal', 'mismatchedbases',
                          'deletedbases', 'adenines', 'thymines', 'cytosines', 'guanines', 'deletions']

    def create(self):
        '''
        Create h5ad database object
        '''
        # Build obs and x dataframes
        with self.phase_tracker.phase():
            x_df, self.size_factors_list = self._build_coverage_matrix_()
            obs_df = self._obs_build_(x_df)
            self._check_metadata_sample_match_(obs_df)
            # obs_df,x_df = self._group_sort_(obs_df,x_df) # Not sure if I need this function
            # Check that the index of the obs and x dataframes are the same
            if not obs_df.index.equals(x_df.index):
                raise ValueError('The index of the obs and x dataframes are not the same. This means somthing went wrong in the sorting process.')
            # Build adata object
            adata = self._adata_build_(obs_df, x_df)

            # Add size factors to adata object as raw layer
            adata.layers['raw'] = adata.X * adata.obs['deseq2_sizefactor'].values[:,None]

            # Add the secondary, all-feature-controlled normalization as an independent layer for comparison
            allfeature_sf = np.array([self.size_factors_allfeatures.get(s, 1.0) for s in obs_df['sample'].values])
            adata.layers['norm_allfeatures'] = adata.layers['raw'] / allfeature_sf[:, None]

        # Determine vst strategy
        vst_strategy = 'vst'
        if self.analysis_args and hasattr(self.analysis_args, 'vst'):
            vst_strategy = str(self.analysis_args.vst).lower()

        if vst_strategy != 'none':
            minfeaturereads = getattr(self.analysis_args, 'minfeaturereads', None)
            minfeaturereads = int(minfeaturereads) if minfeaturereads is not None else 30
            fit_mask = self._vst_fit_mask_(obs_df, minfeaturereads)
            adata.obs['vst_fit_excluded'] = obs_df['vst_fit_excluded'].values
            with self.phase_tracker.phase():
                adata.layers['vst'] = self._compute_vst_(
                    adata.X, adata.layers['raw'], adata.obs['deseq2_sizefactor'].values, adata.obs.get('group', 'all'),
                    vst_strategy, getattr(self.analysis_args, 'dispfittype', 'parametric'), adata.obs_names, adata.var_names,
                    fit_mask=fit_mask,
                )

        # Quality check adata by dropping NaN values and printing summary
        if adata.obs.isna().any(axis=0).any():
            self.logger.warning('NaN values found in obs dataframe this is commonly caused by missing samples in your metadata file or havinge a different number of observations per sample.\n' + \
                  'Another consideration is to make sure your .tsv/.csv is using tabs or commas appropriately.\n' + \
                  'It can also be caused by metadata containing NaN or None values. Please check your metadata file to make sure the following are correct\n' + \
                  f'Observations:\n{str(adata.obs.columns[adata.obs.isna().any(axis=0)].tolist())}\n' + \
                  f'Samples with NaN:\n{str(list(set(adata.obs["sample"][adata.obs.isna().any(axis=1)].tolist())))}\n')
        # Add output name to adata object index
        adata.obs.index = [os.path.basename(self.output).split('.')[0] + '_' + str(x) for x in adata.obs.index]

        # Declared category order, applied before the split so every variant shares one ordering
        # (identity columns are shared across variants; only numeric columns are per-variant).
        # Last, so it can name columns derived late in the build such as 'fragment'.
        order = getattr(self.analysis_args, 'order', None) if self.analysis_args else None
        if order:
            toolsTG.apply_category_order(adata.obs, order)
            for column, levels in order.items():
                self.logger.info(f'Ordered obs column {column!r}: {list(levels)} '
                                 f'(reference level: {levels[0]!r})')

        # Apply read-length split variants (added as layers/obsm/uns onto this SAME object,
        # rather than the old behavior of writing separate _u{N}.h5ad/_o{N}.h5ad files)
        if self.analysis_args and getattr(self.analysis_args, 'readlengthsplit', None):
            self._apply_readlength_split_(adata)

        # Save adata object
        with self.phase_tracker.phase():
            toolsTG.write_h5ad(adata, self.output)
            self.logger.info(f'Writing h5ad database object to {self.output}')

    def _check_metadata_sample_match_(self, obs_df):
        '''
        Validate that sample names in the metadata file line up with the samples actually
        present in the coverage file (i.e. the samples that were mapped/counted by the
        pipeline), before the AnnData object is written.

        A coverage sample missing from the metadata is fatal: `_obs_build_` looks it up with
        `dict.get()`, so it silently fills every metadata column with NaN for that sample's
        rows instead of failing -- this is the exact failure mode this check exists to catch
        early, with an actionable message, instead of downstream as an opaque NaN warning.
        A metadata sample with no matching coverage sample is not fatal (e.g. a metadata file
        intentionally kept as a superset across runs) and is only reported as a warning.
        '''
        coverage_samples = set(obs_df['sample'].unique())
        missing_from_metadata = sorted(coverage_samples - self.metadata_samples)
        unused_in_metadata = sorted(self.metadata_samples - coverage_samples)

        if missing_from_metadata:
            raise ValueError(
                'Metadata Check failed: the following sample(s) appear in the coverage/count '
                f'output but not in the metadata file ({self.metadata_path}), which would '
                f'otherwise silently produce NaN metadata columns for them: {missing_from_metadata}\n'
                'Add these samples to the metadata file, or check that the "sample" column '
                'spelling matches the samples file exactly.'
            )
        if unused_in_metadata:
            self.logger.warning(
                'the following sample(s) in the metadata file have no matching sample '
                f'in the coverage/count output and will not be used: {unused_in_metadata}'
            )

    def _build_coverage_matrix_(self):
        '''
        Build the concatenated, position/coverage-type-labeled coverage matrix -- equivalent
        to what becomes adata.X for the full/default variant -- plus its aligned size-factor
        list. Extracted from create() so it's reusable by compute_variant_contribution() for
        split-variant contributions.
        '''
        x_dfs = []
        size_factors_list = None
        for cov_type in self.cov_types:
            x_df = self._x_build_(cov_type)
            # Build size factors list if it does not exist
            if not size_factors_list:
                size_factors_list = [self.size_factors.get(i) for i in ['_'.join(x.split('_')[1:]) for x in x_df.index.values]]
            # 'adenines', 'thymines', 'cytosines', 'guanines', 'deletions' are already raw counts so they need to be normalized by size factor first
            if cov_type in ['adenines', 'thymines', 'cytosines', 'guanines', 'deletions']:
                x_df = x_df.div(size_factors_list, axis=0)
            x_dfs.append(x_df)
        x_df = pd.concat(x_dfs, axis=1, sort=False)
        x_df = x_df.astype('float64')  # Not sure if this is the dtype I want to use defaults to float64
        # Rename columns of x_df to include position and coverage type
        clist = [[p + '_' + cov for p in self.positions] for cov in self.cov_types]
        clist = [name for per_cov in clist for name in per_cov]
        x_df.columns = clist
        return x_df, size_factors_list

    def _vst_fit_mask_(self, obs_df, minfeaturereads):
        '''
        Marks which rows qualify for the VST dispersion-trend FIT: a tRNA gene's total raw read
        count, summed across all of its samples, must meet `minfeaturereads` (roadmap.md's
        `--mincoverage` -> `--minfeaturereads` item). Sets `obs_df['vst_fit_excluded']` (the
        inverse) so it's visible on the resulting AnnData object, and returns the qualifying mask
        for `_compute_vst_`'s `fit_mask`. This flag affects ONLY `layers['vst']` -- raw counts,
        normalized counts, the obs gene universe, and every counting/DE/volcano/heatmap path are
        governed separately and are unaffected by it.
        '''
        gene_totals = obs_df.groupby('trna')['nreads_total_raw'].transform('sum')
        qualifies = gene_totals >= minfeaturereads
        obs_df['vst_fit_excluded'] = ~qualifies
        return qualifies.values

    def _compute_vst_(self, x_norm, x_raw, sizefactor_values, group_values, vst_strategy, dispfittype, obs_index, var_index, fit_mask=None):
        '''
        Compute a Variance Stabilizing Transformation of `x_norm`/`x_raw`. Extracted from
        create()'s VST block so it's reusable by compute_variant_contribution() for
        split-variant contributions.

        `fit_mask`, if given, is a boolean array (aligned with `obs_index`) restricting which
        rows contribute to the dispersion-trend FIT (roadmap.md's `--minfeaturereads` item) --
        rows outside the mask are excluded from the fit but still get a VST value computed by
        applying the resulting fitted trend to their own (raw, sizefactor-normalized) counts, the
        same as every other row. `None` (or all-True) means "fit on everything", the prior
        behavior.
        '''
        sizefactor_values = np.asarray(sizefactor_values) if sizefactor_values is not None else np.ones(x_raw.shape[0])
        use_fit_mask = fit_mask is not None and not np.all(fit_mask)

        if vst_strategy == 'log1p':
            self.logger.info('Applying Variance Stabilizing Transformation (log1p + StandardScaler)...')
            log_counts = np.log1p(x_norm)
            scaler = sklearn.preprocessing.StandardScaler()
            if use_fit_mask:
                scaler.fit(log_counts[fit_mask])
                return scaler.transform(log_counts)
            return scaler.fit_transform(log_counts)

        self.logger.info('Applying Variance Stabilizing Transformation (PyDESeq2 VST)...')
        try:
            # Prepare raw counts matrix (PyDESeq2 requires integers)
            # Round and cap negatives to 0 (shouldn't be negatives in raw, but safety first)
            raw_integer_counts = np.clip(np.round(x_raw), 0, None).astype(int)

            # Restrict what PyDESeq2 sees to the fit-qualifying rows only -- the dispersion trend
            # curve is fit exclusively from these, so a handful of noisy/low-coverage features
            # can't distort it for every other feature too. The fitted trend is applied to every
            # row (this mask's rows and all others) further below.
            if use_fit_mask:
                fit_counts, fit_obs_index = raw_integer_counts[fit_mask], obs_index[fit_mask]
                fit_sf = sizefactor_values[fit_mask]
                fit_group = np.asarray(group_values)[fit_mask] if not isinstance(group_values, str) else group_values
            else:
                fit_counts, fit_obs_index, fit_sf, fit_group = raw_integer_counts, obs_index, sizefactor_values, group_values

            raw_df = pd.DataFrame(fit_counts, index=fit_obs_index, columns=var_index)

            # Create a minimalist metadata DataFrame
            meta_df = pd.DataFrame({'condition': pd.Series(np.asarray(fit_group), index=fit_obs_index)}, index=fit_obs_index)

            # Initialize PyDESeq2 DeseqDataSet
            # NOTE: Since PyDESeq2's dispersion estimations often use parametric/mean fit
            # based on feature variances, we use the fit_type that was passed previously.
            vst_dds = DeseqDataSet(
                counts=raw_df,
                metadata=meta_df,
                design_factors="condition",
                fit_type=dispfittype,
                quiet=True
            )

            # Attach the pre-calculated (tRNA-control, default) size factors so VST reflects
            # the same scaling as the rest of the pipeline.
            #
            # PyDESeq2's own vst_fit() only skips its internal fit_size_factors() call when
            # BOTH obsm['size_factors'] is set AND self.logmeans is not None (pydeseq2/dds.py,
            # vst_fit). self.logmeans is only ever set inside fit_size_factors() itself, so
            # setting obsm['size_factors'] alone was not enough -- fit_size_factors() silently
            # reran and discarded the size factors set here. Because tRNA coverage data is
            # zero-heavy enough that nearly every feature has at least one zero-count sample,
            # that fallback landed in PyDESeq2's "iterative" size-factor method, which jointly
            # optimizes one size factor per SAMPLE via scipy.optimize.minimize(method="Powell")
            # -- a derivative-free search whose cost blows up non-linearly with sample count,
            # hanging on any dataset with more than roughly 50-100 samples.
            #
            # We avoid both problems by computing logmeans/filtered_genes the same cheap,
            # vectorized way PyDESeq2's own "poscounts" fit_type does (no optimization loop),
            # and deriving normed_counts from our own pre-computed size factors so
            # fit_size_factors() is skipped entirely.
            sf = np.asarray(fit_sf) if fit_sf is not None else np.ones(vst_dds.n_obs)
            nz_log_counts = np.zeros_like(vst_dds.X, dtype=float)
            np.log(vst_dds.X, out=nz_log_counts, where=vst_dds.X != 0)
            logmeans = nz_log_counts.mean(axis=0)
            vst_dds.filtered_genes = (~np.isinf(logmeans)) & (logmeans > 0)
            vst_dds.logmeans = logmeans
            vst_dds.obsm['size_factors'] = sf
            vst_dds.obs['size_factors'] = sf
            vst_dds.layers['normed_counts'] = vst_dds.X / sf[:, None]
            vst_dds.var['_normed_means'] = vst_dds.layers['normed_counts'].mean(axis=0)

            if not use_fit_mask:
                # Calculate vst
                vst_dds.vst(use_design=False)
                return vst_dds.layers['vst_counts']

            # Fit only (no transform yet) -- vst_dds only has the fit-qualifying rows, so its own
            # vst_transform() can't be reused directly for the full row set below: applying it to
            # unseen counts recomputes normalization via PyDESeq2's own median-of-ratios estimate
            # instead of reusing our externally-supplied per-row size factors, which would give
            # excluded rows a different, inconsistent normalization from included ones. Instead,
            # the fitted trend curve/dispersion is read back below and applied by hand to every
            # row's own size-factor-normalized counts, the same formula vst_transform() itself
            # uses internally (pydeseq2/dds.py's vst_transform).
            vst_dds.vst_fit(use_design=False)
            full_normed_counts = raw_integer_counts / sizefactor_values[:, None]
            vst_fit_type = vst_dds.vst_fit_type
            if vst_fit_type == 'parametric':
                a0, a1 = vst_dds.uns['vst_trend_coeffs']
                return np.log2(
                    (1 + a1 + 2 * a0 * full_normed_counts
                     + 2 * np.sqrt(a0 * full_normed_counts * (1 + a1 + a0 * full_normed_counts)))
                    / (4 * a0)
                )
            elif vst_fit_type == 'mean':
                from scipy.stats import trim_mean
                gene_dispersions = vst_dds.var['vst_genewise_dispersions']
                use_for_mean = gene_dispersions > 10 * vst_dds.min_disp
                mean_disp = trim_mean(gene_dispersions[use_for_mean], proportiontocut=0.001)
                return (2 * np.arcsinh(np.sqrt(mean_disp * full_normed_counts))
                        - np.log(mean_disp) - np.log(4)) / np.log(2)
            else:
                raise NotImplementedError(f"Found fit_type '{vst_fit_type}'. Expected 'parametric' or 'mean'.")
        except Exception as e:
            self.logger.warning(f"PyDESeq2 native VST failed ({e}). Falling back to log1p + StandardScaler...")
            log_counts = np.log1p(x_norm)
            scaler = sklearn.preprocessing.StandardScaler()
            if use_fit_mask:
                scaler.fit(log_counts[fit_mask])
                return scaler.transform(log_counts)
            return scaler.fit_transform(log_counts)

    def compute_variant_contribution(self, vst_strategy='vst', dispfittype='parametric', minfeaturereads=30):
        '''
        Compute this loader's data as a "variant contribution" -- coverage matrices, numeric
        per-obs columns, and count/sizefactor summaries -- for merging into an existing target
        AnnData object via merge_variant_into_adata(). Does not build a standalone AnnData.
        Intended to be called on a "loader-only" instance (constructed with analysis_args=None,
        pointed at one split variant's results_dir_name/graphs_dir_name).
        '''
        x_df, size_factors_list = self._build_coverage_matrix_()
        obs_df = self._obs_build_(x_df)

        # Numeric, per-obs columns only -- identity columns (trna, sample, group, amino, ...)
        # are shared across all variants and already live on the target adata's real .obs.
        numeric_cols = [c for c in obs_df.columns if c.startswith('nreads_') or c == 'sizefactor']
        obsm_counts = obs_df[numeric_cols].rename(columns={'sizefactor': 'deseq2_sizefactor'})

        x_norm = x_df
        x_raw = pd.DataFrame(x_norm.values * np.array(size_factors_list)[:, None], index=x_df.index, columns=x_df.columns)
        # Complete-variant only -- see VariantContribution.x_norm_allfeatures.
        x_norm_allfeatures = None
        if self.size_factors_allfeatures is not None:
            allfeature_sf = np.array([self.size_factors_allfeatures.get(s, 1.0) for s in obs_df['sample'].values])
            x_norm_allfeatures = pd.DataFrame(x_raw.values / allfeature_sf[:, None], index=x_df.index, columns=x_df.columns)

        x_vst = None
        if vst_strategy != 'none':
            fit_mask = self._vst_fit_mask_(obs_df, minfeaturereads)
            x_vst = self._compute_vst_(
                x_norm.values, x_raw.values, obsm_counts['deseq2_sizefactor'].values, obs_df.get('group', 'all'),
                vst_strategy, dispfittype, x_df.index, x_df.columns, fit_mask=fit_mask,
            )

        return VariantContribution(
            x_raw=x_raw, x_norm=x_norm, x_norm_allfeatures=x_norm_allfeatures, x_vst=x_vst,
            obsm_counts=obsm_counts, sizefactors_trna=self.size_factors, sizefactors_allfeatures=self.size_factors_allfeatures,
            type_counts=self.type_counts, type_real_counts=self.type_real_counts, amino_counts=self.amino_counts,
            anticodon_counts=self.anticodon_counts, nontrna_counts=self.non_trna_read_counts,
            mismatch_counts=self.mismatch_counts,
        )

    def _apply_readlength_split_(self, adata):
        '''
        Compute and merge the under/over read-length split variants for `self.analysis_args
        .readlengthsplit` into `adata` in place, as new layers/obsm/uns entries (see
        merge_variant_into_adata()) -- replaces the old behavior of writing separate
        `_u{N}.h5ad`/`_o{N}.h5ad` files. On-disk `results/u{N}`/`graphs/u{N}` (and o{N})
        directories are still produced via AnalysisPipeline, unchanged. The split BAM files
        themselves (under `<bamdir>/u{N}`/`o{N}`) are temporary scratch files by default and
        are deleted once this variant has been merged into `adata` -- pass `--savesplitbams`
        to keep them on disk instead.
        '''
        cutoff = self.analysis_args.readlengthsplit
        abs_output = os.path.abspath(self.output)
        base_output_dir = os.path.dirname(abs_output)
        default_bamdir = self.analysis_args.bamdir if self.analysis_args.bamdir else os.path.join("processed", "bam")
        vst_strategy = str(getattr(self.analysis_args, 'vst', 'vst')).lower()
        dispfittype = getattr(self.analysis_args, 'dispfittype', 'parametric')
        minfeaturereads = getattr(self.analysis_args, 'minfeaturereads', None)
        minfeaturereads = int(minfeaturereads) if minfeaturereads is not None else 30
        savesplitbams = getattr(self.analysis_args, 'savesplitbams', False)

        try:
            for direction, prefix in [('under', 'u'), ('over', 'o')]:
                tag = f'{prefix}{cutoff}'
                self.logger.info(f"Running analysis pipeline ({direction.capitalize()} {cutoff})...")

                args_variant = SimpleNamespace(**vars(self.analysis_args))
                args_variant.bamdir = os.path.join(default_bamdir, tag)
                # CRITICAL: Prevent recursive splitting by nullifying readlengthsplit
                args_variant.readlengthsplit = None
                args_variant.results_dir_name, args_variant.graphs_dir_name = toolsTG.variant_dir_names(args_variant, tag=tag)
                args_variant.output = os.path.join(base_output_dir, f"{os.path.splitext(os.path.basename(abs_output))[0]}_{tag}.h5ad")

                pipeline_variant = AnalysisPipeline(
                    args_variant, expname=base_output_dir,
                    phase_tracker=self.phase_tracker, variant_label=f"{direction.capitalize()} {cutoff}",
                    split_tag=tag,
                )
                pipeline_variant.run()

                self.logger.info(f"Building AnnData contribution ({direction.capitalize()} {cutoff})...")
                loader = AnnDataBuilder(base_output_dir, self.metadata_path, None, analysis_args=None,
                                         results_dir_name=args_variant.results_dir_name, graphs_dir_name=args_variant.graphs_dir_name)
                contribution = loader.compute_variant_contribution(vst_strategy=vst_strategy, dispfittype=dispfittype, minfeaturereads=minfeaturereads)

                build_flags = toolsTG.sanitize_flags(args_variant)
                merge_variant_into_adata(adata, contribution, tag=tag, direction=direction, cutoff=cutoff, build_flags=build_flags, overwrite=True)
        finally:
            if not savesplitbams:
                for tag in (f'u{cutoff}', f'o{cutoff}'):
                    if self._split_dirs_preexisted.get(tag, False):
                        # This tag's directory already had every sample's split BAM before this
                        # run started, and this run didn't force a fresh split -- it's not this
                        # run's scratch to clean up, regardless of --savesplitbams here.
                        continue
                    split_dir = os.path.join(default_bamdir, tag)
                    if os.path.isdir(split_dir):
                        shutil.rmtree(split_dir, ignore_errors=True)

    def _seq_build_(self, gap=False):
        # Build reference sequence dataframe
        seq_df = self._x_build_('actualbase')
        if not gap:
            # Drop gap positions from reference sequence dataframe
            seq_df = seq_df.loc[:,~seq_df.columns.str.contains('a')]
            seq_df = seq_df.loc[:,~seq_df.columns.str.contains('b')]
            seq_df = seq_df.loc[:,~seq_df.columns.str.contains('e')]
            seq_df = seq_df.loc[:,~seq_df.columns.str.contains('-1')]
        seqs = [''.join(x) for x in seq_df.values.tolist()]

        return seqs

    def _x_build_(self, cov_type):
        '''
        Build x dataframe from coverage file
        '''
        x_df = self.coverage.pivot(index='uniquefeat', values=[cov_type], columns='position')
        cols = x_df.columns.get_level_values(1).values
        x_df = x_df.T.reset_index(drop=True).T
        x_df.columns = cols
        x_df = x_df.reindex(columns=self.positions)
        
        return x_df
    
    def _obs_build_(self, x_df):
        '''
        Build obs dataframe from coverage file derived x dataframe
        '''
        # Create obs dataframe from coverage file derived x dataframe
        obs_df = pd.DataFrame([[x.split('_')[0], '_'.join(x.split('_')[1:])] for x in x_df.index.values], columns=['trna','sample'], index=x_df.index)
        # Add metadata to obs dataframe
        for i in self.metadata:
            obs_df[i] = [self.metadata.get(i).get(x) for x in obs_df['sample'].values]
        # Add tRNA type, amino acid, iso, refseq, sizefactor, and nreads to obs dataframe
        trna_obs = obs_df['trna'].str.split('-',n=4,expand=True)
        obs_df['pseudogene'] = trna_obs[0]
        obs_df['amino'] = trna_obs[1]
        obs_df['iso'] = trna_obs[2]
        obs_df['refseq'] = self.seqs
        obs_df['refseq_full'] = self.seqs_full
        # Build size factors list, warning about missing samples and defaulting to 1.0
        sizefactor_list = []
        missing_samples = set()
        for x in x_df.index.values:
            sample_key = '_'.join(x.split('_')[1:])
            sf = self.size_factors.get(sample_key)
            if sf is None:
                missing_samples.add(sample_key)
                sf = 1.0  # Default to 1.0 when size factor is not found
            sizefactor_list.append(sf)
        
        if missing_samples:
            self.logger.warning(f"Size factors not found for {len(missing_samples)} samples. Using default value of 1.0.")
            self.logger.warning(f"  Missing samples: {list(missing_samples)[:5]}{'...' if len(missing_samples) > 5 else ''}")
        
        obs_df['sizefactor'] = sizefactor_list
        # obs_df['nreads_unique_raw'] = [self.unique_counts.get(i[0]).get('_'.join(i[1:])) if self.unique_counts.get(i[0]) else 0 for i in [x.split('_') for x in x_df.index.values]] # Some samples may not have any reads for a given tRNA in the unique_counts dictionary might want to double check this
        # Create unique counts for each tRNA split by type
        for rt in ['fiveprime','threeprime','wholecounts','other']:
            obs_df[f'nreads_{rt}_unique_raw'] = [self.unique_counts.get(i[0]+f'_{rt}').get('_'.join(i[1:])) if self.unique_counts.get(i[0]+f'_{rt}') else 0 for i in [x.split('_') for x in x_df.index.values]]
            obs_df[f'nreads_{rt}_unique_norm'] = obs_df[f'nreads_{rt}_unique_raw']/obs_df['sizefactor']
        # Add total unique read counts to obs dataframe
        obs_df['nreads_total_unique_raw'] = obs_df[[f'nreads_{rt}_unique_raw' for rt in ['fiveprime','threeprime','wholecounts','other']]].sum(axis=1)
        obs_df['nreads_total_unique_norm'] = obs_df['nreads_total_unique_raw']/obs_df['sizefactor']
        # Add read counts (from bowtie2) for each tRNA split by type
        for rt in self.read_types:
            obs_df['nreads_' + rt + '_norm'] = [self.normalized_read_counts.get(j[0] + '_' + rt).get('_'.join(j[1:])) if self.normalized_read_counts.get(j[0] + '_' + rt) else 0 for j in [x.split('_') for x in x_df.index.values]]
            obs_df['nreads_' + rt + '_raw'] = obs_df['nreads_' + rt + '_norm']*obs_df['sizefactor']
        # Add total read counts (from bowtie2) to obs dataframe - excluding partial precounts, trailer counts, and antisense counts so that it matches unique counts
        obs_df['nreads_total_raw'] = obs_df[[f'nreads_{rt}_raw' for rt in ['fiveprime','threeprime','wholecounts','other']]].sum(axis=1)
        obs_df['nreads_total_norm'] = obs_df['nreads_total_raw']/obs_df['sizefactor']
        # Add unique feature column to obs dataframe - this is the index of the x dataframe - not sure if this is necessary
        obs_df['uniquefeat'] = obs_df.index.values
        return obs_df

    def _adata_build_(self, obs_df, x_df):
        '''
        Build AnnData object from obs and x dataframes
        '''
        positions = [i.split('_')[0] for i in x_df.columns.values]
        coverage = [i.split('_')[-1] for i in x_df.columns.values]
        # Generate gap list
        gap_list = pd.Series(positions)
        gap_list = gap_list[gap_list.str.contains('a') | gap_list.str.contains('b') | gap_list.str.contains('e') | gap_list.str.contains('-1')]
        gap_list = pd.Series(positions).isin(gap_list).tolist()
        # Build AnnData object
        adata = ad.AnnData(x_df)
        adata.var['gap'] = gap_list
        adata.var['positions'] = positions
        adata.var['coverage'] = coverage
        # Create sprinzl position and fragment type information
        loc_dict, loc_half_dict, loc_code_dict = _sprinzl_location_maps()
        # Add the location data to adata object
        adata.var['location'] = adata.var['positions'].map(loc_dict)
        adata.var['half'] = adata.var['positions'].map(loc_half_dict)
        # Short region codes from SprinzlPos.pm, alongside the descriptive names.
        adata.var['location_code'] = adata.var['positions'].map(loc_code_dict)
        # Add metadata dataframe
        adata.obs['dataset'] = os.path.basename(self.output) # Name of the dataset output file - Usefull for combining multiple datasets if the merge function is used
        adata.obs['trna'] = obs_df['trna'].values
        adata.obs['iso'] = obs_df['iso'].values
        adata.obs['amino'] = obs_df['amino'].values
        # Add sample and group metadata
        adata.obs['sample'] = obs_df['sample'].values
        adata.obs['group'] = obs_df['group'].values
        adata.obs['pseudogene'] = obs_df['pseudogene'].values
        # Add custom dataframe obs
        for ob in self.observations:
            adata.obs[ob] = obs_df[ob].values
        # Add the numer of reads per tRNA as observations-annotation to adata from trna unique counts file
        for rt in ['wholecounts','fiveprime','threeprime','other']:
            adata.obs['nreads_' + rt + '_unique_raw'] = obs_df['nreads_' + rt + '_unique_raw'].values
            adata.obs['nreads_' + rt + '_unique_norm'] = obs_df['nreads_' + rt + '_unique_norm'].values
        adata.obs['nreads_total_unique_raw'] = obs_df['nreads_total_unique_raw'].values
        adata.obs['nreads_total_unique_norm'] = obs_df['nreads_total_unique_norm'].values
        # Add the numer of reads per tRNA as observations-annotation to adata from trna normalized read counts file
        for rt in self.read_types:
            adata.obs['nreads_' + rt + '_raw'] = obs_df['nreads_' + rt + '_raw'].values
            adata.obs['nreads_' + rt + '_norm'] = obs_df['nreads_' + rt + '_norm'].values
        adata.obs['nreads_total_raw'] = obs_df['nreads_total_raw'].values
        adata.obs['nreads_total_norm'] = obs_df['nreads_total_norm'].values
        # Add fragment type classification to the adata object
        fragtype = []
        for i in range(adata.shape[0]):
            fiveprimereadends = adata.obs['nreads_fiveprime_unique_norm'].iloc[i]
            fiveprimemean = adata.X[i][adata.var['half']=='fiveprime'].mean()
            fiveprimestd = np.std(adata.X[i][adata.var['half']=='fiveprime'])
            threeprimereadends = adata.obs['nreads_threeprime_unique_norm'].iloc[i]
            threeprimemean = adata.X[i][adata.var['half']=='threeprime'].mean()
            threeprimestd = np.std(adata.X[i][adata.var['half']=='threeprime'])
            totalreads = adata.obs['nreads_total_unique_norm'].iloc[i]
            totalstd = np.std(adata.X[i])
            centermean = adata.X[i][adata.var['half']=='center'].mean()
            centerstd = np.std(adata.X[i][adata.var['half']=='center'])
            # If the readend counts are more than 2 standard deviations away from the opposite end then it is a fragment
            if totalreads <= 20:
                fragtype.append('low_coverage')
            elif fiveprimereadends > np.abs(threeprimereadends + 2*totalstd):
                if np.abs(centermean - adata.X[i][adata.var['location'].isin(FRAGMENT_FIVEPRIME_REGIONS)].mean()) < fiveprimestd:
                    fragtype.append('fiveprime_half')
                else:
                    fragtype.append('fiveprime_fragment')
            elif threeprimereadends > np.abs(fiveprimereadends + 2*totalstd):
                if np.abs(centermean - adata.X[i][adata.var['location'].isin(FRAGMENT_THREEPRIME_REGIONS)].mean()) < threeprimestd:
                    fragtype.append('threeprime_half')
                else:
                    fragtype.append('threeprime_fragment')
            else:
                if np.abs(fiveprimemean - threeprimemean) > totalstd:
                    fragtype.append('other_fragment')
                elif np.abs(adata.X[i][adata.var['location'].isin(FRAGMENT_FIVEPRIME_REGIONS)].mean() - \
                            adata.X[i][adata.var['location'].isin(FRAGMENT_THREEPRIME_REGIONS)].mean()) > centermean + centerstd:
                    fragtype.append('multiple_fragment')
                else:
                    fragtype.append('whole')
        adata.obs['fragment'] = fragtype
        # Add size factor
        adata.obs['deseq2_sizefactor'] = obs_df['sizefactor'].values
        # Add aligned reference sequence based on values in coverage file
        adata.obs['refseq'] = obs_df['refseq'].values
        adata.obs['refseq_full'] = obs_df['refseq_full'].values
        # Add anticodon counts as uns
        adata.uns['anticodon_counts'] = self.anticodon_counts
        # Add amino acid counts as uns
        adata.uns['amino_counts'] = self.amino_counts
        # Add type counts as uns
        adata.uns['type_counts'] = self.type_counts
        adata.uns['type_real_counts'] = self.type_real_counts
        # Add the read-level mismatch histogram as uns (absent if the count file wasn't written)
        if self.mismatch_counts is not None:
            adata.uns['mismatch_counts'] = self.mismatch_counts
        # Add non tRNA counts as uns
        adata.uns['nontRNA_counts'] = self.non_trna_read_counts
        # Add runinfo as uns
        adata.uns['trnagraphruninfo'] = self.trnagraph_run_info
        # Add both DESeq2 size factor sets as uns (tRNA-controlled default, and all-feature secondary)
        adata.uns['deseq2_sizefactors_trna'] = self.size_factors
        adata.uns['deseq2_sizefactors_allfeatures'] = self.size_factors_allfeatures
        # Replicate-correlation QC. Cheap (one correlation matrix over the count table) and the
        # fastest way to see whether the declared grouping is supported by the data, or whether a
        # single sample disagrees with its own replicates badly enough to warrant dropping.
        try:
            corr = toolsTG.replicate_correlation(adata.obs)
            adata.uns['replicate_correlation'] = corr['per_sample']
            adata.uns['replicate_correlation_pairs'] = corr['pairs']
            adata.uns['replicate_correlation_summary'] = corr['summary']
            self._log_replicate_correlation(corr)
        except Exception as exc:
            # Diagnostic only -- never fail a build over it.
            self.logger.warning(f"Could not compute replicate correlation: {exc}")

        # Add 'group' log2FC value/pval to uns since it is the default for the volcano/heatmap and saves time later
        _precompute_default_log2fc(adata, threads=self.analysis_args.threads)

        return adata

    def _log_replicate_correlation(self, corr):
        '''Reports the QC to the log and to a tab-separated file beside the other results.'''
        summary = corr['summary']
        if summary['n_within_pairs'] == 0:
            self.logger.warning(
                "Replicate correlation: every sample is its own group, so there are no "
                "replicates to correlate. If that is not intended, check the 'group' column of "
                "your metadata -- the generated trim_metadata.tsv template sets group = sample."
            )
        else:
            self.logger.info(
                f"Replicate correlation: within-group r2 {summary['within_mean']:.4f}, "
                f"between-group r2 {summary['between_mean']:.4f}, "
                f"separation {summary['separation']:+.4f}"
            )

        path = self.expinfo.replicatecorrelation
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as out:
            for key, value in summary.items():
                formatted = f"{value:.6f}" if isinstance(value, float) else value
                print(f"# {key}\t{formatted}", file=out)
            corr['per_sample'].to_csv(out, sep='\t', index=False, float_format='%.6f')

    def _handle_preprocessing(self, args):
        '''
        Handles auto-mapping and splitting based on args.
        '''
        bamdir = args.bamdir if args.bamdir else os.path.join("processed", "bam")
        # Two operations, two flags. They were one, which made "re-cut the splits" impossible to
        # ask for without also re-mapping from FASTQ -- unreachable once the FASTQs have moved,
        # which is exactly the situation a stale split is discovered in.
        overwrite = getattr(args, 'overwritebams', False)
        overwrite_splits = getattr(args, 'overwritesplits', False)
        
        # Check if we need to map
        try:
            sf = toolsTG.samplefile(args.input)
            samples = sf.getsamples()
        except Exception as e:
            self.logger.error(f"Could not parse metadata file '{args.input}' to check for existing BAMs: {e}")
            raise

        missing_bams = False
        if not os.path.exists(bamdir):
            missing_bams = True
        else:
            for s in samples:
                if not os.path.isfile(os.path.join(bamdir, f"{s}.bam")):
                    missing_bams = True
                    break
        
        if missing_bams or overwrite:
            self.logger.info("Running mapping step...")
            map_args = SimpleNamespace(**vars(args))
            map_args.output = os.path.basename(os.path.dirname(args.output)) # Experiment name
            map_args.force_remap = overwrite
            map_args.local = False 
            map_args.skipcheck = False
            
            from . import toolsMap
            toolsMap.MapSamples(map_args).main()
        else:
            self.logger.info("BAM files found. Skipping mapping (use --overwritebams to force).")

        # Splitting Logic
        if hasattr(args, 'readlengthsplit') and args.readlengthsplit:
            cutoff = args.readlengthsplit
            # Check if split BAMs exist -- recorded *before* this run's own splitting step runs,
            # so _apply_readlength_split_()'s later cleanup can tell "already there, untouched by
            # this run" apart from "this run just created (or --overwritesplits-regenerated) it"
            # (see roadmap.md's "BAM deletion safety" item). --overwritesplits forces a fresh
            # split even for an already-complete tag, so that tag is this run's own output
            # regardless.
            preexisting = _split_bam_dirs_preexisting(bamdir, cutoff, samples)
            self._split_dirs_preexisted = {tag: complete and not overwrite_splits
                                           for tag, complete in preexisting.items()}
            missing_split = not all(preexisting.values())

            # Before deciding to REUSE a cached split, check it still belongs to the BAMs beside
            # it. Skipped when re-cutting anyway, since the offending files are about to go.
            if not (missing_split or overwrite_splits):
                _refuse_stale_split_bams(bamdir, cutoff, samples, self.logger)

            if missing_split or overwrite_splits:
                self.logger.info(f"Running split step (cutoff {cutoff})...")
                from . import toolsSplit
                split_args = SimpleNamespace(
                    input=args.input,
                    readlengthsplit=cutoff,
                    bamdir=bamdir,
                    overwritebams=overwrite_splits,
                    threads=args.threads
                )
                toolsSplit.BamSplitter(split_args).process()
            else:
                self.logger.info(f"Split BAM files found for cutoff {cutoff}, and each is newer "
                                 f"than the BAM it was cut from. Skipping split (use "
                                 f"--overwritesplits to force).")


def _get_git_version_():
    '''
    Return (git_version, git_hash) for the installed trnagraph package, with a fallback for
    non-git installations. Run from this file's own directory, letting git find the .git
    folder in a parent directory.
    '''
    trnagraphdir = os.path.dirname(os.path.abspath(__file__))
    try:
        git_version = subprocess.check_output(
            ['git', 'describe', '--always'], cwd=trnagraphdir, stderr=subprocess.DEVNULL
        ).decode().strip()
        git_hash = subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=trnagraphdir, stderr=subprocess.DEVNULL
        ).decode().strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        git_version = 'unknown (not a git repository)'
        git_hash = 'unknown'
    return git_version, git_hash


def _split_bam_dirs_preexisting(bamdir, cutoff, samples):
    '''
    For each of u{cutoff}/o{cutoff}, checks whether every sample's split BAM was already present
    under `bamdir` *before* this run's own splitting step runs. Used to scope split-BAM cleanup
    (`--savesplitbams` off) to directories this run actually created -- a directory that already
    had every sample's split BAM present (e.g. kept on disk by a prior `--savesplitbams` run) is
    never deleted, regardless of whether this run's own splitting step ran or was skipped.
    '''
    preexisting = {}
    for tag in (f'u{cutoff}', f'o{cutoff}'):
        tag_dir = os.path.join(bamdir, tag)
        preexisting[tag] = os.path.isdir(tag_dir) and all(
            os.path.isfile(os.path.join(tag_dir, f"{s}.bam")) for s in samples
        )
    return preexisting


def _stale_split_bams(bamdir, cutoff, samples):
    '''
    Split BAMs that are older than the unsplit BAM they were supposedly cut from.

    A split BAM is written by reading its parent, so it can only ever be newer. Older means the
    parent has been replaced since -- and the split now describes reads that are no longer in it.

    This is not hypothetical. On the hg38 c305 dataset, `preprocess map --dedup` wrote the
    deduplicated BAMs into `processed/bam_dedup/` on top of the pre-dedup ones, leaving the
    `u60/`/`o60/` subdirectories cut from the PRE-dedup parents in place beside them. The next
    build found every split BAM present, logged "Split BAM files found for cutoff 60. Skipping
    split", and counted pre-dedup reads into both split variants: `u60 + o60` came to 68.7M
    reads against the full variant's 6.97M, a tenfold overstatement that looked entirely
    plausible in every figure downstream. The giveaway was a `u60` split BAM of 1.44 GB whose
    parent was 621 MB -- a length-filtered subset larger than the whole.

    Returns [(tag, sample, split_path)] for every offender, so the caller can name them.
    '''
    stale = []
    for tag in (f'u{cutoff}', f'o{cutoff}'):
        tag_dir = os.path.join(bamdir, tag)
        if not os.path.isdir(tag_dir):
            continue
        for sample in samples:
            split_path = os.path.join(tag_dir, f"{sample}.bam")
            parent_path = os.path.join(bamdir, f"{sample}.bam")
            if not (os.path.isfile(split_path) and os.path.isfile(parent_path)):
                continue
            if os.path.getmtime(split_path) < os.path.getmtime(parent_path):
                stale.append((tag, sample, split_path))
    return stale


def _refuse_stale_split_bams(bamdir, cutoff, samples, logger):
    '''
    Abort, naming the offending files, when a cached split predates its parent.

    Refuses rather than silently re-splitting: re-splitting rewrites gigabytes of the user's
    BAMs off the back of an mtime comparison, and -- more importantly -- a silent fix would have
    produced correct numbers here while leaving the user believing their splits had been sound
    all along. The staleness is itself the finding worth surfacing.
    '''
    stale = _stale_split_bams(bamdir, cutoff, samples)
    if not stale:
        return
    shown = '\n'.join(f'  {path}' for _, _, path in stale[:8])
    if len(stale) > 8:
        shown += f'\n  ... and {len(stale) - 8} more'
    raise toolsTG.InvalidParameterError(
        f"{len(stale)} split BAM file(s) for cutoff {cutoff} are OLDER than the unsplit BAM they "
        f"were cut from, so they describe reads that are no longer in it -- most often because "
        f"the parent BAMs were re-mapped or deduplicated after the split was taken:\n{shown}\n"
        f"Counting these would silently mix the two read sets. Re-cut them with "
        f"--overwritesplits, or delete {os.path.join(bamdir, f'u{cutoff}')} and "
        f"{os.path.join(bamdir, f'o{cutoff}')} and run again."
    )


def _empty_split_nontrna_counts(target_adata):
    '''
    The non-tRNA counts a split variant carries: an empty frame over the object's sample axis.

    A read-length cutoff partitions tRNA reads by design, but non-tRNA features are not being
    classified by that criterion at all and span a far wider size range, so a split's non-tRNA
    numbers record where the cutoff fell rather than anything about the data. They are therefore
    excluded from split variants outright rather than recomputed per split (which produced the
    "o60 got the non-tRNA plots, u60 didn't" bug) or copied from the full variant (which left a
    duplicate that nothing was allowed to read).

    This must be an explicit empty value rather than an omission: build_variant_view() overlays
    uns keys onto a split's view only when the key is present, so leaving it out would show the
    parent object's full-variant counts on the split instead.

    Columns come from the full variant's own frame where it has any, preserving its exact sample
    labels and order, and otherwise from the object's obs -- `--gtf` is optional, and without one
    the full variant's frame is itself empty and may carry no columns at all.
    '''
    existing = target_adata.uns.get('nontRNA_counts')
    if existing is not None and len(existing.columns) > 0:
        columns = list(existing.columns)
    else:
        columns = list(dict.fromkeys(target_adata.obs['sample'].astype(str)))
    return pd.DataFrame(columns=columns)


def merge_variant_into_adata(target_adata, contribution: VariantContribution, tag, direction, cutoff, build_flags, overwrite=False):
    '''
    Merge a VariantContribution (from AnnDataBuilder.compute_variant_contribution()) into
    `target_adata` as a new size-split variant, under the naming scheme documented in
    tRNAgraph/docs/data_structure.md:
      - layers[f'raw_{tag}'] / layers[f'norm_{tag}'] / layers[f'vst_{tag}'] (deliberately no
        norm_allfeatures_{tag}: all-feature normalization is complete-variant only)
      - obsm[f'size_split_{tag}'] (all numeric per-obs columns, unsuffixed names)
      - uns['size_splits'][tag] (sizefactors, counts, log2FC, build provenance)
    Mutates and returns `target_adata` in place; does not write it to disk -- the caller
    decides when/where to persist. `target_adata.obs.index` must already be in its final,
    prefixed form (i.e. this must be called AFTER AnnDataBuilder.create()'s obs-index-prefix
    step) so the prefix used to align `contribution`'s raw `{trna}_{sample}` index matches.
    '''
    existing = target_adata.uns.get('size_splits', {})
    if tag in existing and not overwrite:
        raise ValueError(f"Split variant '{tag}' already exists in this AnnData object. Pass --overwrite to replace it.")

    # contribution frames are indexed by the raw '{trna}_{sample}' form (no dataset-basename
    # prefix, since compute_variant_contribution() never calls _adata_build_()); reconstruct
    # the same prefix the target's own obs index carries so rows align.
    dataset_basename = str(target_adata.obs['dataset'].iloc[0]).split('.')[0]
    prefixed_index = pd.Index([f'{dataset_basename}_{raw_idx}' for raw_idx in contribution.obsm_counts.index])

    def _reindexed(df):
        df = df.copy()
        df.index = prefixed_index
        aligned = df.reindex(target_adata.obs.index)
        missing = aligned.index[aligned.isna().all(axis=1)]
        if len(missing) > 0:
            logger.warning(f"{len(missing)} observations in the target object had no matching row in the '{tag}' split contribution (e.g. {list(missing[:3])}); filling with 0.")
        return aligned.fillna(0.0)

    target_adata.layers[f'raw_{tag}'] = _reindexed(contribution.x_raw).values
    target_adata.layers[f'norm_{tag}'] = _reindexed(contribution.x_norm).values
    # No norm_allfeatures_{tag} layer: all-feature normalization is complete-variant only, so
    # `--variant allfeatures:<tag>` deliberately has nothing to resolve to.
    if contribution.x_vst is not None:
        x_vst_df = pd.DataFrame(contribution.x_vst, index=contribution.x_norm.index, columns=contribution.x_norm.columns)
        target_adata.layers[f'vst_{tag}'] = _reindexed(x_vst_df).values

    target_adata.obsm[f'size_split_{tag}'] = _reindexed(contribution.obsm_counts)

    git_version, git_hash = _get_git_version_()
    target_adata.uns.setdefault('size_splits', {})[tag] = {
        'cutoff': cutoff,
        'direction': direction,
        'date_added': datetime.datetime.now().isoformat(),
        'trnagraph_git_version': git_version,
        'trnagraph_git_hash': git_hash,
        'results_dir_name': build_flags.get('results_dir_name'),
        'graphs_dir_name': build_flags.get('graphs_dir_name'),
        'build_flags': build_flags,
        'sizefactors_trna': contribution.sizefactors_trna,
        'type_counts': contribution.type_counts,
        'type_real_counts': contribution.type_real_counts,
        'amino_counts': contribution.amino_counts,
        'anticodon_counts': contribution.anticodon_counts,
        'nontRNA_counts': _empty_split_nontrna_counts(target_adata),
        # tRNA rows only: a split variant has its non-tRNA features excluded, so keeping the
        # 'nontrna' rows here would let the mismatch plot draw a series nothing else in the
        # variant reports.
        'mismatch_counts': _trna_only_mismatch_counts(contribution.mismatch_counts),
    }

    # Precompute default log2FC for common cutoffs, mirroring _adata_build_'s equivalent block
    # for the full variant, via a temporary resolved view so adataLog2FC's obs-column lookups
    # work completely unchanged.
    temp_spec = toolsTG.VariantTag(raw=f'norm:{tag}', norm='norm', tag=tag)
    temp_view = toolsTG.build_variant_view(target_adata, temp_spec)
    variant_threads = build_flags.get('threads')
    _precompute_default_log2fc(temp_view, threads=variant_threads if isinstance(variant_threads, int) else None)
    target_adata.uns['size_splits'][tag]['log2FC'] = temp_view.uns.get('log2FC', {})

    return target_adata


def apply_order(args):
    '''
    Write a declared category order into an object that already exists.

    Sibling of add_split(): an explicit, in-place mutation of a built object, so that an object
    predating the declaration can gain it without a rebuild -- which for the human dataset is
    not practical, its inputs being treated as read-only.

    Writes back over the input unless args.output names somewhere else.
    '''
    logger = logging.getLogger(__name__)
    adata = ad.read_h5ad(args.anndata)
    toolsTG.apply_category_order(adata.obs, args.order)
    for column, levels in (args.order or {}).items():
        logger.info(f'Ordered obs column {column!r}: {list(levels)} (reference level: {levels[0]!r})')
    output_path = args.output or args.anndata
    adata.write_h5ad(output_path)
    logger.info(f'Wrote {output_path}')
    return output_path


def add_split(args):
    '''
    Add a new read-length split variant (an under/over cutoff pair) to an EXISTING h5ad
    AnnData object -- e.g. add u50/o50 to an object that already has u60/o60 -- without
    disturbing any variant already present. Implements `trnagraph analyze addsplit`. Uses
    the same compute_variant_contribution()/merge_variant_into_adata() unit that
    AnnDataBuilder._apply_readlength_split_() uses at initial build time, so both paths
    produce identical results for the same cutoff/data. The split BAM files this generates
    (under `<bamdir>/u<N>`/`o<N>`) are temporary scratch files by default and are deleted
    once merged into `adata` -- pass `--savesplitbams` to keep them on disk instead.
    '''
    adata = ad.read_h5ad(args.anndata)

    # Recover original build parameters from this object's own provenance record, reversing
    # the None -> 'None' sentinel used when it was written (see __init__'s run_flags).
    recorded_flags = dict(adata.uns.get('trnagraphruninfo', {}).get('flags', {}))
    recorded_flags = {k: (None if v == 'None' else v) for k, v in recorded_flags.items()}

    def _effective(cli_value, flag_key, default=None):
        value = cli_value if cli_value is not None else recorded_flags.get(flag_key)
        return value if value is not None else default

    effective_input = _effective(args.metadata, 'input')
    effective_bamdir = _effective(args.bamdir, 'bamdir')
    effective_database = _effective(args.database, 'database')
    effective_gtf = _effective(args.gtf, 'gtf')
    effective_dispfittype = _effective(args.dispfittype, 'dispfittype', 'parametric')
    effective_vst = _effective(args.vst, 'vst', 'vst')
    effective_minfeaturereads = int(_effective(getattr(args, 'minfeaturereads', None), 'minfeaturereads', 30))

    if not effective_input or not os.path.isfile(effective_input):
        raise ValueError(f"Could not recover a valid metadata file from this object's build provenance. Pass --metadata explicitly (got: {effective_input!r}).")
    if not effective_bamdir or not os.path.isdir(effective_bamdir):
        raise ValueError(f"Could not recover a valid BAM directory from this object's build provenance. Pass --bamdir explicitly (got: {effective_bamdir!r}).")

    # Light conflicting-run-info validation: refuse (unless --force) if an explicitly-overridden
    # database/gtf differs from what this object was originally built with -- the narrowly-scoped
    # seed for the roadmap's "prevent merging AnnData objects with conflicting run info" item.
    conflicts = []
    for cli_value, flag_key, label in [(args.database, 'database', 'database'), (args.gtf, 'gtf', 'gtf')]:
        recorded_value = recorded_flags.get(flag_key)
        if cli_value is not None and recorded_value is not None and cli_value != recorded_value:
            conflicts.append(f"--{label} '{cli_value}' was given, but this object was originally built with {label}='{recorded_value}'")
    if conflicts:
        message = "Detected parameters that conflict with this object's original build provenance:\n  " + "\n  ".join(conflicts)
        if not args.force:
            raise ValueError(message + "\nPass --force to proceed anyway.")
        logger.warning(f"{message}\nProceeding anyway due to --force.")

    cutoff = args.readlengthsplit
    existing = adata.uns.get('size_splits', {})
    for tag in (f'u{cutoff}', f'o{cutoff}'):
        if tag in existing and not args.overwrite:
            raise ValueError(f"Split variant '{tag}' already exists in this AnnData object. Pass --overwrite to replace it.")

    base_output_dir = os.path.dirname(os.path.abspath(args.anndata))
    savesplitbams = getattr(args, 'savesplitbams', False)

    # Recorded *before* this run's own splitting step runs, so the cleanup below can tell
    # "already there, untouched by this run" apart from "this run just created (or
    # --overwritesplits-regenerated) it" (see roadmap.md's "BAM deletion safety" item).
    try:
        samples = toolsTG.samplefile(effective_input).getsamples()
    except Exception as e:
        logger.error(f"Could not parse metadata file '{effective_input}' to check for existing split BAMs: {e}")
        raise
    preexisting = _split_bam_dirs_preexisting(effective_bamdir, cutoff, samples)
    overwrite_splits = getattr(args, 'overwritesplits', False)
    split_dirs_preexisted = {tag: complete and not overwrite_splits
                             for tag, complete in preexisting.items()}
    if not overwrite_splits:
        _refuse_stale_split_bams(effective_bamdir, cutoff, samples, logger)

    from . import toolsSplit
    split_args = SimpleNamespace(input=effective_input, readlengthsplit=cutoff, bamdir=effective_bamdir,
                                  overwritebams=overwrite_splits, threads=args.threads)
    toolsSplit.BamSplitter(split_args).process()

    try:
        for direction, prefix in [('under', 'u'), ('over', 'o')]:
            tag = f'{prefix}{cutoff}'
            logger.info(f"Running analysis pipeline ({direction.capitalize()} {cutoff})...")

            args_variant = SimpleNamespace(**recorded_flags)
            args_variant.input = effective_input
            args_variant.database = effective_database
            args_variant.gtf = effective_gtf
            args_variant.dispfittype = effective_dispfittype
            args_variant.vst = effective_vst
            args_variant.bamdir = os.path.join(effective_bamdir, tag)
            args_variant.readlengthsplit = None
            args_variant.overwritebams = overwrite_splits
            args_variant.threads = args.threads
            args_variant.results_dir_name, args_variant.graphs_dir_name = toolsTG.variant_dir_names(args_variant, tag=tag)
            args_variant.output = os.path.join(base_output_dir, f"{os.path.splitext(os.path.basename(args.anndata))[0]}_{tag}.h5ad")

            pipeline_variant = AnalysisPipeline(args_variant, expname=base_output_dir, split_tag=tag)
            pipeline_variant.run()

            logger.info(f"Building AnnData contribution ({direction.capitalize()} {cutoff})...")
            loader = AnnDataBuilder(base_output_dir, effective_input, None, analysis_args=None,
                                     results_dir_name=args_variant.results_dir_name, graphs_dir_name=args_variant.graphs_dir_name)
            contribution = loader.compute_variant_contribution(vst_strategy=str(effective_vst).lower(), dispfittype=effective_dispfittype, minfeaturereads=effective_minfeaturereads)

            build_flags = toolsTG.sanitize_flags(args_variant)
            merge_variant_into_adata(adata, contribution, tag=tag, direction=direction, cutoff=cutoff, build_flags=build_flags, overwrite=args.overwrite)
    finally:
        if not savesplitbams:
            for tag in (f'u{cutoff}', f'o{cutoff}'):
                if split_dirs_preexisted.get(tag, False):
                    # This tag's directory already had every sample's split BAM before this run
                    # started, and this run didn't force a fresh split -- it's not this run's
                    # scratch to clean up, regardless of --savesplitbams here.
                    continue
                split_dir = os.path.join(effective_bamdir, tag)
                if os.path.isdir(split_dir):
                    shutil.rmtree(split_dir, ignore_errors=True)

    output_path = os.path.abspath(args.output) if args.output else os.path.abspath(args.anndata)
    toolsTG.write_h5ad(adata, output_path)
    logger.info(f'Writing h5ad database object to {output_path}')
    return output_path
