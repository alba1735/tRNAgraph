#!/usr/bin/env python3

import itertools
import sys
import os
import subprocess
import time
import re
import random
import string
import tempfile
import logging
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple
import pysam

from . import toolsTG
from . import toolsDedup

logger = logging.getLogger(__name__)


def isprimarymapping(mapping):
    return not (mapping.flag & 0x0100 > 0)

def process_mappings(trnafile_path, infile, outfile, logfile=sys.stderr, progname=None, fqname=None, libname=None, minnontrnasize=20):
    trnadata = toolsTG.transcriptfile(trnafile_path)
    trnatranscripts = set(trnadata.gettranscripts())
    
    totalreads = 0
    multimaps = 0
    duperemove = 0
    shortened = 0
    mapsremoved = 0
    totalmaps = 0
    
    trnareads = 0
    maxreads = 0
    diffreads = 0
    
    uniquenontrnas = 0
    nonuniquenontrnas = 0
    
    ambanticodon = 0
    ambamino = 0
    ambtrna = 0
    
    imperfect = 0
    extraimperfect = 0
    
    maxmaps = 50 # Default from choosemappings.py

    try:
        bamfile = pysam.AlignmentFile(infile, "r")
    except Exception as e:
        print(f"Error opening input stream: {e}", file=logfile)
        return None

    newheader = bamfile.header.to_dict()
    newheader["RG"] = list()
    newheader["RG"].append(dict())
    
    if "PG" not in newheader:
        newheader["PG"] = list()
    if progname is not None:
        newheader["PG"].append({"PN" :progname, "ID": progname})
    if fqname is not None:
        newheader["RG"][0]["ID"] = fqname 
    if libname is not None:
        newheader["RG"][0]["LB"] = libname
        
    try:
        outbam = pysam.AlignmentFile(outfile, "wb", header=newheader)
    except Exception as e:
        print(f"Error opening output stream: {e}", file=logfile)
        return None

    for pairedname, allmaps in itertools.groupby(bamfile, lambda x: x.qname):
        allmaps = list(allmaps)
        if sum(curr.flag & 0x004 > 0 for curr in allmaps):
            continue
        totalreads += 1
        
        readlength = None
        clipsize = 50
        mappings = 0
        currscore = None
        newset = set()
        
        for currmap in allmaps:
            tagdict = dict()
            for curr in currmap.tags:
                tagdict[curr[0]] = curr[1]
            totalmaps += 1
            if currmap.tid == -1:
                continue
            
            readlength = len(currmap.seq)
            mappings += 1
            
            if currscore is None or currscore < tagdict["AS"]:
                newset = set()
                newset.add(currmap)
                currscore = tagdict["AS"]
            elif currscore == tagdict["AS"]:
                newset.add(currmap)
            else:
                pass
        
        if mappings > 1:
            multimaps += 1
        if len(newset) < mappings:
            shortened += 1
            
        if len(newset) >= 50:
            maxreads += 1
            
        finalset = list()
        
        if sum(bamfile.getrname(curr.tid) in trnatranscripts for curr in newset) > 0:
            trnareads += 1
            diff = len(newset) - sum(bamfile.getrname(curr.tid) in trnatranscripts for curr in newset)
            anticodons = frozenset(trnadata.getanticodon(bamfile.getrname(curr.tid)) for curr in newset if bamfile.getrname(curr.tid) in trnatranscripts)
            aminos = frozenset(trnadata.getamino(bamfile.getrname(curr.tid)) for curr in newset if bamfile.getrname(curr.tid) in trnatranscripts)
            trnamappings = list(curr for curr in newset if bamfile.getrname(curr.tid) in trnatranscripts)
            locusmaps = list(itertools.chain.from_iterable(trnadata.transcriptdict[bamfile.getrname(curr.tid)] for curr in trnamappings))
            
            if trnamappings[0].get_tag("XM") + trnamappings[0].get_tag("XO") > 0:
                imperfect += 1
            if trnamappings[0].get_tag("XM") + trnamappings[0].get_tag("XO") > 2:
                extraimperfect += 1
            
            if len(anticodons - frozenset(['NNN'])) > 1:
                ambanticodon += 1
            
            if len(aminos - frozenset(['Und'])) > 1:
                ambamino += 1
                
            if diff > 0:
                diffreads += 1
            
            if len(trnamappings) > 1:
                ambtrna += 1
            
            for currtrnamap in trnamappings:
                currtrnamap.tags = currtrnamap.tags + [("YA",len(anticodons))] + [("YM",len(aminos))]  + [("YR",len(trnamappings))] +  [("YL",len(locusmaps))]
            finalset = trnamappings
            
        else:
            #for non-tRNA, remove reads that are too small
            if readlength < minnontrnasize:
                continue
            
            for curr in newset:
                finalset.append(curr)
                
            if len(finalset) > maxmaps:
                duperemove += 1
                continue
            if len(newset) > 1:
                nonuniquenontrnas += 1
            else:
                uniquenontrnas += 1
        
        mapsremoved += mappings - len(finalset)
        if sum(isprimarymapping(curr) for curr in finalset) < 1:
            for i, curr in enumerate(finalset):
                if i == 0:
                    curr.flag &= ~ 0x0100
                    outbam.write(curr)
                else:
                    outbam.write(curr)
        else:
            for curr in finalset:
                outbam.write(curr)
    
    outbam.close()
    bamfile.close()
    
    return TRNAMapInfo(ambtrna, ambanticodon, ambamino, trnareads, uniquenontrnas, nonuniquenontrnas)

class TRNAMapInfo:
    def __init__(self, multtrans, multac, multamino, trna, singlenon, multiplenon):
        self.multtrans = int(multtrans)
        self.multac = int(multac)
        self.multamino = int(multamino)
        self.trna = int(trna)
        self.singlenon = int(singlenon)
        self.multiplenon = int(multiplenon)
        self.multitrna = (self.multtrans + self.multac + self.multamino)
        self.singletrna = self.trna - self.multitrna
        
    def uniquereads(self):
       return self.singletrna + self.singlenon
    def nonuniquereads(self):
       return self.multitrna + self.multiplenon

class MapInfo:
    def __init__(self, singlemap, multimap, unmap, totalreads, bowtietext, samplename, failedrun=False, bowtiecommand=None, trnamapinfo=None):
        self.unmaps = unmap
        self.bowtiesinglemaps = singlemap
        self.bowtiemultimaps = multimap
        self.totalreads = totalreads
        self.bowtietext = bowtietext
        self.samplename = samplename
        self.failedrun = failedrun
        self.bowtiecommand = bowtiecommand
        
        self.trnamapinfo = trnamapinfo
        if trnamapinfo is not None:
            self.singlemaps = trnamapinfo.uniquereads()
            self.multimaps = trnamapinfo.nonuniquereads()
        else: #has to use the raw bowtie2 output if no tRNA data
            self.singlemaps = singlemap
            self.multimaps = multimap
        self.unmap = int(self.totalreads) - (int(self.multimaps) + int(self.singlemaps))
        
    def printbowtie(self, logfile=sys.stderr):
        print("******************************************************************", file=logfile)
        print(self.bowtiecommand, file=logfile)
        print(self.bowtietext, file=logfile)

class MapReads:
    def __init__(self, bowtiedb, trnafile, minnontrnasize=20, local=False, maxmaps=100, program='bowtie2', threads=None):
        self.logger = logging.getLogger(__name__)
        self.bowtiedb = bowtiedb
        self.trnafile = trnafile
        self.minnontrnasize = minnontrnasize
        self.local = local
        self.maxmaps = maxmaps
        self.program = program
        if threads:
            self.threads = threads
        else:
            try:
                self.threads = min(8, cpu_count())
            except Exception:
                self.threads = 1

    def map_sample(self, samplename, fastqfile, bamfile, expname, logfile=None):
        localmode = []
        if self.local:
            localmode = ["--local"]
        
        # Construct bowtie2 command list
        bowtie_args = [self.program] + localmode + ["-x", self.bowtiedb, "-k", str(self.maxmaps), "--very-sensitive", "--ignore-quals", "--np", "5", "--reorder", "-p", str(self.threads), "-U", fastqfile]
        
        temploc = os.path.basename(bamfile) + ''.join(random.choice(string.ascii_lowercase) for i in range(8))
        self.logger.info(temploc)
        
        # Construct samtools command list
        samtools_args = ["samtools", "sort", "-T", f"{toolsTG.sort_temp_dir(bamfile)}/{temploc}temp", "-", "-o", f"{bamfile}.bam"]
        
        bowtie_cmd_str = " ".join(bowtie_args)
        self.logger.info(bowtie_cmd_str)
        if logfile:
            print(bowtie_cmd_str, file=logfile)
            logfile.flush()
            
        # Use a temporary file for bowtie2 stderr to avoid deadlock
        with tempfile.TemporaryFile(mode='w+') as bowtie_err_file:
            try:
                # Start bowtie2 (stdout is text SAM)
                bowtie_proc = subprocess.Popen(bowtie_args, stdout=subprocess.PIPE, stderr=bowtie_err_file, universal_newlines=True)
                
                # Start samtools (stdin is binary BAM)
                samtools_proc = subprocess.Popen(samtools_args, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
                
                # Process mappings
                trnamapinfo = process_mappings(self.trnafile, bowtie_proc.stdout, samtools_proc.stdin, 
                                               logfile=logfile if logfile else sys.stderr, 
                                               progname="tRNAgraph", fqname=fastqfile, libname=expname, 
                                               minnontrnasize=self.minnontrnasize)
                
                # Wait for processes
                bowtie_proc.wait()
                # samtools_proc.stdin is closed by process_mappings via pysam
                samtools_out, samtools_err = samtools_proc.communicate()
                
                # Read bowtie2 stderr
                bowtie_err_file.seek(0)
                errinfo = bowtie_err_file.read()
                
                if logfile is not None:
                    print(errinfo, file=logfile)
                    if samtools_err:
                        print(samtools_err.decode('utf-8', errors='replace'), file=logfile)
                    logfile.flush()
            
                if bowtie_proc.returncode != 0:
                    return MapInfo(0, 0, 0, 0, errinfo, samplename, failedrun=True, bowtiecommand=bowtie_cmd_str)

                # Parse bowtie2 output
                rereadtotal = re.search(r'(\d+).*reads', errinfo)
                rereadunmap = re.search(r'\s*(\d+).*0 times', errinfo)
                rereadsingle = re.search(r'\s*(\d+).*exactly 1 time', errinfo)
                rereadmult = re.search(r'\s*(\d+).*>1 times', errinfo)
                
                if rereadtotal and rereadunmap and rereadsingle and rereadmult:
                    totalreads = rereadtotal.group(1)
                    unmappedreads = rereadunmap.group(1)
                    singlemaps = rereadsingle.group(1)
                    multmaps = rereadmult.group(1)
                    return MapInfo(singlemaps, multmaps, unmappedreads, totalreads, errinfo, samplename, bowtiecommand=bowtie_cmd_str, trnamapinfo=trnamapinfo)
                else:
                    self.logger.error(f"Could not map {fastqfile}, check mapstats file")
                    self.logger.error("Exiting...")
                    self.logger.error(errinfo)
                    return MapInfo(0, 0, 0, 0, errinfo, samplename, failedrun=True, bowtiecommand=bowtie_cmd_str)

            except Exception as e:
                self.logger.error(f"Error during mapping: {e}")
                if logfile:
                    print(f"Error during mapping: {e}", file=logfile)
                return MapInfo(0, 0, 0, 0, str(e), samplename, failedrun=True, bowtiecommand=bowtie_cmd_str)

    @staticmethod
    def checkheaders(bamname, fqname):
        try:
            bamfile = pysam.AlignmentFile(bamname, "rb")
        except ValueError:
            return True
        except IOError as e:
            logger.error(f"Failed to read {bamname}")
            logger.error(e)
            sys.exit(1)
        newheader = bamfile.header
        if 'PG' in newheader and len(newheader["PG"]) > 1 and newheader["PG"][1]["PN"] == "tRNAgraph":
            if newheader["RG"][0]["ID"] != fqname:
                return False
        return True

# Concurrency ceiling for the deduplication phase. Lower than the mapping pool's because
# umi_tools holds a position's reads in memory while resolving its UMI network; the manual
# tRAX-era workflow serialised samples entirely for that reason.
DEDUP_MAX_CONCURRENCY = 4


def dedup_sample_wrapper(args):
    """Runs one sample's deduplication in a pool worker.

    Returns None on success or an error string, rather than raising: an exception crossing a
    multiprocessing boundary loses its traceback, and a partial failure needs to name which
    sample failed. dedup_sample() restores the original bam before it raises, so a failed
    sample leaves its mapping intact.
    """
    bam_path, method, keep_prededup = args
    try:
        toolsDedup.dedup_sample(bam_path, method=method, keep_prededup=keep_prededup)
    except Exception as exc:
        return f"Deduplication failed for {os.path.basename(bam_path)}: {exc}"
    return None


def map_sample_wrapper(args):
    mapper, samplename, fastqfile, bamfile, expname = args
    return mapper.map_sample(samplename, fastqfile, bamfile, expname)


class trnadatabase:
    def __init__(self, dbname):
        self.dbname = dbname
        self.trnatable = dbname+"-trnatable.txt"
        self.bowtiedb = dbname+"-tRNAgenome"
        self.locifile = dbname+"-trnaloci.bed"
        self.maturetrnas=dbname+"-maturetRNAs.bed"
        self.trnaalign = dbname+"-trnaalign.stk"
        self.locialign = dbname+"-trnaloci.stk"
        self.trnanums = dbname+"-alignnum.txt"
        self.locinums = dbname+"-locusnum.txt" 
        self.trnafasta = dbname+"-maturetRNAs.fa"
        self.modomics = dbname+"-modomics.txt"
        self.otherseqs = dbname+"-otherseqs.txt"
        self.dbinfo = dbname+"-dbinfo.txt"
    
    def getorgtype(self):
        """Organism mode recorded by `preprocess makedb`, defaulting to euk.

        Blank lines are skipped -- `fields[0]` on an empty split() raised
        IndexError. A recorded mode that is not one tRNAgraph knows about is an
        error rather than a silent fall back to eukaryotic positions, which would
        otherwise produce plausible-looking but wrong Sprinzl numbering.
        """
        from . import toolsGetCoverage

        orgtype = "euk"
        if os.path.exists(self.dbinfo):
            for currline in open(self.dbinfo):
                fields = currline.split()
                if len(fields) >= 2 and fields[0] == "orgmode":
                    orgtype = fields[1]
        if orgtype not in toolsGetCoverage.POSITION_TABLES:
            raise ValueError(
                f"Database {self.dbinfo} records an unknown organism mode "
                f"{orgtype!r}. Expected one of: "
                + ", ".join(sorted(toolsGetCoverage.POSITION_TABLES))
            )
        return orgtype

class expdatabase:
    def __init__(self, expname, results_dir_name="results", graphs_dir_name="graphs"):
        self.expname = expname.rstrip(os.sep)
        basename = os.path.basename(self.expname)
        self.resultsdir = os.path.join(self.expname, results_dir_name)
        self.graphsdir = os.path.join(self.expname, graphs_dir_name)
        
        # Helper to join paths
        def res_path(path): return os.path.join(self.resultsdir, path)
        def graph_path(path): return os.path.join(self.graphsdir, path)

        self.uniquename = res_path("unique/"+basename+"-unique")
        self.allfeats = res_path(basename+"-allfeats.bed")
        
        self.mapinfo = res_path(basename+"-mapinfo.txt")
        self.mapplot = graph_path(basename+"-mapinfo.pdf")

        self.trnamapfile = res_path(basename+"-trnamapinfo.txt")
        self.trnamapplot = graph_path(basename+"-trnamapinfo.pdf")
        
        self.maplog = res_path(basename+"-mapstats.txt")
        self.dedupinfo = res_path(basename+"-dedupinfo.txt")
        self.dedupstats = res_path(basename+"-dedupstats.txt")
        self.replicatecorrelation = res_path(basename+"-replicatecorrelation.txt")
        self.genetypes = res_path(basename+"-genetypes.txt")
        self.genecounts = res_path(basename+"-readcounts.txt")
        self.trnacounts = res_path(basename+"-trnacounts.txt")
        
        self.normalizedcounts = res_path(basename+"-normalizedreadcounts.txt")
        self.normalizedtrnacounts = res_path("trna/"+basename+"-trna_normalizedreadcounts.txt")
        self.sizefactors = res_path(basename+"-SizeFactors.txt")
        self.trnasizefactors = res_path("trna/"+basename+"-trna_SizeFactors.txt")
        self.allfeaturesizefactors = res_path("allfeature/"+basename+"-allfeature_SizeFactors.txt")
        self.normalizedcounts_allfeatures = res_path("allfeature/"+basename+"-allfeature_normalizedreadcounts.txt")

        self.genetypecounts=res_path(basename+"-typecounts.txt")
        self.genetypeplot=graph_path(basename+"-typecounts.pdf")

        self.genetyperealcounts=res_path(basename+"-typerealcounts.txt")
        self.genetyperealplot=graph_path(basename+"-typerealcounts.pdf")
        
        self.trnaaminofile=res_path(basename+"-aminocounts.txt")
        self.trnauniqaminofile=res_path("unique/"+basename+"-unique-aminos.txt") 

        self.trnaaminoplot=graph_path(basename+"-aminocounts.pdf")
        self.trnaaminorealplot=graph_path(basename+"-aminorealcounts.pdf")
        
        self.trnaanticodonfile=res_path(basename+"-anticodoncounts.txt")
        self.trnauniqanticodonfile=res_path("unique/"+basename+"-unique-anticodons.txt")
        self.trnauniqcountsfile=res_path("unique/"+basename+"-unique-trnas.txt")
        
        self.trnalengthfile=res_path(basename+"-readlengths.txt")
        self.trnalengthplot=graph_path(basename+"-readlengths.pdf")
        
        self.mismatchcountfile=res_path(basename+"-mismatches.txt")

        self.trnacoveragefile=res_path(basename+"-coverage.txt")
        self.trnacoverageplot=graph_path(basename+"-coverage.pdf")
        self.trnacombinecoverageplot=graph_path(basename+"-combinecoverage.pdf")
        
        self.trnauniqcoveragefile=res_path(basename+"-uniqcoverage.txt")

        self.locicoveragefile=res_path("pretRNAs/"+basename+"-pretRNAcoverage.txt")
        self.locicoverageplot=graph_path("pretRNAs/"+basename+"-pretRNAcoverage.pdf")
        self.locicombinecoverageplot=graph_path("pretRNAs/"+basename+"-pretRNAcombinecoverage.pdf")
        
        # Two separate artifacts, deliberately: the table is tRAX's own
        # `newcoverageplots.R:581` output, the BED is `getgenomicmismatches.py`'s. See
        # toolsGetCoverage's SIGMISMATCH_* constants for why they don't share thresholds.
        self.sigmismatchfile = res_path("mismatch/"+basename+"-sigmismatch.txt")
        self.sigmismatchbed = res_path("mismatch/"+basename+"-sigmismatch.bed")

        self.trnauniquefile=res_path("unique/"+basename+"-trnauniquecounts.txt")
        self.trnaendfile=res_path(basename+"-trnaendcounts.txt")
        
        self.pcaplot = graph_path(basename+"-pca.pdf")
        self.pcatrnaplot = graph_path(basename+"-pcatrna.pdf")
        self.pcaacplot = graph_path("unique/"+basename+"-pcaac.pdf")

        self.qaoutputname = res_path(basename+"-qa.html")

class MapSamples:
    def __init__(self, args):
        self.logger = logging.getLogger(__name__)
        self.args = args
        self.dbname = args.database
        self.expname = args.output
        self.samplefilename = args.input
        self.force_remap = getattr(args, 'force_remap', False)
        self.bamdir = args.bamdir if args.bamdir else os.path.join("processed", "bam")
        if args.threads:
            self.cores = args.threads
        else:
            try:
                self.cores = min(8, cpu_count())
            except Exception:
                self.cores = 8
        self.minnontrnasize = args.minnontrnasize
        self.local = args.local
        self.skipfqcheck = args.skipcheck
        self.quiet = getattr(args, 'quiet', False)
        # Opt-in: deduplication changes counts and cannot be reversed without remapping, and on
        # a dataset without UMIs it would remove genuine reads. Read defensively because several
        # callers (the test suite, addsplit) build this namespace by hand.
        self.dedup = getattr(args, 'dedup', False)
        self.keep_prededup = getattr(args, 'keep_prededup', False)
        self.dedup_method = getattr(args, 'dedup_method', toolsDedup.DEFAULT_DEDUP_METHOD)

        self.trnainfo = trnadatabase(self.dbname)
        self.expinfo = expdatabase(self.expname)
        
    def main(self):
        # Create directories
        if not os.path.exists(self.expname):
            os.makedirs(self.expname)
        if not os.path.exists(self.bamdir):
            os.makedirs(self.bamdir)
            
        # Create results directory (graphsdir is never written to by this command -- nothing
        # to create it for)
        if not os.path.exists(self.expinfo.resultsdir):
            os.makedirs(self.expinfo.resultsdir)

        # Expand dbname
        self.dbname = os.path.expanduser(self.dbname)
        
        # Mapping Reads
        self.logger.info("Mapping Reads")
        self.mapsamples()

        # Deduplication is a separate phase rather than part of mapsamples() because it operates
        # on the finished, sorted bams -- including any that mapsamples() skipped as already
        # present -- and because umi_tools is memory-hungry enough that it should not share the
        # mapping pool.
        if self.dedup:
            self.dedupsamples()

    def dedupsamples(self):
        sampledata = toolsTG.samplefile(self.samplefilename)
        samples = sampledata.getsamples()

        # umi_tools is single-threaded, so the only parallelism available is running samples
        # concurrently. The cap is deliberately below the mapping pool's: umi_tools holds a
        # position's reads in memory while resolving its UMI network, and the manual workflow this
        # replaces (tRAX-era dedup.bash) serialised the samples explicitly for that reason. Four
        # concurrent jobs takes a measured 2.4h nine-sample human run to roughly the length of its
        # slowest single sample (~32 min) without approaching bowtie2's footprint.
        pool_size = min(DEDUP_MAX_CONCURRENCY, max(1, self.cores // 8), len(samples))

        pending = []
        for samplename in samples:
            bamfile = os.path.join(self.bamdir, samplename + ".bam")
            if not os.path.isfile(bamfile):
                self.logger.warning(f"No bam file found for {samplename}, skipping deduplication")
                continue
            pending.append((samplename, bamfile))

        if not pending:
            self.write_dedupinfo([])
            return

        self.logger.info(
            f"Deduplicating {len(pending)} samples by UMI "
            f"(method={self.dedup_method}, {pool_size} concurrent jobs)"
        )

        # Separators are detected up front, in the parent: it is a cheap read of the first reads of
        # each bam, and doing it here means a missing UMI is refused before any sample has been
        # modified, rather than after some have already been deduplicated.
        records = []
        for samplename, bamfile in pending:
            records.append((samplename, bamfile, toolsDedup.detect_umi_separator(bamfile)))

        args = [(bamfile, self.dedup_method, self.keep_prededup) for _, bamfile, _ in records]
        if pool_size == 1:
            results = [dedup_sample_wrapper(a) for a in toolsTG.progress_iterator(
                args, total=len(args), desc="Deduplicating samples",
                logger=self.logger, quiet=self.quiet)]
        else:
            with Pool(processes=pool_size) as pool:
                results = list(toolsTG.progress_iterator(
                    pool.imap_unordered(dedup_sample_wrapper, args),
                    total=len(args), desc="Deduplicating samples",
                    logger=self.logger, quiet=self.quiet,
                ))

        failures = [r for r in results if r is not None]
        if failures:
            for message in failures:
                self.logger.error(message)
            sys.exit(1)

        self.write_dedupinfo(records)

    def write_dedupinfo(self, records):
        '''
        Records that deduplication happened, and how.

        `map` writes no runinfo of its own -- that is an `analyze build` artifact -- so without
        this a deduplicated bam directory is indistinguishable from a non-deduplicated one, even
        though every downstream count differs. Per-sample read totals are deliberately not
        recomputed here: they are in umi_tools' own per-sample log next to each bam, and
        recovering them at this level would mean a second full pass over every file.
        '''
        version = toolsDedup.umi_tools_version()
        os.makedirs(self.expinfo.resultsdir, exist_ok=True)
        with open(self.expinfo.dedupinfo, 'w') as info:
            print(f"umi_tools_version\t{version}", file=info)
            print(f"method\t{self.dedup_method}", file=info)
            print(f"keep_prededup\t{self.keep_prededup}", file=info)
            print("sample\tbam\tumi_separator\tumi_tools_log", file=info)
            for samplename, bamfile, separator in records:
                logpath = toolsDedup.dedup_log_path(bamfile)
                print(f"{samplename}\t{bamfile}\t{separator}\t{logpath}", file=info)

        self.write_dedupstats(records)

    def write_dedupstats(self, records):
        '''
        Per-sample deduplication statistics, as a tab-separated table.

        Separate from dedupinfo because they answer a different question: dedupinfo records what
        was run, this records whether it worked. Read across samples it is the fastest way to
        spot one that should be dropped -- an outlying `reads_per_molecule` means uneven PCR
        amplification, a low `reads_per_position` with few `positions` means a low-complexity
        library, and `max_umis_per_position` approaching 4^n for an n-base UMI means the tag
        space is saturating and counts at the deepest features are compressed.

        Every number here comes from umi_tools' own end-of-run log, so this costs a small text
        parse rather than another pass over the bams.
        '''
        rows = [toolsDedup.dedup_stats_row(name, toolsDedup.dedup_log_path(bam))
                for name, bam, _ in records]

        def fmt(value):
            if value is None:
                return 'NA'
            return f"{value:.2f}" if isinstance(value, float) else str(value)

        with open(self.expinfo.dedupstats, 'w') as out:
            print('\t'.join(toolsDedup.DEDUP_STATS_COLUMNS), file=out)
            for row in rows:
                print('\t'.join(fmt(row[c]) for c in toolsDedup.DEDUP_STATS_COLUMNS), file=out)

        for row in rows:
            if row['reads_per_molecule'] is not None:
                self.logger.info(
                    f"  {row['sample']}: {row['retained_pct']:.1f}% retained, "
                    f"{row['reads_per_molecule']:.2f} reads/molecule, "
                    f"{row['reads_per_position']:.2f} reads/position"
                )

    def mapsamples(self):
        # Calculate resource allocation
        total_cores = self.cores
        if total_cores > 1:
             pool_size = max(1, total_cores // 8)
             bowtie_threads = total_cores // pool_size
        else:
             pool_size = 1
             bowtie_threads = 1
        self.logger.info(f"Mapping with {pool_size} concurrent jobs, {bowtie_threads} threads per job.")

        # Initialize MapReads
        mapper = MapReads(bowtiedb=self.trnainfo.bowtiedb, trnafile=self.trnainfo.trnatable, 
                          minnontrnasize=self.minnontrnasize, 
                          local=self.local, threads=bowtie_threads)
        
        # Get samples
        sampledata = toolsTG.samplefile(self.samplefilename)
        samples = sampledata.getsamples()
        
        # Prepare arguments for multiprocessing
        map_args = []
        for samplename in samples:
            fastqfile = sampledata.getfastq(samplename)
            bamfile = os.path.join(self.bamdir, samplename)
            
            if not self.force_remap and os.path.isfile(bamfile + ".bam"):
                if not MapReads.checkheaders(bamfile + ".bam", fastqfile):
                    self.logger.warning(f"Bam file {bamfile}.bam does not match fq file {fastqfile}")
                    if not self.skipfqcheck:
                        sys.exit(1)
                self.logger.warning(f"Skipping {samplename}")
            else:
                map_args.append((mapper, samplename, fastqfile, bamfile, self.expname))

        # Run mapping
        # Open mapstats log file
        maplog_file = None
        if self.expinfo.maplog:
            try:
                maplog_file = open(self.expinfo.maplog, 'w')
            except IOError as e:
                self.logger.warning(f"Could not open mapstats file {self.expinfo.maplog}: {e}")

        # Run mapping
        mapresults = {}
        if map_args:
            with Pool(processes=pool_size) as pool:
                for result in toolsTG.progress_iterator(
                    pool.imap_unordered(map_sample_wrapper, map_args),
                    total=len(map_args), desc="Mapping samples", logger=self.logger,
                    quiet=self.quiet,
                ):
                    if result.failedrun:
                        self.logger.error("Failure to Bowtie2 map")
                        result.printbowtie()
                        if maplog_file:
                            result.printbowtie(logfile=maplog_file)
                        sys.exit(1)
                    mapresults[result.samplename] = result
                    result.printbowtie()
                    if maplog_file:
                        result.printbowtie(logfile=maplog_file)
        
        if maplog_file:
            maplog_file.close()
        
        # Write mapinfo
        if self.expinfo.mapinfo:
            try:
                with open(self.expinfo.mapinfo, 'w') as mapinfo:
                    print("\t".join(samples), file=mapinfo)
                    print("unmap\t" + "\t".join(str(mapresults[s].unmaps) if s in mapresults else "0" for s in samples), file=mapinfo)
                    print("single\t" + "\t".join(str(mapresults[s].singlemaps) if s in mapresults else "0" for s in samples), file=mapinfo)
                    print("multi\t" + "\t".join(str(mapresults[s].multimaps) if s in mapresults else "0" for s in samples), file=mapinfo)
                    # TODO: Uncomment if totalreads and bowtiecommand needed in mapinfo - Currently commented out to match trax format
                    # print("total\t" + "\t".join(str(mapresults[s].totalreads) if s in mapresults else "0" for s in samples), file=mapinfo)
                    # print("bowtiecommand\t" + "\t".join(str(mapresults[s].bowtiecommand) if s in mapresults else "" for s in samples), file=mapinfo)
            except IOError as e:
                self.logger.warning(f"Could not write mapinfo file {self.expinfo.mapinfo}: {e}")
        
        # Write trna mapinfo
        if self.expinfo.trnamapfile and self.expinfo.mapinfo:
            try:
                with open(self.expinfo.trnamapfile, 'w') as trnamapinfo:
                    print("\t".join(samples), file=trnamapinfo)
                    print("multtrans\t" + "\t".join(str(mapresults[s].trnamapinfo.multtrans) if s in mapresults and mapresults[s].trnamapinfo is not None else "0" for s in samples), file=trnamapinfo)
                    print("multac\t" + "\t".join(str(mapresults[s].trnamapinfo.multac) if s in mapresults and mapresults[s].trnamapinfo is not None else "0" for s in samples), file=trnamapinfo)
                    print("multamino\t" + "\t".join(str(mapresults[s].trnamapinfo.multamino) if s in mapresults and mapresults[s].trnamapinfo is not None else "0" for s in samples), file=trnamapinfo)
                    print("trna\t" + "\t".join(str(mapresults[s].trnamapinfo.trna) if s in mapresults and mapresults[s].trnamapinfo is not None else "0" for s in samples), file=trnamapinfo)
                    print("singlenon\t" + "\t".join(str(mapresults[s].trnamapinfo.singlenon) if s in mapresults and mapresults[s].trnamapinfo is not None else "0" for s in samples), file=trnamapinfo)
                    print("multiplenon\t" + "\t".join(str(mapresults[s].trnamapinfo.multiplenon) if s in mapresults and mapresults[s].trnamapinfo is not None else "0" for s in samples), file=trnamapinfo)
                    print("total\t" + "\t".join(str(mapresults[s].trnamapinfo.uniquereads() + mapresults[s].trnamapinfo.nonuniquereads()) if s in mapresults and mapresults[s].trnamapinfo is not None else "0" for s in samples), file=trnamapinfo)
            except IOError as e:
                self.logger.warning(f"Could not write trnamapinfo file {self.expinfo.trnamapfile}: {e}")
