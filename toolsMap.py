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
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple
import pysam

from . import toolsTG



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
        print(temploc, file=sys.stderr)
        
        # Construct samtools command list
        samtools_args = ["samtools", "sort", "-T", f"{tempfile.gettempdir()}/{temploc}temp", "-", "-o", f"{bamfile}.bam"]
        
        bowtie_cmd_str = " ".join(bowtie_args)
        print(bowtie_cmd_str, file=sys.stderr)
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
                    print(f"Could not map {fastqfile}, check mapstats file", file=sys.stderr)
                    print("Exiting...", file=sys.stderr)
                    print(errinfo, file=sys.stderr)
                    return MapInfo(0, 0, 0, 0, errinfo, samplename, failedrun=True, bowtiecommand=bowtie_cmd_str)

            except Exception as e:
                print(f"Error during mapping: {e}", file=sys.stderr)
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
            print(f"Failed to read {bamname}", file=sys.stderr)
            print(e, file=sys.stderr)
            sys.exit(1)
        newheader = bamfile.header
        if 'PG' in newheader and len(newheader["PG"]) > 1 and newheader["PG"][1]["PN"] == "tRNAgraph":
            if newheader["RG"][0]["ID"] != fqname:
                return False
        return True

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
        orgtype = "euk"
        if os.path.exists(self.dbinfo):
            for currline in open(self.dbinfo):
                fields = currline.split()
                if fields[0] == "orgmode":
                    orgtype = fields[1]
        return orgtype

class expdatabase:
    def __init__(self, expname):
        self.expname = expname
        self.resultsdir = os.path.join(expname, "results")
        self.graphsdir = os.path.join(expname, "graphs")
        
        # Helper to join paths
        def res_path(path): return os.path.join(self.resultsdir, path)
        def graph_path(path): return os.path.join(self.graphsdir, path)

        self.uniquename = res_path("unique/"+expname+"-unique")
        self.allfeats = res_path(expname+"-allfeats.bed")
        
        self.mapinfo = res_path(expname+"-mapinfo.txt")
        self.mapplot = graph_path(expname+"-mapinfo.pdf")

        self.trnamapfile = res_path(expname+"-trnamapinfo.txt")
        self.trnamapplot = graph_path(expname+"-trnamapinfo.pdf")
        
        self.maplog = res_path(expname+"-mapstats.txt")
        self.genetypes = res_path(expname+"-genetypes.txt")
        self.genecounts = res_path(expname+"-readcounts.txt")
        self.trnacounts = res_path(expname+"-trnacounts.txt")
        
        self.normalizedcounts = res_path(expname+"-normalizedreadcounts.txt")
        self.normalizedtrnacounts = res_path("trna/"+expname+"-trna_normalizedreadcounts.txt")
        self.sizefactors = res_path(expname+"-SizeFactors.txt")
        self.trnasizefactors = res_path("trna/"+expname+"-trna_SizeFactors.txt")

        self.genetypecounts=res_path(expname+"-typecounts.txt")
        self.genetypeplot=graph_path(expname+"-typecounts.pdf")

        self.genetyperealcounts=res_path(expname+"-typerealcounts.txt")
        self.genetyperealplot=graph_path(expname+"-typerealcounts.pdf")
        
        self.trnaaminofile=res_path(expname+"-aminocounts.txt")
        self.trnauniqaminofile=res_path("unique/"+expname+"-unique-aminos.txt") 

        self.trnaaminoplot=graph_path(expname+"-aminocounts.pdf")
        self.trnaaminorealplot=graph_path(expname+"-aminorealcounts.pdf")
        
        self.trnaanticodonfile=res_path(expname+"-anticodoncounts.txt")
        self.trnauniqanticodonfile=res_path("unique/"+expname+"-unique-anticodons.txt")
        self.trnauniqcountsfile=res_path("unique/"+expname+"-unique-trnas.txt")
        
        self.trnalengthfile=res_path(expname+"-readlengths.txt")
        self.trnalengthplot=graph_path(expname+"-readlengths.pdf")
        
        self.mismatchcountfile=res_path(expname+"-mismatches.txt")
        self.mismatchcountplot=graph_path(expname+"-mismatches.pdf")
        
        self.trnacoveragefile=res_path(expname+"-coverage.txt")
        self.trnacoverageplot=graph_path(expname+"-coverage.pdf")
        self.trnacombinecoverageplot=graph_path(expname+"-combinecoverage.pdf")
        
        self.trnauniqcoveragefile=res_path(expname+"-uniqcoverage.txt")

        self.locicoveragefile=res_path("pretRNAs/"+expname+"-pretRNAcoverage.txt")
        self.locicoverageplot=graph_path("pretRNAs/"+expname+"-pretRNAcoverage.pdf")
        self.locicombinecoverageplot=graph_path("pretRNAs/"+expname+"-pretRNAcombinecoverage.pdf")
        
        self.trnamismatchfile = res_path("mismatch/"+expname+"-mismatchcoverage.txt")
        self.sigmismatchfile = res_path("mismatch/"+expname+"-sigmismatch.txt")
        self.trnamismatchplot = graph_path("mismatch/"+expname+"-mismatchcoverage.pdf")
        
        self.trnadeletefile = res_path("mismatch/"+expname+"-deletecoverage.txt")
        self.trnadeleteplot = graph_path("mismatch/"+expname+"-deletecoverage.pdf")
        
        self.trnamismatchreport = res_path("mismatch/"+expname+"-mismatchreport.txt")
        self.trnauniquefile=res_path("unique/"+expname+"-trnauniquecounts.txt")
        self.trnaendfile=res_path(expname+"-trnaendcounts.txt")
        
        self.pcaplot = graph_path(expname+"-pca.pdf")
        self.pcatrnaplot = graph_path(expname+"-pcatrna.pdf")
        self.pcaacplot = graph_path("unique/"+expname+"-pcaac.pdf")

        self.qaoutputname = res_path(expname+"-qa.html")

class MapSamples:
    def __init__(self, args):
        self.args = args
        self.dbname = args.database
        self.expname = args.experiment
        self.samplefilename = args.samples
        self.lazyremap = args.lazy
        self.bamdir = args.bamdir if args.bamdir else os.path.join("bam", self.expname)
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
        
        self.trnainfo = trnadatabase(self.dbname)
        self.expinfo = expdatabase(self.expname)
        
    def main(self):
        # Create directories
        if not os.path.exists(self.expname):
            os.makedirs(self.expname)
        if not os.path.exists(self.bamdir):
            os.makedirs(self.bamdir)
            
        # Create results and graphs directories
        if not os.path.exists(self.expinfo.resultsdir):
            os.makedirs(self.expinfo.resultsdir)
        if not os.path.exists(self.expinfo.graphsdir):
            os.makedirs(self.expinfo.graphsdir)

        # Expand dbname
        self.dbname = os.path.expanduser(self.dbname)
        
        # Mapping Reads
        print("Mapping Reads", file=sys.stderr)
        self.mapsamples()

    def mapsamples(self):
        # Calculate resource allocation
        total_cores = self.cores
        if total_cores > 1:
             pool_size = max(1, total_cores // 8)
             bowtie_threads = total_cores // pool_size
        else:
             pool_size = 1
             bowtie_threads = 1
        print(f"Mapping with {pool_size} concurrent jobs, {bowtie_threads} threads per job.", file=sys.stderr)

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
            
            if self.lazyremap and os.path.isfile(bamfile + ".bam"):
                if not MapReads.checkheaders(bamfile + ".bam", fastqfile):
                    print(f"Bam file {bamfile}.bam does not match fq file {fastqfile}", file=sys.stderr)
                    if not self.skipfqcheck:
                        sys.exit(1)
                print(f"Skipping {samplename}", file=sys.stderr)
            else:
                map_args.append((mapper, samplename, fastqfile, bamfile, self.expname))

        # Run mapping
        # Open mapstats log file
        maplog_file = None
        if self.expinfo.maplog:
            try:
                maplog_file = open(self.expinfo.maplog, 'w')
            except IOError as e:
                print(f"Could not open mapstats file {self.expinfo.maplog}: {e}", file=sys.stderr)

        # Run mapping
        mapresults = {}
        if map_args:
            with Pool(processes=pool_size) as pool:
                for result in pool.imap_unordered(map_sample_wrapper, map_args):
                    if result.failedrun:
                        print("Failure to Bowtie2 map", file=sys.stderr)
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
                    print("total\t" + "\t".join(str(mapresults[s].totalreads) if s in mapresults else "0" for s in samples), file=mapinfo)
                    print("bowtiecommand\t" + "\t".join(str(mapresults[s].bowtiecommand) if s in mapresults else "" for s in samples), file=mapinfo)
            except IOError as e:
                print(f"Could not write mapinfo file {self.expinfo.mapinfo}: {e}", file=sys.stderr)
        
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
                print(f"Could not write trnamapinfo file {self.expinfo.trnamapfile}: {e}", file=sys.stderr)
