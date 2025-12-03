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
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pysam

import toolsTG
import toolsCountReads
import toolsGetCoverage
import toolsTrackHub



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
        self.uniquename = expname+"/unique/"+expname+"-unique"
        self.allfeats = expname+"/"+expname+"-allfeats.bed"
        
        self.mapinfo = expname+"/"+expname+"-mapinfo.txt"
        self.mapplot = expname+"/"+expname+"-mapinfo.pdf"

        self.trnamapfile = expname+"/"+expname+"-trnamapinfo.txt"
        self.trnamapplot = expname+"/"+expname+"-trnamapinfo.pdf"
        
        self.maplog = expname+"/"+expname+"-mapstats.txt"
        self.genetypes = expname+"/"+expname+"-genetypes.txt"
        self.genecounts = expname+"/"+expname+"-readcounts.txt"
        self.trnacounts = expname+"/"+expname+"-trnacounts.txt"
        
        self.normalizedcounts = expname+"/"+expname+"-normalizedreadcounts.txt"
        self.normalizedtrnacounts = expname+"/trna/"+expname+"-trna_normalizedreadcounts.txt"
        self.sizefactors = expname+"/"+expname+"-SizeFactors.txt"
        self.trnasizefactors = expname+"/trna/"+expname+"-SizeFactors.txt"

        self.genetypecounts=expname+"/"+expname+"-typecounts.txt"
        self.genetypeplot=expname+"/"+expname+"-typecounts.pdf"

        self.genetyperealcounts=expname+"/"+expname+"-typerealcounts.txt"
        self.genetyperealplot=expname+"/"+expname+"-typerealcounts.pdf"
        
        self.trnaaminofile=expname+"/"+expname+"-aminocounts.txt"
        self.trnauniqaminofile=expname+"/unique/"+expname+"-unique-aminos.txt" 

        self.trnaaminoplot=expname+"/"+expname+"-aminocounts.pdf"
        self.trnaaminorealplot=expname+"/"+expname+"-aminorealcounts.pdf"
        
        self.trnaanticodonfile=expname+"/"+expname+"-anticodoncounts.txt"
        self.trnauniqanticodonfile=expname+"/unique/"+expname+"-unique-anticodons.txt"
        self.trnauniqcountsfile=expname+"/unique/"+expname+"-unique-trnas.txt"
        
        self.trnalengthfile=expname+"/"+expname+"-readlengths.txt"
        self.trnalengthplot=expname+"/"+expname+"-readlengths.pdf"
        
        self.mismatchcountfile=expname+"/"+expname+"-mismatches.txt"
        self.mismatchcountplot=expname+"/"+expname+"-mismatches.pdf"
        
        self.trnacoveragefile=expname+"/"+expname+"-coverage.txt"
        self.trnacoverageplot=expname+"/"+expname+"-coverage.pdf"
        self.trnacombinecoverageplot=expname+"/"+expname+"-combinecoverage.pdf"
        
        self.trnauniqcoveragefile=expname+"/"+expname+"-uniqcoverage.txt"

        self.locicoveragefile=expname+"/pretRNAs/"+expname+"-pretRNAcoverage.txt"
        self.locicoverageplot=expname+"/pretRNAs/"+expname+"-pretRNAcoverage.pdf"
        self.locicombinecoverageplot=expname+"/pretRNAs/"+expname+"-pretRNAcombinecoverage.pdf"
        
        self.trnamismatchfile = expname+"/mismatch/"+expname+"-mismatchcoverage.txt"
        self.trnamismatchplot = expname+"/mismatch/"+expname+"-mismatchcoverage.pdf"
        
        self.trnadeletefile = expname+"/mismatch/"+expname+"-deletecoverage.txt"
        self.trnadeleteplot = expname+"/mismatch/"+expname+"-deletecoverage.pdf"
        
        self.trnamismatchreport = expname+"/mismatch/"+expname+"-mismatchreport.txt"
        self.trnauniquefile=expname+"/unique/"+expname+"-trnauniquecounts.txt"
        self.trnaendfile=expname+"/"+expname+"-trnaendcounts.txt"
        
        self.pcaplot = expname+"/"+expname+"-pca.pdf"
        self.pcatrnaplot = expname+"/"+expname+"-pcatrna.pdf"
        self.pcaacplot = expname+"/unique/"+expname+"-pcaac.pdf"

        self.qaoutputname = expname+"/"+expname+"-qa.html"

class MapSamples:
    def __init__(self, args):
        self.args = args
        self.dbname = args.database
        self.expname = args.experiment
        self.samplefilename = args.samples
        self.ensgtf = args.gtf
        self.bedfiles = args.bed
        self.lazyremap = args.lazy
        self.nofrag = args.nofrag
        self.nosizefactors = args.nosizefactors
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
        self.maxmismatches = args.maxmismatches
        self.mincoverage = args.mincoverage
        self.uniqueonlycov = args.uniqueonly
        self.pairfile = args.pairs
        self.paironly = args.paironly
        self.hubonly = args.hubonly
        self.makehubs = args.hub
        self.maponly = args.maponly
        self.dumpother = args.dumpother
        
        self.trnainfo = trnadatabase(self.dbname)
        self.expinfo = expdatabase(self.expname)
        
    def main(self):
        # Create directories
        if not os.path.exists(self.expname):
            os.makedirs(self.expname)
        if not os.path.exists(self.expname+"/mismatch"):
            os.makedirs(self.expname+"/mismatch")
        if not os.path.exists(self.expname+"/pretRNAs"):
            os.makedirs(self.expname+"/pretRNAs")
        if not os.path.exists(self.expname+"/unique"):
            os.makedirs(self.expname+"/unique")
        if not os.path.exists(self.expname+"/trna"):
            os.makedirs(self.expname+"/trna")
        if not os.path.exists(self.bamdir):
            os.makedirs(self.bamdir)

        # Expand dbname
        self.dbname = os.path.expanduser(self.dbname)
        
        if self.hubonly:
            print("Generating Track Hub...", file=sys.stderr)
            self.createtrackhub()
            return

        # Mapping Reads
        print("Mapping Reads", file=sys.stderr)
        self.mapsamples()
        
        if self.maponly:
            return

        self.makefeaturebed()
        
        # Counting Reads
        print("Counting Reads", file=sys.stderr)
        self.countfeatures()
        
        print("Analyzing counts", file=sys.stderr)
        # DESeq2 analysis using PyDESeq2
        if not self.nosizefactors:
            self.run_deseq2()
            
        # Counting Read Types
        print("Counting Read Types", file=sys.stderr)
        self.counttypes()
        
        # Coverage plots
        print("Generating Read Coverage plots", file=sys.stderr)
        orgtype = self.trnainfo.getorgtype()
        self.gettrnacoverage(orgtype)
        
        if self.makehubs:
            self.createtrackhub()

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
        mapresults = {}
        if map_args:
            with Pool(processes=pool_size) as pool:
                for result in pool.imap_unordered(map_sample_wrapper, map_args):
                    if result.failedrun:
                        print("Failure to Bowtie2 map", file=sys.stderr)
                        result.printbowtie()
                        sys.exit(1)
                    mapresults[result.samplename] = result
                    result.printbowtie()
        
        # Write mapinfo
        if self.expinfo.mapinfo:
            with open(self.expinfo.mapinfo, 'w') as mapinfo:
                print("\t".join(samples), file=mapinfo)
                print("unmap\t" + "\t".join(str(mapresults[s].unmaps) if s in mapresults else "0" for s in samples), file=mapinfo)
                print("single\t" + "\t".join(str(mapresults[s].singlemaps) if s in mapresults else "0" for s in samples), file=mapinfo)
                print("multi\t" + "\t".join(str(mapresults[s].multimaps) if s in mapresults else "0" for s in samples), file=mapinfo)
                print("total\t" + "\t".join(str(mapresults[s].totalreads) if s in mapresults else "0" for s in samples), file=mapinfo)
                print("bowtiecommand\t" + "\t".join(str(mapresults[s].bowtiecommand) if s in mapresults else "" for s in samples), file=mapinfo)
        
        # Write trna mapinfo
        if self.expinfo.trnamapfile and self.expinfo.mapinfo:
            with open(self.expinfo.trnamapfile, 'w') as trnamapinfo:
                print("\t".join(samples), file=trnamapinfo)
                print("multtrans\t" + "\t".join(str(mapresults[s].trnamapinfo.multtrans) if s in mapresults and mapresults[s].trnamapinfo is not None else "0" for s in samples), file=trnamapinfo)
                print("multac\t" + "\t".join(str(mapresults[s].trnamapinfo.multac) if s in mapresults and mapresults[s].trnamapinfo is not None else "0" for s in samples), file=trnamapinfo)
                print("multamino\t" + "\t".join(str(mapresults[s].trnamapinfo.multamino) if s in mapresults and mapresults[s].trnamapinfo is not None else "0" for s in samples), file=trnamapinfo)
                print("trna\t" + "\t".join(str(mapresults[s].trnamapinfo.trna) if s in mapresults and mapresults[s].trnamapinfo is not None else "0" for s in samples), file=trnamapinfo)
                print("singlenon\t" + "\t".join(str(mapresults[s].trnamapinfo.singlenon) if s in mapresults and mapresults[s].trnamapinfo is not None else "0" for s in samples), file=trnamapinfo)
                print("multiplenon\t" + "\t".join(str(mapresults[s].trnamapinfo.multiplenon) if s in mapresults and mapresults[s].trnamapinfo is not None else "0" for s in samples), file=trnamapinfo)
                print("total\t" + "\t".join(str(mapresults[s].trnamapinfo.uniquereads() + mapresults[s].trnamapinfo.nonuniquereads()) if s in mapresults and mapresults[s].trnamapinfo is not None else "0" for s in samples), file=trnamapinfo)

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
                            nofrag=self.nofrag, cores=self.cores, maxmismatches=self.maxmismatches)

    def run_deseq2(self):
        # Load counts
        counts_df = pd.read_csv(self.expinfo.genecounts, sep='\t', index_col=0)
        # Transpose because PyDESeq2 expects samples as rows
        counts_df = counts_df.T
        
        # Load sample info
        try:
            # Use sep=None to auto-detect delimiter (handles tabs or spaces)
            sample_df = pd.read_csv(self.samplefilename, sep=None, engine='python', header=None, names=['sample', 'condition', 'replicate'])
            sample_df.set_index('sample', inplace=True)
        except Exception as e:
            print(f"Error reading sample file {self.samplefilename}: {e}", file=sys.stderr)
            return
        
        # Filter samples that are in counts
        sample_df = sample_df.loc[counts_df.index]
        
        # Run DESeq2
        if counts_df.empty or (counts_df.sum().sum() == 0):
            print("Warning: Counts matrix is empty or all zeros. Skipping DESeq2.", file=sys.stderr)
            # Create dummy size factors
            pd.DataFrame(1.0, index=sample_df.index, columns=['sizeFactor']).T.to_csv(self.expinfo.sizefactors, sep=' ', index=False)
            # Create empty normalized counts
            counts_df.T.to_csv(self.expinfo.normalizedcounts, sep='\t')
            return

        try:
            dds = DeseqDataSet(counts=counts_df, metadata=sample_df, design_factors="condition")
            dds.deseq2()
        except Exception as e:
            print(f"Warning: DESeq2 failed: {e}", file=sys.stderr)
            # Create dummy size factors
            pd.DataFrame(1.0, index=sample_df.index, columns=['sizeFactor']).T.to_csv(self.expinfo.sizefactors, sep=' ', index=False)
            # Save raw counts as normalized counts (fallback)
            counts_df.T.to_csv(self.expinfo.normalizedcounts, sep='\t')
            return
        
        # Save normalized counts
        if 'normed_counts' in dds.layers:
            norm_counts = dds.layers['normed_counts']
            norm_counts = pd.DataFrame(norm_counts, index=dds.obs_names, columns=dds.var_names)
            norm_counts = norm_counts.T # Transpose back
            norm_counts.to_csv(self.expinfo.normalizedcounts, sep='\t')
        else:
            print("Warning: 'normed_counts' not found in dds.layers. Using raw counts.", file=sys.stderr)
            counts_df.T.to_csv(self.expinfo.normalizedcounts, sep='\t')
        
        # Save size factors
        if 'size_factors' in dds.obs:
            size_factors = dds.obs['size_factors']
            pd.DataFrame(size_factors.values, index=dds.obs_names, columns=['sizeFactor']).T.to_csv(self.expinfo.sizefactors, sep=' ', index=False)
        elif 'size_factors' in dds.obsm:
            size_factors = dds.obsm['size_factors']
            pd.DataFrame(size_factors.values, index=dds.obs_names, columns=['sizeFactor']).T.to_csv(self.expinfo.sizefactors, sep=' ', index=False)
        else:
             print("Warning: 'size_factors' not found in dds.obs or dds.obsm. Using 1.0.", file=sys.stderr)
             pd.DataFrame(1.0, index=dds.obs_names, columns=['sizeFactor']).T.to_csv(self.expinfo.sizefactors, sep=' ', index=False)

        # Run pairwise comparisons if pairs file is provided
        if self.pairfile:
            self.run_pairwise_de(dds, sample_df)

    def run_pairwise_de(self, dds, sample_df):
        try:
            # Use sep='\s+' to handle any whitespace
            pairs_df = pd.read_csv(self.pairfile, sep=r'\s+', engine='python', header=None, names=['Sample1', 'Sample2'])
        except Exception as e:
            print(f"Error reading pairs file: {e}", file=sys.stderr)
            return

        # Ensure output directory exists
        os.makedirs(os.path.join(self.expinfo.expname, "de_results"), exist_ok=True)

        for index, row in pairs_df.iterrows():
            sample1 = row['Sample1']
            sample2 = row['Sample2']
            
            # Check if samples exist in metadata or are valid conditions
            is_sample1 = sample1 in sample_df.index
            is_sample2 = sample2 in sample_df.index
            
            is_cond1 = sample1 in sample_df['condition'].values
            is_cond2 = sample2 in sample_df['condition'].values
            
            if not (is_sample1 or is_cond1) or not (is_sample2 or is_cond2):
                print(f"Warning: {sample1} or {sample2} not found in metadata (as sample or condition). Skipping pair.", file=sys.stderr)
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
                print(f"Warning: Comparison between same condition {cond1}. Skipping DE.", file=sys.stderr)
                continue

            print(f"Running DE for {cond1} vs {cond2}...", file=sys.stderr)
            
            try:
                stat_res = DeseqStats(dds, contrast=["condition", cond1, cond2])
                stat_res.summary()
                res_df = stat_res.results_df
                
                # Save results
                out_file = os.path.join(self.expinfo.expname, "de_results", f"{cond1}_vs_{cond2}.txt")
                res_df.to_csv(out_file, sep='\t')
                print(f"Saved DE results to {out_file}", file=sys.stderr)
                
            except Exception as e:
                print(f"Error running DE for {cond1} vs {cond2}: {e}", file=sys.stderr)

    def counttypes(self):
        if not self.nosizefactors:
            toolsCountReads.main(sizefactors=self.expinfo.sizefactors, combinereps=True,
                                bamdir=self.bamdir, otherseqs=self.trnainfo.otherseqs,
                                samplefile=self.samplefilename, maturetrnas=[self.trnainfo.maturetrnas],
                                trnatable=self.trnainfo.trnatable, trnaaminofile=self.expinfo.trnaaminofile,
                                trnaanticodonfile=self.expinfo.trnaanticodonfile, ensemblgtf=self.ensgtf,
                                trnaloci=[self.trnainfo.locifile], countfile=self.expinfo.genetypecounts,
                                realcountfile=self.expinfo.genetyperealcounts, mismatchfile=self.expinfo.mismatchcountfile,
                                bedfile=self.bedfiles, readlengthfile=self.expinfo.trnalengthfile,
                                countfrags=False, bamnofeature=self.dumpother,
                                uniquename=self.expinfo.uniquename, fraguniq=not self.nofrag, cores=self.cores)
        else:
            toolsCountReads.main(combinereps=True, samplefile=self.samplefilename,
                                    maturetrnas=[self.trnainfo.maturetrnas], otherseqs=self.trnainfo.otherseqs,
                                    bamdir=self.bamdir, trnatable=self.trnainfo.trnatable,
                                    trnaaminofile=self.expinfo.trnaaminofile, trnaanticodonfile=self.expinfo.trnaanticodonfile,
                                    ensemblgtf=self.ensgtf, trnaloci=[self.trnainfo.locifile],
                                    countfile=self.expinfo.genetypecounts, realcountfile=self.expinfo.genetyperealcounts,
                                    bedfile=self.bedfiles, readlengthfile=self.expinfo.trnalengthfile,
                                    countfrags=False, uniquename=self.expinfo.uniquename, cores=self.cores)

    def gettrnacoverage(self, orgtype):
        if not self.nosizefactors:
            toolsGetCoverage.main(samplefile=self.samplefilename, bedfile=[self.trnainfo.maturetrnas],
                                 locibed=[self.trnainfo.locifile], locistk=self.trnainfo.locialign,
                                 bamdir=self.bamdir, lociedgemargin=30, sizefactors=self.expinfo.sizefactors,
                                 orgtype=orgtype, locicoverage=self.expinfo.locicoveragefile,
                                 stkfile=self.trnainfo.trnaalign, numfile=self.trnainfo.trnanums,
                                 locinums=self.trnainfo.locinums, allcoverage=self.expinfo.trnacoveragefile,
                                 trnafasta=self.trnainfo.trnafasta, cores=self.cores,
                                 uniqcoverage=self.expinfo.trnauniqcoveragefile, mincoverage=self.mincoverage,
                                 uniqueonly=self.uniqueonlycov)
        else:
            toolsGetCoverage.main(samplefile=self.samplefilename, bedfile=[self.trnainfo.maturetrnas],
                                 stkfile=self.trnainfo.trnaalign, uniquename=self.expname+"/"+self.expname,
                                 orgtype=orgtype, bamdir=self.bamdir, allcoverage=self.expinfo.trnacoveragefile,
                                 trnafasta=self.trnainfo.trnafasta, cores=self.cores,
                                 uniqcoverage=self.expinfo.trnauniqcoveragefile, mincoverage=self.mincoverage,
                                 uniqueonly=self.uniqueonlycov, locibed=[], locistk=self.trnainfo.locialign)

    def createtrackhub(self):
        hub_builder = toolsTrackHub.TrackHubBuilder(
            genomedatabase=self.trnainfo.dbname,
            samplefilename=self.samplefilename,
            expname=self.expinfo.expname,
            threads=self.cores
        )
        hub_builder.run()
