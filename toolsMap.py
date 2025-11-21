#!/usr/bin/env python3

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

import toolsQC

import toolsTG
import toolsCountReads
import toolsGetCoverage
import toolsTrackHub

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
    def __init__(self, bowtiedb, trnafile, scriptdir, minnontrnasize=20, local=False, maxmaps=100, program='bowtie2'):
        self.bowtiedb = bowtiedb
        self.trnafile = trnafile
        self.scriptdir = scriptdir
        self.minnontrnasize = minnontrnasize
        self.local = local
        self.maxmaps = maxmaps
        self.program = program

    def map_sample(self, samplename, fastqfile, bamfile, expname, logfile=None):
        localmode = " "
        if self.local:
            localmode = " --local "
        
        # Construct bowtie2 command
        # Note: choosemappings.py is expected to be in scriptdir
        choosemappings_script = os.path.join(self.scriptdir, 'choosemappings.py')
        
        bowtiecommand = f"{self.program}{localmode} -x {self.bowtiedb} -k {self.maxmaps} --very-sensitive --ignore-quals --np 5 --reorder -p 1 -U {fastqfile}"
        
        temploc = os.path.basename(bamfile) + ''.join(random.choice(string.ascii_lowercase) for i in range(8))
        print(temploc, file=sys.stderr)
        
        bowtiecommand += f" | {choosemappings_script} {self.trnafile} --progname=TRAX --fqname={fastqfile} --expname={expname} --minnontrnasize={self.minnontrnasize}"
        bowtiecommand += f" | samtools sort -T {tempfile.gettempdir()}/{temploc}temp - -o {bamfile}.bam"
        
        print(bowtiecommand, file=sys.stderr)
        if logfile:
            print(bowtiecommand, file=logfile)
            logfile.flush()
            
        bowtierun = subprocess.Popen(bowtiecommand, shell=True, stderr=subprocess.PIPE, universal_newlines=True)
        output = bowtierun.communicate()
        errinfo = output[1]
        
        if logfile is not None:
            print(errinfo, file=logfile) 
            logfile.flush()
            
        if bowtierun.returncode:
            return MapInfo(0, 0, 0, 0, errinfo, samplename, failedrun=True, bowtiecommand=bowtiecommand)

        # Parse bowtie2 output
        rereadmulttrans = re.search(r'tRNA Reads with multiple transcripts:(\d+)', errinfo)
        rereadmultac = re.search(r'tRNA Reads with multiple anticodons:(\d+)', errinfo)
        rereadmultamino = re.search(r'tRNA Reads with multiple aminos:(\d+)', errinfo)
        rereadtrna = re.search(r'Total tRNA Reads:(\d+)', errinfo)
        rereadsinglenon = re.search(r'Single mapped non-tRNAs:(\d+)', errinfo)
        rereadmultiplenon = re.search(r'Multiply mapped non-tRNAs:(\d+)', errinfo)
        
        trnamapinfo = None
        if rereadmulttrans and rereadmultac and rereadmultamino and rereadtrna and rereadsinglenon and rereadmultiplenon:
             trnamapinfo = TRNAMapInfo(rereadmulttrans.group(1), rereadmultac.group(1), rereadmultamino.group(1), 
                                       rereadtrna.group(1), rereadsinglenon.group(1), rereadmultiplenon.group(1)) 

        rereadtotal = re.search(r'(\d+).*reads', errinfo)
        rereadunmap = re.search(r'\s*(\d+).*0 times', errinfo)
        rereadsingle = re.search(r'\s*(\d+).*exactly 1 time', errinfo)
        rereadmult = re.search(r'\s*(\d+).*>1 times', errinfo)
        
        if rereadtotal and rereadunmap and rereadsingle and rereadmult:
            totalreads = rereadtotal.group(1)
            unmappedreads = rereadunmap.group(1)
            singlemaps = rereadsingle.group(1)
            multmaps = rereadmult.group(1)
            return MapInfo(singlemaps, multmaps, unmappedreads, totalreads, errinfo, samplename, bowtiecommand=bowtiecommand, trnamapinfo=trnamapinfo)
        else:
            print(f"Could not map {fastqfile}, check mapstats file", file=sys.stderr)
            print("Exiting...", file=sys.stderr)
            print(errinfo, file=sys.stderr)
            return MapInfo(0, 0, 0, 0, errinfo, samplename, failedrun=True, bowtiecommand=bowtiecommand)

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
        if 'PG' in newheader and len(newheader["PG"]) > 1 and newheader["PG"][1]["PN"] == "TRAX":
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
        self.bamdir = args.bamdir if args.bamdir else "./"
        self.cores = args.threads if args.threads else 8
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

        # Expand dbname
        self.dbname = os.path.expanduser(self.dbname)
        
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
        
        self.gettraxqc()
        
        if self.makehubs:
            self.createtrackhub()

    def mapsamples(self):
        # Initialize MapReads
        scriptdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../tRAX'))
        mapper = MapReads(bowtiedb=self.trnainfo.bowtiedb, trnafile=self.trnainfo.trnatable, 
                          scriptdir=scriptdir, minnontrnasize=self.minnontrnasize, 
                          local=self.local)
        
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
            with Pool(processes=self.cores) as pool:
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

        # Write trnamapinfo
        if self.expinfo.trnamapfile:
            with open(self.expinfo.trnamapfile, 'w') as trnamapinfo:
                print("\t".join(samples), file=trnamapinfo)
                print("multi_nontRNA\t" + "\t".join(str(mapresults[s].trnamapinfo.multiplenon) if s in mapresults and mapresults[s].trnamapinfo else "0" for s in samples), file=trnamapinfo)
                print("unique_nontRNA\t" + "\t".join(str(mapresults[s].trnamapinfo.singlenon) if s in mapresults and mapresults[s].trnamapinfo else "0" for s in samples), file=trnamapinfo)
                print("multi_amino\t" + "\t".join(str(mapresults[s].trnamapinfo.multamino) if s in mapresults and mapresults[s].trnamapinfo else "0" for s in samples), file=trnamapinfo)
                print("unique_amino\t" + "\t".join(str(mapresults[s].trnamapinfo.multac) if s in mapresults and mapresults[s].trnamapinfo else "0" for s in samples), file=trnamapinfo)
                print("unique_anticodon\t" + "\t".join(str(mapresults[s].trnamapinfo.multtrans) if s in mapresults and mapresults[s].trnamapinfo else "0" for s in samples), file=trnamapinfo)
                print("unique_tRNA\t" + "\t".join(str(mapresults[s].trnamapinfo.singletrna) if s in mapresults and mapresults[s].trnamapinfo else "0" for s in samples), file=trnamapinfo)


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
        # R plotting commented out
        # runrscript(scriptdir+"/pcareadcounts.R",expinfo.normalizedcounts,samplefile,expinfo.pcaplot)
        # runrscript(scriptdir+"/pcareadcounts.R",expinfo.trnacounts,samplefile,expinfo.pcatrnaplot)

    def run_deseq2(self):
        # Load counts
        counts_df = pd.read_csv(self.expinfo.genecounts, sep='\t', index_col=0)
        # Transpose because PyDESeq2 expects samples as rows
        counts_df = counts_df.T
        
        # Load sample info
        sample_df = pd.read_csv(self.samplefilename, sep='\t', header=None, names=['sample', 'condition', 'replicate'])
        sample_df.set_index('sample', inplace=True)
        
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
        if 'size_factors' in dds.obsm:
            size_factors = dds.obsm['size_factors']
            pd.DataFrame(size_factors, index=dds.obs_names, columns=['sizeFactor']).T.to_csv(self.expinfo.sizefactors, sep=' ', index=False)
        else:
             print("Warning: 'size_factors' not found in dds.obsm. Using 1.0.", file=sys.stderr)
             pd.DataFrame(1.0, index=dds.obs_names, columns=['sizeFactor']).T.to_csv(self.expinfo.sizefactors, sep=' ', index=False)

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
        
        # R plotting commented out
        # runrscript(scriptdir+"/genefeatures.R",expinfo.genetypecounts,expinfo.genetypeplot)
        # runrscript(scriptdir+"/featuretypes.R",expinfo.trnaaminofile,expinfo.trnaaminoplot, "all")
        # ...

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
                                 uniqueonly=self.uniqueonlycov, locibed=[], locistk=self.trnainfo.locialign) # Added missing required args with defaults/empty if needed
        
        # R plotting commented out
        # runrscript(scriptdir+"/newcoverageplots.R", ...)

    def gettraxqc(self):
        toolsQC.main(samplefile=self.samplefilename, databasename=self.trnainfo.dbname,
                    experimentname=self.expinfo.expname, tgirt=self.nofrag, output=self.expinfo.qaoutputname)

    def createtrackhub(self):
        hub_builder = toolsTrackHub.TrackHubBuilder(
            genomedatabase=self.trnainfo.dbname,
            samplefilename=self.samplefilename,
            expname=self.expinfo.expname,
            threads=self.cores
        )
        hub_builder.run()
