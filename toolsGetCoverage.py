#!/usr/bin/env python3

import sys
import os
import re
import gzip
import string
import itertools
import multiprocessing
from collections import defaultdict
from typing import List, Dict, Set, Tuple, Optional, Union, Generator, Any, Iterable
import pysam
import numpy as np
from dataclasses import dataclass

# --- Helper functions and classes from trnasequtils (modernized) ---
try:
    import toolsTG
except ImportError:
    from . import toolsTG



# --- getcoverage logic (modernized) ---

bam_match = 0
bam_cins = 1
bam_cdel = 2

def cigarreflength(cigar: List[Tuple[int, int]]) -> int:
    return sum(length for op, length in cigar if op in [bam_match, bam_cdel])

def cigarreadlength(cigar: List[Tuple[int, int]]) -> int:
    return sum(length for op, length in cigar if op in [bam_match, bam_cins])

def cigarrefcoverage(cigar: List[Tuple[int, int]]) -> Generator[int, None, None]:
    nextsum = 1
    for op, length in cigar:
        if op == bam_cins:
            nextsum += length
        elif op == bam_cdel:
            for _ in range(length):
                yield 0
        elif op == bam_match:
            for _ in range(length):
                yield nextsum
                nextsum = 1

gapchars = set("-._~")

class ReadCoverage:
    def __init__(self, region: toolsTG.GenomeRange):
        self.region = region
        self.coverage = [0] * region.length()
        self.totalreads = 0
        self.length = region.length()

    def coveragealign(self, alignment: str, sizefactor: float = 1.0) -> Generator[Optional[float], None, None]:
        # Check alignment length
        align_no_gaps = "".join(c for c in alignment if c not in gapchars)
        if len(self.coverage) != len(align_no_gaps):
             # This might happen if the bed file and alignment don't match
             # For now, we'll just print a warning and return None or 0s?
             # The original code exits.
             print(f"Alignment length mismatch: Cov={len(self.coverage)}, Align={len(align_no_gaps)}", file=sys.stderr)
             return

        i = 0
        lastcoverage = 0.0
        for char in alignment:
            if char in gapchars:
                yield lastcoverage # Repeat last coverage for gaps? Original code does this.
            else:
                val = self.coverage[i] / sizefactor
                lastcoverage = val
                yield val
                i += 1

    def addread(self, read: toolsTG.GenomeRange):
        self.totalreads += 1
        if self.region.strand == "+":
            start = max(0, read.start - self.region.start)
            end = min(self.length, read.end - self.region.start)
        else:
            start = max(0, self.region.end - read.end)
            end = min(self.length, self.region.end - read.start)
        
        # Python slicing is faster than loop
        # But we need to increment.
        # self.coverage[start:end] += 1 # This doesn't work on lists
        for i in range(start, end):
            self.coverage[i] += 1

    def addbase(self, base: int):
        posbase = base - self.region.start
        if 0 <= posbase < len(self.coverage):
            self.coverage[posbase] += 1

def nasum(operands: Iterable[Any]) -> Any:
    ops = list(operands)
    if all(curr == "NA" for curr in ops):
        return "NA"
    if any(curr == "NA" for curr in ops):
        print("Trying to add incompatible alignments", file=sys.stderr)
        sys.exit(1)
    return sum(ops)

def sumsamples(allcoverages: Dict, sampledata: toolsTG.samplefile, repname: str, currfeat: toolsTG.GenomeRange, 
               trnastk: toolsTG.RnaAlignment, sizefactors: Dict[str, float]) -> Generator[float, None, None]:
    
    samples = sampledata.getrepsamples(repname)
    # Generator of generators
    generators = [allcoverages[s][currfeat.name].coveragealign(trnastk.aligned_sequences[currfeat.name], sizefactor=sizefactors[s]) for s in samples]
    
    for values in zip(*generators):
        yield nasum(values)

eukpositions = [-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,'17a',18,19,20,'20a','20b',21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14','e15','e16','e17','e18','e19',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76]
archpositions = [-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,'17a',18,19,20,'20a','20b',21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76]
bactpositions = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14','e15','e16','e17',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76]
mitopositions = [-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14','e15','e16','e17',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76]


def gettnanums(trnaalign: toolsTG.RnaAlignment, margin: int = 0, orgtype: str = "euk") -> List[str]:
    trnanum = []
    currcount = 0
    enum = 1
    gapnum = 1
    intronnum = 1
    if orgtype == "bact":
        positions = bactpositions
    elif orgtype == "arch":
        positions = archpositions
    elif orgtype == "mito":
        positions = mitopositions
    else:
        positions = eukpositions
    
    for i in range(margin):
        trnanum.append(f'head{margin - i}')
        
    for i, struct in enumerate(trnaalign.consensus):
        if currcount >= len(positions):
            trnanum.append(f'gap{gapnum}')
            gapnum += 1
            currcount += 1
        elif struct in "+=*":
             if currcount == 0 and struct == '=' and orgtype != "bact":
                 currcount = 1
                 gapnum = 1
             
             pos = positions[currcount]
             if isinstance(pos, str) and pos.startswith('e'):
                 trnanum.append(f'e{enum}')
                 enum += 1
                 currcount += 1
                 gapnum = 1
             elif pos == '-':
                 trnanum.append(f'{currcount}.gap{gapnum}')
                 gapnum += 1
                 currcount += 1
             else:
                 trnanum.append(str(pos))
                 currcount += 1
                 gapnum = 1
        else:
            if currcount < len(positions) and positions[currcount] == 38:
                trnanum.append(f'intron{intronnum}')
                intronnum += 1
            else:
                trnanum.append(f'{currcount}.gap{gapnum}')
                gapnum += 1
                
    for i in range(margin):
        trnanum.append(f'tail{i+1}')
    return trnanum

def readtrnanums(numfile: str, margin: int = 0) -> Generator[str, None, None]:
    for i in range(margin):
        yield f'head{margin - i}'
    with open(numfile) as f:
        for line in f:
            fields = line.rstrip().split("\t")
            yield fields[3]
    for i in range(margin):
        yield f'tail{i+1}'

@dataclass
class CoverageInfo:
    readcounts: Dict[str, int]
    allcoverages: Dict[str, ReadCoverage]
    readstarts: Dict[str, ReadCoverage]
    readends: Dict[str, ReadCoverage]
    multaminocoverages: Dict[str, ReadCoverage]
    multaccoverages: Dict[str, ReadCoverage]
    multtrnacoverages: Dict[str, ReadCoverage]
    uniquecoverages: Dict[str, ReadCoverage]
    uniquegenomecoverages: Dict[str, ReadCoverage]
    multigenomecoverages: Dict[str, ReadCoverage]
    readmismatches: Dict[str, ReadCoverage]
    adeninemismatches: Dict[str, ReadCoverage]
    thyminemismatches: Dict[str, ReadCoverage]
    cytosinemismatches: Dict[str, ReadCoverage]
    guanosinemismatches: Dict[str, ReadCoverage]
    readskips: Dict[str, ReadCoverage]
    trimcoverage: Dict[str, ReadCoverage]
    trimmismatches: Dict[str, ReadCoverage]

@dataclass
class LociCoverageInfo:
    readcounts: Dict[str, int]
    allcoverages: Dict[str, ReadCoverage]

def getlocicoverage(currsample: str, sampledata: toolsTG.samplefile, trnaloci: List[toolsTG.GenomeRange], 
                    maxmismatches: Optional[int] = None, minextend: Optional[int] = None) -> LociCoverageInfo:
    currbam = sampledata.getbam(currsample)
    allcoverages = {}
    readcounts = {}
    
    try:
        if not os.path.isfile(currbam + ".bai"):
            pysam.index(currbam)
        bamfile = pysam.AlignmentFile(currbam, "rb")
    except IOError as e:
        print(f"{e}", file=sys.stderr)
        sys.exit(1)
        
    for currfeat in trnaloci:
        allcoverages[currfeat.name] = ReadCoverage(currfeat)
        readcounts[currfeat.name] = 0
        
        for currread in toolsTG.getbam(bamfile, currfeat, primaryonly=False, singleonly=False, allowindels=True):
            if maxmismatches is not None:
                xm = currread.getmismatches()
                if xm is not None and xm > maxmismatches:
                    continue
            if currfeat.coverage(currread) > 10:
                if minextend is not None:
                    if not (currread.start + minextend <= currfeat.start or currread.end - minextend >= currfeat.end):
                        continue
                
                readcounts[currfeat.name] += 1
                allcoverages[currfeat.name].addread(currread)
                
    return LociCoverageInfo(readcounts, allcoverages)

def getsamplecoverage(currsample: str, sampledata: toolsTG.samplefile, trnalist: List[toolsTG.GenomeRange], 
                      trnaseqs: Dict[str, str], maxmismatches: Optional[int] = None, 
                      minextend: Optional[int] = None, removestart: bool = True, 
                      uniqueonly: bool = False) -> CoverageInfo:
    
    currbam = sampledata.getbam(currsample)
    
    # Initialize dictionaries
    allcoverages = {}
    multaminocoverages = {}
    multaccoverages = {}
    multtrnacoverages = {}
    uniquecoverages = {}
    uniquegenomecoverages = {}
    multigenomecoverages = {}
    
    readmismatches = {}
    adeninemismatches = {}
    thyminemismatches = {}
    cytosinemismatches = {}
    guanosinemismatches = {}
    readstarts = {}
    readends = {}
    readskips = {}
    
    trimreadcoverage = {}
    trimreadmismatches = {}
    readcounts = {}
    
    try:
        if not os.path.isfile(currbam + ".bai"):
            pysam.index(currbam)
        bamfile = pysam.AlignmentFile(currbam, "rb")
    except IOError as e:
        print(f"{e}", file=sys.stderr)
        sys.exit(1)
        
    for currfeat in trnalist:
        name = currfeat.name
        allcoverages[name] = ReadCoverage(currfeat)
        multaminocoverages[name] = ReadCoverage(currfeat)
        multaccoverages[name] = ReadCoverage(currfeat)
        multtrnacoverages[name] = ReadCoverage(currfeat)
        uniquecoverages[name] = ReadCoverage(currfeat)
        uniquegenomecoverages[name] = ReadCoverage(currfeat)
        multigenomecoverages[name] = ReadCoverage(currfeat)
        
        readstarts[name] = ReadCoverage(currfeat)
        readends[name] = ReadCoverage(currfeat)
        readmismatches[name] = ReadCoverage(currfeat)
        adeninemismatches[name] = ReadCoverage(currfeat)
        thyminemismatches[name] = ReadCoverage(currfeat)
        cytosinemismatches[name] = ReadCoverage(currfeat)
        guanosinemismatches[name] = ReadCoverage(currfeat)
        readskips[name] = ReadCoverage(currfeat)
        
        trimreadcoverage[name] = ReadCoverage(currfeat)
        trimreadmismatches[name] = ReadCoverage(currfeat)
        readcounts[name] = 0
        
        for currread in toolsTG.getbam(bamfile, currfeat, primaryonly=False, singleonly=False, allowindels=True):
            if maxmismatches is not None:
                xm = currread.getmismatches()
                if xm is not None and xm > maxmismatches:
                    continue
            if currfeat.coverage(currread) > 10:
                if minextend is not None:
                    if not (currread.start + minextend <= currfeat.start or currread.end - minextend >= currfeat.end):
                        continue
                
                if uniqueonly and not currread.issinglemapped():
                    continue
                    
                readstart = currread.getfirst(1)
                readend = currread.getlast(1)
                
                readcounts[name] += 1
                allcoverages[name].addread(currread)
                readstarts[name].addread(readstart)
                readends[name].addread(readend)
                
                if not currread.isuniqueaminomapping():
                    multaminocoverages[name].addread(currread)
                elif not currread.isuniqueacmapping():
                    multaccoverages[name].addread(currread)
                elif not currread.isuniquetrnamapping():
                    multtrnacoverages[name].addread(currread)
                else:
                    uniquecoverages[name].addread(currread)
                    
                if currread.issinglemapped():
                    uniquegenomecoverages[name].addread(currread)
                else:
                    multigenomecoverages[name].addread(currread)
                
                currseq = currread.getseq()
                readcov = list(cigarrefcoverage(currread.getcigar()))
                
                # Construct aligned sequence
                alignseq_chars = []
                cigar_ref_len = cigarreflength(currread.getcigar())
                
                # This logic from original code seems a bit complex/fragile regarding indices
                # "currseq[sum(readcov[0:i])] if readcov[i] > 0 else '-'"
                # Let's try to replicate it carefully.
                current_read_idx = 0
                for i in range(cigar_ref_len):
                    if i < len(readcov):
                        cov = readcov[i]
                        if cov > 0:
                            if current_read_idx < len(currseq):
                                alignseq_chars.append(currseq[current_read_idx])
                                current_read_idx += 1 # Assuming cov is 1 for match? 
                                # Wait, readcov returns nextsum which is accumulated length?
                                # No, cigarrefcoverage yields:
                                # match: nextsum (accumulated insertions), then nextsum=1
                                # del: 0
                                # ins: accumulates to nextsum
                                
                                # Let's re-read cigarrefcoverage in original code.
                                # if match: yield nextsum; nextsum=1
                                # if del: yield 0
                                # if ins: nextsum += len
                                
                                # So readcov[i] is the number of bases in read corresponding to ref pos i?
                                # If > 0, it means there is a base. If > 1, it means insertion + base?
                                # The original code: currseq[sum(readcov[0:i])]
                                # sum(readcov[0:i]) is the index in read.
                                pass
                            else:
                                alignseq_chars.append("N") # Should not happen
                        else:
                            alignseq_chars.append("-")
                    else:
                         alignseq_chars.append("-")

                # Re-implementing alignseq construction properly based on original logic
                # Original: alignseq = "".join(currseq[sum(readcov[0:i])] if readcov[i] > 0 else "-" for i in range(cigarreflength(currread.getcigar())))
                
                # Optimization: calculate cumulative sum once
                readcov_cumsum = [0]
                current_sum = 0
                for x in readcov:
                    current_sum += x
                    readcov_cumsum.append(current_sum)
                
                alignseq_list = []
                for i in range(cigar_ref_len):
                    if i < len(readcov) and readcov[i] > 0:
                        idx = readcov_cumsum[i]
                        if idx < len(currseq):
                            alignseq_list.append(currseq[idx])
                        else:
                            alignseq_list.append("N")
                    else:
                        alignseq_list.append("-")
                alignseq = "".join(alignseq_list)

                refseq = trnaseqs[name][currread.start : currread.start + cigar_ref_len]
                
                for currpos in range(cigar_ref_len):
                    if currpos >= len(alignseq) or currpos >= len(refseq):
                        break
                        
                    currbase = alignseq[currpos]
                    refbase = refseq[currpos]
                    
                    if refbase not in gapchars:
                        if currbase == "-":
                            readskips[name].addbase(currread.start + currpos)
                        
                        if refbase != currbase:
                            readmismatches[name].addbase(currread.start + currpos)
                            
                        if currpos > 3:
                            trimreadcoverage[name].addbase(currread.start + currpos)
                            if refbase != currbase:
                                trimreadmismatches[name].addbase(currread.start + currpos)
                                
                        if currpos > 3 and removestart:
                            if currbase == "A":
                                adeninemismatches[name].addbase(currread.start + currpos)
                            elif currbase == "T":
                                thyminemismatches[name].addbase(currread.start + currpos)
                            elif currbase == "C":
                                cytosinemismatches[name].addbase(currread.start + currpos)
                            elif currbase == "G":
                                guanosinemismatches[name].addbase(currread.start + currpos)

    return CoverageInfo(readcounts, allcoverages, readstarts, readends, multaminocoverages, 
                        multaccoverages, multtrnacoverages, uniquecoverages, uniquegenomecoverages, 
                        multigenomecoverages, readmismatches, adeninemismatches, thyminemismatches, 
                        cytosinemismatches, guanosinemismatches, readskips, trimreadcoverage, trimreadmismatches)

def transcriptcoverage(samplecoverages: Dict[str, CoverageInfo], mismatchreport: Any, 
                       trnalist: List[toolsTG.GenomeRange], sampledata: toolsTG.samplefile, sizefactor: Dict[str, float], 
                       mincoverage: int, trnastk: toolsTG.RnaAlignment, positionnums: List[str], 
                       skipgaps: bool = True, sigmismatch: Any = None):
    
    samples = sampledata.getsamples()
    mismatchthreshold = 0.05

    for currfeat in trnalist:
        name = currfeat.name
        totalreads = sum(samplecoverages[s].allcoverages[name].totalreads for s in samples)
        ambigreads = sum(samplecoverages[s].multaminocoverages[name].totalreads for s in samples)
        
        if totalreads - ambigreads < mincoverage:
            continue
            
        if sigmismatch:
            mismatchpos = {}
            coveragepos = {}
            trnalen = 0
            
            for currsample in samples:
                sc = samplecoverages[currsample]
                align = trnastk.aligned_sequences[name]
                
                # Use raw counts (sizefactor=1) for percentage calculation
                covcounts_raw = list(sc.allcoverages[name].coveragealign(align, sizefactor=1))
                mismatches_raw = list(sc.readmismatches[name].coveragealign(align, sizefactor=1))
                
                trnalen = len(mismatches_raw)
                mismatchpos[currsample] = mismatches_raw
                coveragepos[currsample] = covcounts_raw
            
            for currpos in range(trnalen):
                if skipgaps and "gap" in str(positionnums[currpos]):
                    continue
                
                maxpercent = max((0 + 1.0 * mismatchpos[s][currpos]) / (10 + coveragepos[s][currpos]) for s in samples)
                
                if maxpercent > mismatchthreshold:
                    try:
                        # Note: currpos is alignment index. If alignment has gaps, this might be off for genomic position.
                        # However, tRAX uses the same logic.
                        print(currfeat.getbase(currpos).bedstring(name = currfeat.name+"_"+str(currpos)+"pos", score = int(maxpercent * 1000)), file=sigmismatch)
                    except Exception:
                        pass

        for currsample in samples:
            sf = sizefactor[currsample]
            sc = samplecoverages[currsample]
            align = trnastk.aligned_sequences[name]
            
            covcounts = list(sc.allcoverages[name].coveragealign(align, sizefactor=sf))
            mismatches = list(sc.readmismatches[name].coveragealign(align, sizefactor=sf))
            deletions = list(sc.readskips[name].coveragealign(align, sizefactor=sf))
            uniquecounts = list(sc.uniquecoverages[name].coveragealign(align, sizefactor=sf))
            multitrna = list(sc.multtrnacoverages[name].coveragealign(align, sizefactor=sf))
            multaccounts = list(sc.multaccoverages[name].coveragealign(align, sizefactor=sf))
            multaminocounts = list(sc.multaminocoverages[name].coveragealign(align, sizefactor=sf))
            allstarts = list(sc.readstarts[name].coveragealign(align, sizefactor=sf))
            allends = list(sc.readends[name].coveragealign(align, sizefactor=sf))
            allcovcount = list(sc.allcoverages[name].coveragealign(align, sizefactor=sf))
            adeninecount = list(sc.adeninemismatches[name].coveragealign(align, sizefactor=1)) # Ints in original
            thyminecount = list(sc.thyminemismatches[name].coveragealign(align, sizefactor=1))
            cytosinecount = list(sc.cytosinemismatches[name].coveragealign(align, sizefactor=1))
            guanosinecount = list(sc.guanosinemismatches[name].coveragealign(align, sizefactor=1))
            readskipcount = list(sc.readskips[name].coveragealign(align, sizefactor=1))
            
            for i, currcount in enumerate(allcovcount):
                if skipgaps and "gap" in str(positionnums[i]):
                    continue
                
                realbase = align[i].upper()
                if realbase in gapchars:
                    realbase = "-"
                if realbase == "U":
                    realbase = "T"
                
                row = [
                    name, currsample, str(positionnums[i]), str(covcounts[i]), str(allstarts[i]), 
                    str(allends[i]), str(uniquecounts[i]), str(multitrna[i]), str(multaccounts[i]), 
                    str(multaminocounts[i]), str(sc.readcounts[name]/sf), realbase, str(mismatches[i]), 
                    str(deletions[i]), str(int(adeninecount[i])), str(int(thyminecount[i])), str(int(cytosinecount[i])), 
                    str(int(guanosinecount[i])), str(int(readskipcount[i]))
                ]
                print("\t".join(row), file=mismatchreport)

def locuscoverage(locicoverages: Dict[str, LociCoverageInfo], locicoveragetable: Any, 
                  locilist: List[toolsTG.GenomeRange], sampledata: toolsTG.samplefile, sizefactor: Dict[str, float], 
                  mincoverage: int, locistk: toolsTG.RnaAlignment, locipositionnums: List[str], 
                  skipgaps: bool = True):
    
    samples = sampledata.getsamples()
    for currfeat in locilist:
        name = currfeat.name
        totalreads = sum(locicoverages[s].allcoverages[name].totalreads for s in samples)
        
        if totalreads < mincoverage:
            continue
            
        for currsample in samples:
            sf = sizefactor[currsample]
            lc = locicoverages[currsample]
            align = locistk.aligned_sequences[name]
            
            allcovcount = list(lc.allcoverages[name].coveragealign(align, sizefactor=sf))
            
            for i, currcount in enumerate(allcovcount):
                if skipgaps and "gap" in str(locipositionnums[i]):
                    continue
                
                row = [name, currsample, str(locipositionnums[i]), str(currcount), str(lc.readcounts[name])]
                print("\t".join(row), file=locicoveragetable)

# Wrapper functions for multiprocessing
def makecoveragepool(args):
    return getsamplecoverage(*args[0], **args[1])

def makelocicoveragepool(args):
    return getlocicoverage(*args[0], **args[1])

def compressargs(*args, **kwargs):
    return (args, kwargs)

def main(samplefile: str, bedfile: List[str], stkfile: str, 
         bamdir: str = "./", edgemargin: int = 0, mincoverage: Optional[int] = 30, 
         sizefactors: Optional[str] = None, orgtype: str = "euk", 
         locicoverage: str = "locicoverage.txt", numfile: Optional[str] = None, 
         locinums: Optional[str] = None, allcoverage: str = "stdout", 
         trnafasta: Optional[str] = None, cores: int = 1, 
         uniqcoverage: Optional[str] = None, uniqueonly: bool = False, 
         maxmismatches: Optional[int] = None, minextend: Optional[int] = None, 
         combinereps: bool = False, uniquename: Optional[str] = None, 
         uniquegenome: Optional[str] = None, lociedgemargin: int = 30,
         locibed: Optional[List[str]] = None, locistk: Optional[str] = None,
         sigmismatch: Optional[str] = None):
    

    if mincoverage is None:
        mincoverage = 30

    sampledata = toolsTG.samplefile(samplefile, bamdir=bamdir)
    
    if trnafasta:
        trnaseqs = {name: seq for name, seq in toolsTG.read_multi_fasta(trnafasta)}
        for currname in trnaseqs:
            trnaseqs[currname] = ("N"*edgemargin) + trnaseqs[currname] + ("N"*edgemargin)
    else:
        trnaseqs = {} 

    sizefactor_dict = defaultdict(lambda: 1.0)
    if sizefactors:
        sizefactor_dict.update(toolsTG.getsizefactors(sizefactors))
        for currsample in sampledata.getsamples():
            if currsample not in sizefactor_dict:
                print(f"Size factor file {sizefactors} missing {currsample}", file=sys.stderr)
                sys.exit(1)

    with open(stkfile, "r") as f:
        trnastk = list(toolsTG.read_rna_stk(f))[0]
    
    locistk_obj = None
    if locistk:
        with open(locistk, "r") as f:
            locistk_obj = list(toolsTG.read_rna_stk(f))[0]
        locistk_obj = locistk_obj.add_margin(lociedgemargin)

    if orgtype != "euk" and numfile and os.path.isfile(numfile):
        positionnums = list(readtrnanums(numfile, margin=edgemargin))
        locipositionnums = list(readtrnanums(locinums, margin=lociedgemargin)) if locinums else []
    else:
        positionnums = gettnanums(trnastk, margin=edgemargin, orgtype=orgtype)
        locipositionnums = gettnanums(locistk_obj, margin=lociedgemargin, orgtype=orgtype) if locistk_obj else []

    trnastk = trnastk.add_margin(edgemargin)

    basetrnas = []
    for currfile in bedfile:
        basetrnas.extend(list(toolsTG.readbed(currfile)))
    
    locitrnas = []
    if locibed:
        for currfile in locibed:
            locitrnas.extend(list(toolsTG.readbed(currfile)))

    trnalist = [curr.addmargin(edgemargin) for curr in basetrnas]
    locilist = [curr.addmargin(lociedgemargin) for curr in locitrnas]

    if allcoverage == "stdout":
        coveragetable = sys.stdout
    else:
        coveragetable = open(allcoverage, "w")
    
    print("\t".join(["Feature","Sample","position","coverage","readstarts","readends","uniquecoverage",
                     "multitrnacoverage","multianticodoncoverage","multiaminocoverage","tRNAreadstotal",
                     "actualbase","mismatchedbases","deletedbases","adenines","thymines","cytosines",
                     "guanines","deletions"]), file=coveragetable)

    locicoveragetable_file = open(locicoverage, "w")
    print("\t".join(["tRNA_name","sample","position","coverage", "total"]), file=locicoveragetable_file)
    
    mismatchcomparetable = open("mismatchcompare.txt", "w")
    print("\t".join(["pos","firsample","secsample","firmismatches","firtotal","secmismatches","sectotal",
                     "firmismatchestrim","firtotaltrim","secmismatchestrim","sectotaltrim"]), file=mismatchcomparetable)

    sigmismatch_file = None
    if sigmismatch:
        sigmismatch_file = open(sigmismatch, "w")

    samples = sampledata.getsamples()
    
    # Chunking
    trnachunksize = 50
    
    # Multiprocessing setup
    if cores > 1:
        pool = multiprocessing.Pool(processes=cores)
    
    for trnanum in range(0, len(trnalist), trnachunksize):
        endpoint = min(len(trnalist), trnanum + trnachunksize)
        trnasubset = trnalist[trnanum:endpoint]
        locisubset = locilist[trnanum:endpoint] if trnanum < len(locilist) else []
        
        samplecoverages = {}
        locicoverages = {}
        
        if cores == 1:
            for currsample in samples:
                samplecoverages[currsample] = getsamplecoverage(currsample, sampledata, trnasubset, trnaseqs, 
                                                                maxmismatches=maxmismatches, minextend=minextend, 
                                                                uniqueonly=uniqueonly)
                if locisubset:
                    locicoverages[currsample] = getlocicoverage(currsample, sampledata, locisubset, 
                                                            maxmismatches=maxmismatches, minextend=minextend)
        else:
            trackargs = []
            lociargs = []
            for currsample in samples:
                trackargs.append(compressargs(currsample, sampledata, trnasubset, trnaseqs, 
                                              maxmismatches=maxmismatches, minextend=minextend, uniqueonly=uniqueonly))
                if locisubset:
                    lociargs.append(compressargs(currsample, sampledata, locisubset, 
                                             maxmismatches=maxmismatches, minextend=minextend))
            
            results = pool.map(makecoveragepool, trackargs)
            for i, currsample in enumerate(samples):
                samplecoverages[currsample] = results[i]
            
            if locisubset:
                lociresults = pool.map(makelocicoveragepool, lociargs)
                for i, currsample in enumerate(samples):
                    locicoverages[currsample] = lociresults[i]

        transcriptcoverage(samplecoverages, coveragetable, trnasubset, sampledata, sizefactor_dict, 
                           mincoverage, trnastk, positionnums, sigmismatch=sigmismatch_file)
        
        if locistk_obj and locisubset:
            locuscoverage(locicoverages, locicoveragetable_file, locisubset, sampledata, sizefactor_dict, 
                      mincoverage, locistk_obj, locipositionnums)

    if cores > 1:
        pool.close()
        pool.join()

    if allcoverage != "stdout":
        coveragetable.close()
    locicoveragetable_file.close()
    mismatchcomparetable.close()
    if sigmismatch_file:
        sigmismatch_file.close()



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Calculate tRNA coverage')
    parser.add_argument('-i', '--samplefile', required=True, help='Sample file')
    parser.add_argument('-b', '--bedfile', required=True, nargs='+', help='Bed file(s)')
    parser.add_argument('-s', '--stkfile', required=True, help='Stockholm alignment file')
    parser.add_argument('--bamdir', default='./', help='Directory containing BAM files')
    parser.add_argument('--mincoverage', type=int, default=30, help='Minimum coverage')
    parser.add_argument('--sizefactors', help='Size factors file')
    parser.add_argument('--orgtype', default='euk', help='Organism type')
    parser.add_argument('--locicoverage', default='locicoverage.txt', help='Loci coverage output file')
    parser.add_argument('--allcoverage', default='stdout', help='All coverage output file')
    parser.add_argument('--trnafasta', help='tRNA fasta file')
    parser.add_argument('--locibed', nargs='+', help='Loci bed file(s)')
    parser.add_argument('--locistk', help='Loci stockholm file')
    parser.add_argument('--cores', type=int, default=1, help='Number of cores')
    parser.add_argument('--uniqueonly', action='store_true', help='Only use unique reads')
    
    args = parser.parse_args()
    
    main(samplefile=args.samplefile, bedfile=args.bedfile, stkfile=args.stkfile, 
         bamdir=args.bamdir, mincoverage=args.mincoverage, sizefactors=args.sizefactors, 
         orgtype=args.orgtype, locicoverage=args.locicoverage, allcoverage=args.allcoverage, 
         trnafasta=args.trnafasta, locibed=args.locibed, locistk=args.locistk, cores=args.cores,
         uniqueonly=args.uniqueonly)
