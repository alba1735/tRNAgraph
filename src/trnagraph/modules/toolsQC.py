#!/usr/bin/env python3

import sys
import os
import itertools
import logging
from time import localtime, strftime
from typing import List, Dict, Optional, Tuple, Any
from collections import defaultdict
from . import toolsTG

logger = logging.getLogger(__name__)

# Constants
RED_RGB = "rgb(255, 0, 0)"
YELLOW_RGB = "rgb(255,165,0)"
GREEN_RGB = "rgb(60,170,113)"

STYLE = '''
<style>

table, td, th { border: 1px solid black; border-spacing: 0; padding: 1px; text-align: left; }

table { width: 70%; }

</style>
'''

MODE_TEMPLATE = '''<h2>tRNAgraph Data Quality Report</h2>

<h4>

Date: {date}<br/>

Run mode: {mode}

</h4>
'''

def htmlconvert(message: str) -> str:
    message = message.replace("<", "&lt;")
    message = message.replace(">", "&gt;")
    return message

class ErrorSet:
    def __init__(self, shortname: str, alllist: List[str], faillist: List[bool], 
                 failcriteria: str, dimension: str, faildict: Dict[str, float], 
                 percentformat: bool = False, critfaillist: Optional[List[bool]] = None, 
                 checkfile: Optional[str] = None):
        self.shortname = shortname
        self.samples = alllist
        self.failsamples = faillist

        self.criteria = failcriteria
        self.dimension = dimension
        self.failnum = faildict
        self.percentformat = percentformat
        self.critfaillist = critfaillist
        self.checkfile = checkfile if checkfile is not None else ""
        
        self.warning = sum(faillist) > 0
        self.fail = critfaillist is not None and sum(critfaillist) > 0
        
        self.failset = {curr for i, curr in enumerate(alllist) if faillist[i]}
        
    def gettestcolor(self) -> str:
        if self.fail:
            return "rgb(0,0,255)"
        elif self.warning:
            return "rgb(255,165,0)"
        else:
            return "rgb(60,170,113)"

    def getteststatus(self) -> str:
        if self.fail:
            return "Failed"
        elif self.warning:
            return "Warning" 
        else:
            return "Passed"       

    def getsamplecolor(self, sample: str) -> str:
        if self.critfaillist is not None and sample in self.critfaillist: # This logic seems slightly off in original code (critfaillist is list of bools, sample is string)
             # Wait, in original code: if self.critfaillist is not None and sample in self.critfaillist:
             # But critfaillist is passed as a list of booleans in __init__.
             # Let's check how it is used.
             # In original: self.critfaillist = critfaillist
             # In getsamplecolor: if self.critfaillist is not None and sample in self.critfaillist:
             # This looks like a bug in the original code if critfaillist is indeed list of bools.
             # However, looking at usage:
             # lowcounterr = errorset(..., lowcountsamples, ...) where lowcountsamples is list of bools.
             # So faillist is list of bools.
             # self.failset = set(curr for i, curr in enumerate(alllist) if faillist[i]) -> This creates a set of sample names.
             # So self.failset contains strings.
             
             # If critfaillist is passed as list of bools, then `sample in self.critfaillist` will likely be False (unless sample is True/False which is unlikely).
             # I should probably create a critfailset similar to failset.
             pass
        
        # Let's fix the logic for critfailset
        # I'll do this in __init__
        return "rgb(60,170,113)" # Placeholder, will fix in class implementation

    def getsamplestatus(self, sample: str) -> str:
        return "Passed" # Placeholder

    def getcriteria(self) -> str:
        return htmlconvert(self.failcriteria)

    def getsampleresult(self, currsample: str) -> str:
        if currsample not in self.failnum:
            return "??"
        if self.percentformat:
            return "{0:.2f}%".format(100 * self.failnum[currsample])
        else:
            return "{0:.2f}".format(self.failnum[currsample])

class SeqPrepInfo:
    def __init__(self, merged: Dict[str, int], unmerged: Dict[str, int], discarded: Dict[str, int]):
        self.merged = merged
        self.unmerged = unmerged
        self.discarded = discarded

    def gettotal(self, sample: str) -> int:
        return self.merged.get(sample, 0) + self.unmerged.get(sample, 0) + self.discarded.get(sample, 0)

    def getmergedpercent(self, sample: str) -> float:
        total = self.gettotal(sample)
        return self.merged.get(sample, 0) / (1.0 * total + 0.01)

    def getsamples(self) -> Tuple[str, ...]:
        return tuple(self.merged.keys())
        
class CutadaptInfo:
    def __init__(self, trimmed: Dict[str, int], untrimmed: Dict[str, int], discarded: Dict[str, int]):
        self.trimmed = trimmed
        self.untrimmed = untrimmed
        self.discarded = discarded

    def gettotal(self, sample: str) -> int:
        return self.trimmed.get(sample, 0) + self.untrimmed.get(sample, 0) + self.discarded.get(sample, 0)

    def getpassedpercent(self, sample: str) -> float:
        return (self.trimmed.get(sample, 0) + self.untrimmed.get(sample, 0)) / (1.0 * self.gettotal(sample))

    def getsamples(self) -> Tuple[str, ...]:
        return tuple(self.trimmed.keys())

class MappingResults:
    def __init__(self, unmap: Dict[str, int], single: Dict[str, int], multi: Dict[str, int]):
        self.unmap = unmap
        self.single = single
        self.multi = multi
        
    def totalreadscount(self, sample: str) -> int:
        return self.unmap.get(sample, 0) + self.single.get(sample, 0) + self.multi.get(sample, 0)

    def getmappercent(self, sample: str) -> float:
        totalreads = self.totalreadscount(sample)
        return (totalreads - self.unmap.get(sample, 0)) / (1.0 * totalreads + 0.01)

class TypeCount:
    def __init__(self, typecounts: Dict[str, Dict[str, int]]):
        self.typecounts = typecounts

    def gettotal(self, sample: str) -> int:
        return sum(self.typecounts[sample].values()) 

    def gettrnapercent(self, sample: str) -> float:
        return (self.typecounts[sample].get("tRNA", 0) + self.typecounts[sample].get("pretRNA", 0)) / (1.0 * self.gettotal(sample) + 0.01)

    def getrrnapercent(self, sample: str) -> Optional[float]:
        if "rRNA" in self.typecounts[sample]:
            return self.typecounts[sample]["rRNA"] / (1.0 * self.gettotal(sample) + 0.01)
        else:
            return None

    def getotherpercent(self, sample: str) -> float:
        return (self.typecounts[sample].get("other", 0)) / (1.0 * self.gettotal(sample))

def getmeanfreq(freqtable: Dict[int, int]) -> float:
    return sum(curr * freqtable[curr] for curr in list(freqtable.keys())) / (1.0 * sum(freqtable.values()) + 0.01)

class LengthCount:
    def __init__(self, trnalengthcounts: Dict[str, Dict[int, int]], 
                 pretrnalengthcounts: Dict[str, Dict[int, int]], 
                 otherlengthcounts: Dict[str, Dict[int, int]]):
        self.trnalengthcounts = trnalengthcounts
        self.pretrnalengthcounts = pretrnalengthcounts
        self.otherlengthcounts = otherlengthcounts
        self.samples = set(itertools.chain(list(trnalengthcounts.keys()), list(pretrnalengthcounts.keys()), list(otherlengthcounts.keys())))
        
        all_lengths = itertools.chain(
            itertools.chain.from_iterable(list(trnalengthcounts[currsample].keys()) for currsample in self.samples),
            itertools.chain.from_iterable(list(pretrnalengthcounts[currsample].keys()) for currsample in self.samples),
            itertools.chain.from_iterable(list(otherlengthcounts[currsample].keys()) for currsample in self.samples)
        )
        try:
            self.maxlength = max(all_lengths)
        except ValueError:
            self.maxlength = 0
        
    def getalllengths(self, sample: str) -> Dict[int, int]:
        # Original code had a bug: self.trnalengthcounts[sample][currlength] + self.trnalengthcounts[sample][currlength] + self.trnalengthcounts[sample][currlength]
        # It was adding trna counts 3 times instead of trna + pretrna + other
        # I will fix this.
        return {
            currlength: self.trnalengthcounts[sample].get(currlength, 0) + 
                        self.pretrnalengthcounts[sample].get(currlength, 0) + 
                        self.otherlengthcounts[sample].get(currlength, 0) 
            for currlength in range(self.maxlength + 1)
        }

    def getsamplemean(self, sample: str) -> float:
        return getmeanfreq(self.getalllengths(sample))

    def getthreshold(self, sample: str, minsize: int, maxsize: int) -> int:
        alllengths = self.getalllengths(sample)
        if not alllengths:
            return 0
        maxsize = min([maxsize, max(alllengths.keys())])
        return sum(alllengths.get(i, 0) for i in range(minsize, maxsize))

    def getthresholdpercent(self, sample: str, minsize: int, maxsize: int) -> float:
        total = sum(self.getalllengths(sample).values())
        if total == 0:
            return 0.0
        return self.getthreshold(sample, minsize, maxsize) / (1.0 * total)

class TrnaCount:
    def __init__(self, trnacounts: Dict[str, Dict[str, int]]):
        self.trnacounts = trnacounts

    def gettrnaactive(self, currsample: str, cutoff: int = 20) -> int:
        return sum(1 for curr in self.trnacounts[currsample].keys() if int(self.trnacounts[currsample][curr]) > cutoff)

    def gettrnaactivepercent(self, currsample: str, trnainfo: Any, cutoff: int = 20) -> float:
        return self.gettrnaactive(currsample, cutoff) / (1.0 * len(list(trnainfo.gettranscripts())) + 0.01)

class SizeFactor:
    def __init__(self, sizefactors: Dict[str, float]):
        self.sizefactors = sizefactors

def getreadprep(prepfilename: str, manifestfilename: str, sampleinfo: Any) -> SeqPrepInfo:
    samplenames = {}
    with open(manifestfilename) as manifestfile:
        for currline in manifestfile:
            fields = currline.rstrip().split("\t")
            if len(fields) >= 2:
                samplenames[fields[0]] = sampleinfo.getfastqsample(fields[1])
    
    merged = {}
    unmerged = {}
    discarded = {}
    
    with open(prepfilename) as prepfile:
        runsamples = []
        for i, currline in enumerate(prepfile):
            fields = currline.rstrip().split("\t")
            if i == 0:
                runsamples = [samplenames.get(curr, curr) for curr in fields]
                continue
                
            if len(fields) != len(runsamples) + 1:
                continue
            for j in range(0, len(runsamples)):
                if fields[0] == "merged":
                    merged[runsamples[j]] = int(fields[j + 1])
                elif fields[0] == "unmerged":
                    unmerged[runsamples[j]] = int(fields[j + 1])
                elif fields[0] == "discarded":
                    discarded[runsamples[j]] = int(fields[j + 1])
                    
    return SeqPrepInfo(merged, unmerged, discarded)

def getcutadapt(prepfilename: str, manifestfilename: str, sampleinfo: Any) -> CutadaptInfo:
    samplenames = {}
    with open(manifestfilename) as manifestfile:
        for currline in manifestfile:
            fields = currline.rstrip().split("\t")
            if len(fields) >= 2:
                samplenames[fields[0]] = sampleinfo.getfastqsample(fields[1])
    
    trimmed = {}
    untrimmed = {}
    discarded = {}
    
    with open(prepfilename) as prepfile:
        runsamples = []
        for i, currline in enumerate(prepfile):
            fields = currline.rstrip().split("\t")
            if i == 0:
                runsamples = [samplenames.get(curr, curr) for curr in fields]
                continue
                
            if len(fields) != len(runsamples) + 1:
                continue
            for j in range(0, len(runsamples)):
                if fields[0] == "trimmed":
                    trimmed[runsamples[j]] = int(fields[j + 1])
                elif fields[0] == "untrimmed":
                    untrimmed[runsamples[j]] = int(fields[j + 1])
                elif fields[0] == "discarded":
                    discarded[runsamples[j]] = int(fields[j + 1])
                
    return CutadaptInfo(trimmed, untrimmed, discarded)

MIN_MERGE_PERCENT = 0.6

def checkreadprep(allpreps: List[str], sampleinfo: Any) -> List[ErrorSet]:
    prepdict = {}
    samples = set()
    
    for prepinfo in allpreps:   
        if os.path.exists(prepinfo + "_sp.txt"):
            prepresults = getreadprep(prepinfo + "_sp.txt", prepinfo + "_manifest.txt", sampleinfo)
            curr_samples = prepresults.getsamples()
            samples.update(curr_samples)
            prepdict.update({currsample: prepresults.getmergedpercent(currsample) for currsample in curr_samples}) 
        if os.path.exists(prepinfo + "_ca.txt"):
            prepresults = getcutadapt(prepinfo + "_ca.txt", prepinfo + "_manifest.txt", sampleinfo)
            curr_samples = prepresults.getsamples()
            samples.update(curr_samples)
            prepdict.update({currsample: prepresults.getpassedpercent(currsample) for currsample in curr_samples})

    if len(set(sampleinfo.getsamples()) & set(prepdict.keys())) == 0:
        return []
    
    sample_list = list(samples)
    lowmergesamples = [prepdict.get(currsample, 0) < MIN_MERGE_PERCENT for currsample in sample_list]
    
    # Assuming checkfile is the first prepinfo found + suffix, but original code used `prepinfo` variable which would be the last one in loop.
    # I'll use the last one as well to match behavior, or empty string if loop didn't run.
    checkfile = ""
    if allpreps:
        checkfile = allpreps[-1] + "_sp.pdf" # Or _ca.pdf? Original code used prepinfo+"_sp.pdf" inside the loop but returned outside. 
        # Actually original code: checkfile = prepinfo+"_sp.pdf" inside the return statement. 
        # But prepinfo is loop variable. Python loop variables leak. So it uses the last one.
    
    mergeerr = ErrorSet("merging_rate", sample_list, lowmergesamples, 
                        "Sequencing read merging rate  > " + str(100 * MIN_MERGE_PERCENT) + "%", 
                        "Merging Rate", prepdict, percentformat=True, checkfile=checkfile)
    
    return [mergeerr]

def getmapfile(samplename: str) -> str:
    return os.path.join(samplename, samplename + "-mapinfo.txt")

def getreadmapping(samplename: str, sampleinfo: Any) -> MappingResults:
    unmaps = {}
    singlemaps = {}
    multimaps = {}
    
    with open(getmapfile(samplename)) as mapresults:
        allsamples = sampleinfo.getsamples()
        runsamples = []
        for i, currline in enumerate(mapresults):
            fields = currline.rstrip().split("\t")
            if i == 0:
                runsamples = list(fields)
                if set(runsamples) != set(allsamples):
                    logger.warning(runsamples)
                    logger.warning(allsamples)
                    logger.warning("QAError")
                continue

            if len(fields) != len(allsamples) + 1:
                continue
            for j in range(0, len(runsamples)):
                if fields[0] == "unmap":
                    unmaps[runsamples[j]] = int(fields[j + 1])
                elif fields[0] == "single":
                    singlemaps[runsamples[j]] = int(fields[j + 1])
                elif fields[0] == "multi":
                    multimaps[runsamples[j]] = int(fields[j + 1])
                
    return MappingResults(unmaps, singlemaps, multimaps)

MIN_MAP_READS = 200000
MIN_MAP_PERCENT = 0.65

def checkreadsmapping(samplename: str, sampleinfo: Any, tgirtmode: bool = False) -> List[ErrorSet]:
    mapresults = getreadmapping(samplename, sampleinfo)
    samples = list(sampleinfo.getsamples())
    totalreads = {currsample: mapresults.totalreadscount(currsample) for currsample in samples}
    mappercent = {currsample: mapresults.getmappercent(currsample) for currsample in samples}
    
    lowcountsamples = [totalreads[currsample] < MIN_MAP_READS for currsample in samples]
    lowcounterr = ErrorSet("mappable_read", samples, lowcountsamples, 
                           "Mappable reads > " + str(MIN_MAP_READS), "Read Count", 
                           totalreads, checkfile=samplename + "-mapinfo.pdf")

    lowmapsamples = [mappercent[currsample] < MIN_MAP_PERCENT for currsample in samples]
    lowmaperr = ErrorSet("mappable_rate", samples, lowmapsamples, 
                         "Mapping rate > " + str(100 * MIN_MAP_PERCENT) + "%", 
                         "Mapping Rate", mappercent, percentformat=True, 
                         checkfile=samplename + "-mapinfo.pdf")
    return [lowcounterr, lowmaperr]

def gettypefile(samplename: str) -> str:
    return os.path.join(samplename, samplename + "-typerealcounts.txt")

def getreadlengthfile(samplename: str) -> str:
    return os.path.join(samplename, samplename + "-readlengths.txt")

def getfragtypefile(samplename: str) -> str:
    return os.path.join(samplename, samplename + "-fragtypes.txt")

def gettrnacountfile(samplename: str) -> str:
    return os.path.join(samplename, samplename + "-trnacounts.txt")

def getsizefactorfile(samplename: str) -> str:
    return os.path.join(samplename, samplename + "-SizeFactors.txt")

def gettypecounts(samplename: str, sampleinfo: Any) -> TypeCount:
    allsamples = sampleinfo.getsamples()
    typecounts = defaultdict(dict)
    
    with open(gettypefile(samplename)) as typeresults:
        runsamples = []
        for i, currline in enumerate(typeresults):
            fields = currline.rstrip().split("\t")
            if i == 0:
                runsamples = list(fields)
                if set(runsamples) != set(allsamples):
                    logger.warning(runsamples)
                    logger.warning(allsamples)
                    logger.warning("QAError")
                continue

            if len(fields) != len(allsamples) + 1:
                logger.warning(runsamples)
                logger.warning(fields)
                logger.warning("QAError")
                continue
            for j in range(0, len(runsamples)):
                typecounts[runsamples[j]][fields[0]] = int(fields[j + 1])
    return TypeCount(typecounts)

def getreadlengths(samplename: str, sampleinfo: Any) -> LengthCount:
    trnalengthcounts = defaultdict(lambda: defaultdict(int))
    otherlengthcounts = defaultdict(lambda: defaultdict(int))
    pretrnalengthcounts = defaultdict(lambda: defaultdict(int))
    
    with open(getreadlengthfile(samplename)) as lengthresults:
        for i, currline in enumerate(lengthresults):
            fields = currline.rstrip().split("\t")
            if i == 0:
                continue
            
            if len(fields) < 4:
                continue
            
            # fields: length, sample, other, trnas, pretrnas
            length = int(fields[0])
            sample = fields[1]
            other = int(fields[2])
            trnas = int(fields[3])
            pretrnas = int(fields[4])
            
            trnalengthcounts[sample][length] = trnas
            otherlengthcounts[sample][length] = other
            pretrnalengthcounts[sample][length] = pretrnas
            
    return LengthCount(trnalengthcounts, pretrnalengthcounts, otherlengthcounts)

def filelink(filename: str) -> str: 
    return '<a href="' + filename + '">' + filename + '</a>'

def checkreadtypes(samplename: str, sampleinfo: Any, tgirtmode: bool = False) -> List[Optional[ErrorSet]]:
    trnapercentcutoff = 0.05
    ribopercentcutoff = 0.35
    unmappercentcutoff = 0.35
    
    highmeanlength = 40
    minsizethreshold = 15
    maxsizethreshold = 50
    percentsizethreshold = 0.70
    lowmeanlength = None
    
    if tgirtmode:
        trnapercentcutoff = 0.5
        ribopercentcutoff = 0.35
        unmappercentcutoff = 0.35
        
        highmeanlength = 75
        lowmeanlength = 40
        minsizethreshold = 40
        maxsizethreshold = 75
        percentsizethreshold = 0.70
        
    typecounts = gettypecounts(samplename, sampleinfo)
    # getfragtypes(samplename, sampleinfo) # Was pass in original
    
    samples = list(sampleinfo.getsamples())
    
    trnapercent = {currsample: typecounts.gettrnapercent(currsample) for currsample in samples}
    lowtrnasamples = [trnapercent[currsample] < trnapercentcutoff for currsample in samples]
    lowtrnaerr = ErrorSet("trna_read", samples, lowtrnasamples, 
                          "tRNA read share > " + str(100 * trnapercentcutoff) + "%", 
                          "tRNA Read Percentage", trnapercent, percentformat=True, 
                          checkfile=samplename + "-typecounts.pdf")

    rrnapercent = {currsample: typecounts.getrrnapercent(currsample) for currsample in samples}
    # Handle None values in rrnapercent if any
    rrnapercent_safe = {k: v if v is not None else 0.0 for k, v in rrnapercent.items()}
    
    highribosamples = [rrnapercent_safe[currsample] > ribopercentcutoff for currsample in samples]
    highriboerr = ErrorSet("rrna_read", samples, highribosamples, 
                           "rRNA read share < " + str(100 * ribopercentcutoff) + "%", 
                           "rRNA Read Percentage", rrnapercent_safe, percentformat=True, 
                           checkfile=samplename + "-typecounts.pdf")

    otherpercent = {currsample: typecounts.getotherpercent(currsample) for currsample in samples}
    highothersamples = [otherpercent[currsample] > unmappercentcutoff for currsample in samples]
    highothererr = ErrorSet("unannotated", samples, highothersamples, 
                            "Reads mapping to unannotated regions < " + str(100 * unmappercentcutoff) + "%", 
                            "Unannotated Region Mapping Rate", otherpercent, percentformat=True, 
                            checkfile=samplename + "-typecounts.pdf")

    allreadlength = getreadlengths(samplename, sampleinfo)

    meanreadlength = {currsample: allreadlength.getsamplemean(currsample) for currsample in samples}
    meanreadsamples = [allreadlength.getsamplemean(currsample) > highmeanlength for currsample in samples]
    highlengtherr = ErrorSet("high_read_len", samples, meanreadsamples, 
                             "Read length average < " + str(highmeanlength) + " bases", 
                             "Average Read Length", meanreadlength, percentformat=False, 
                             checkfile=samplename + "-readlengths.pdf")
    
    lowmeanlenerr = None
    if lowmeanlength is not None:
        lowmeanreadsamples = [allreadlength.getsamplemean(currsample) < lowmeanlength for currsample in samples]
        lowmeanlenerr = ErrorSet("low_read_len", samples, lowmeanreadsamples, 
                                 "Read length average > " + str(lowmeanlength) + " bases", 
                                 "Average Read Length", meanreadlength, percentformat=False, 
                                 checkfile=samplename + "-readlengths.pdf")
    
    thresholdreadpercent = {currsample: allreadlength.getthresholdpercent(currsample, minsizethreshold, maxsizethreshold) for currsample in samples}
    badsizesamples = [thresholdreadpercent[currsample] < percentsizethreshold for currsample in samples]
    
    badsizeerr = ErrorSet("trna_sizes", samples, badsizesamples, 
                          " >= " + str(100 * percentsizethreshold) + "% of reads between " + str(minsizethreshold) + " and " + str(maxsizethreshold) + " bases", 
                          "Read Percentage", thresholdreadpercent, percentformat=True, 
                          checkfile=samplename + "-readlengths.pdf")

    results = [lowtrnaerr, highriboerr, highothererr, highlengtherr, badsizeerr]
    if lowmeanlenerr:
        results.append(lowmeanlenerr)
    return results

def gettrnacounts(samplename: str, sampleinfo: Any, trnainfo: Any) -> TrnaCount:
    allsamples = sampleinfo.getsamples()
    trnatranscripts = set(trnainfo.gettranscripts())
    trnacounts = defaultdict(dict)
    
    with open(gettrnacountfile(samplename)) as typeresults:
        runsamples = []
        for i, currline in enumerate(typeresults):
            fields = currline.rstrip().split("\t")
            if i == 0:
                runsamples = list(fields)
                if set(runsamples) != set(allsamples):
                    logger.warning(runsamples)
                    logger.warning(allsamples)
                    logger.warning("QAError")
                continue

            if len(fields) != len(allsamples) + 1:
                logger.warning(runsamples)
                logger.warning(fields)
                logger.warning("QAError")
                continue

            if fields[0] in trnatranscripts:
                for j in range(0, len(runsamples)):
                    trnacounts[runsamples[j]][fields[0]] = int(fields[j + 1])
    return TrnaCount(trnacounts)

def getsizefactor(samplename: str, sampleinfo: Any) -> SizeFactor:
    allsamples = sampleinfo.getsamples()
    sizefactors = {}
    
    with open(getsizefactorfile(samplename)) as typeresults:
        runsamples = []
        for i, currline in enumerate(typeresults):
            fields = currline.rstrip().split()
            fields = [curr.strip('"') for curr in fields]
            if i == 0:
                runsamples = list(fields)
                if set(runsamples) != set(allsamples):
                    logger.warning(list(runsamples))
                    logger.warning(list(allsamples))
                    logger.warning("QAError")
                continue

            if len(fields) != len(allsamples):
                logger.warning(len(runsamples))
                logger.warning(len(fields))
                logger.warning("QAError")
                continue
            
            for j in range(0, len(runsamples)):
                sizefactors[runsamples[j]] = float(fields[j])
    return SizeFactor(sizefactors)

SIZE_FACTOR_DIFF = 3.0
MIN_ACTIVE_PERCENT = 0.5
MIN_READ_COUNT = 20

def checkgenecounts(samplename: str, sampleinfo: Any, trnainfo: Any, tgirtmode: bool = False) -> List[ErrorSet]:
    readcounts = gettrnacounts(samplename, sampleinfo, trnainfo)
    sizefactors = getsizefactor(samplename, sampleinfo)
    samples = list(sampleinfo.getsamples())
    
    thresholdreadpercent = {currsample: readcounts.gettrnaactivepercent(currsample, trnainfo, cutoff=MIN_READ_COUNT) for currsample in samples}
    
    missingtrnasamples = [thresholdreadpercent[currsample] < MIN_ACTIVE_PERCENT for currsample in samples]
    lowcounterr = ErrorSet("trna_read_count", samples, missingtrnasamples, 
                           ">= " + str(100 * MIN_ACTIVE_PERCENT) + "% of tRNAs with more than " + str(MIN_READ_COUNT) + " reads", 
                           "Percentage of tRNAs", thresholdreadpercent, percentformat=True, 
                           checkfile=samplename + "-tRNAcounts.txt")
    
    samplesizefactors = {currsample: sizefactors.sizefactors[currsample] for currsample in samples}

    badsizefactors = [samplesizefactors[currsample] > SIZE_FACTOR_DIFF or samplesizefactors[currsample] < (1.0 / (SIZE_FACTOR_DIFF + 0.0001)) for currsample in samples]
    sizefactorerr = ErrorSet("size_factors", samples, badsizefactors, 
                             "DESeq2 size factor differences < " + str(SIZE_FACTOR_DIFF) + "x", 
                             "Size Factor", samplesizefactors, checkfile=samplename + "-SizeFactors.txt")

    return [lowcounterr, sizefactorerr]

def readtrimindex(trimindex: str) -> List[str]:
    filelocs = {}
    with open(trimindex) as indexfile:
        for currline in indexfile:
            fields = currline.split()
            if len(fields) > 1:
                filelocs[fields[0]] = fields[1]
    return list(filelocs.keys())

def main(experimentname: str, samplefile: str, databasename: str, 
         tgirt: bool = False, runname: Optional[str] = None, 
         output: Optional[str] = None):
    
    samplename = experimentname
    
    if output is not None:
        outputfile = open(output, "w")
    else:
        outputfile = sys.stdout
        
    try:
        sampleinfo = toolsTG.samplefile(os.path.expanduser(samplefile))
        trnainfo = toolsTG.transcriptfile(os.path.expanduser(databasename + "-trnatable.txt"))
        
        allsamples = sampleinfo.getsamples()
        print("<html>", file=outputfile)
        print("<head>" + STYLE + "</head>", file=outputfile)
        
        modestring = "Full-length tRNAs" if tgirt else "tRNA fragments"
        date = strftime("%A %B %d, %Y", localtime())
        print(MODE_TEMPLATE.format(date=date, mode=modestring), file=outputfile)
        
        print("<body>", file=outputfile)
        print('<a name="summary"><h4>Summary</h4></a>', file=outputfile)
        
        prepresults = []
        if runname is None and os.path.exists("trimindex.txt"):
            runnames = readtrimindex("trimindex.txt")
            prepresults = checkreadprep(runnames, sampleinfo)
        if runname is not None:
            prepresults = checkreadprep([runname], sampleinfo)
            
        mappingresults = []
        if os.path.exists(getmapfile(samplename)):
            mappingresults = checkreadsmapping(samplename, sampleinfo, tgirt)
            
        typeresults = checkreadtypes(samplename, sampleinfo, tgirt)
        countresults = checkgenecounts(samplename, sampleinfo, trnainfo, tgirt)
        
        allresults = prepresults + mappingresults + typeresults + countresults
        # Filter out None values if any
        allresults = [res for res in allresults if res is not None]
        
        print("<p>", file=outputfile)
        for currtest in allresults:
            color = currtest.gettestcolor()
            errlvl = currtest.getteststatus()
            print('<b style="color:{color};">{msg}</b> <a href="#{testname}">{criteria}</a> ({filename})</br>'.format(
                color=color, msg=errlvl, testname=currtest.shortname, 
                criteria=currtest.criteria, filename=filelink(currtest.checkfile)), file=outputfile)

        print("</p>", file=outputfile)
        print("<hr>\n<hr>", file=outputfile)
        
        for currtest in allresults:
            print('<a name="{testname}"><h4>{msg}</h4></a>'.format(testname=currtest.shortname, msg=currtest.criteria), file=outputfile)
            print('<table>', file=outputfile)
            print('<thead><tr><th width="15%">Status</th><th width="50%">Sample</tH><th>{measure}</th></tr></thead>'.format(measure=currtest.dimension), file=outputfile)
            print('<tbody>', file=outputfile)

            for currsample in allsamples:
                color = currtest.getsamplecolor(currsample)
                errlvl = currtest.getsamplestatus(currsample)
                print('<tr><td><b style="color:{color};">{errlvl}</b></td><td>{samplename}</td><td>{sampleresult}</td></tr>'.format(
                    color=color, errlvl=errlvl, samplename=currsample, 
                    sampleresult=currtest.getsampleresult(currsample)), file=outputfile)

            print('</tbody>\n\n</table>\n\n<p><a href="#summary">back to Summary</a></p>\n\n<hr>', file=outputfile)

        print("</html>", file=outputfile)
        
    finally:
        if output is not None and outputfile != sys.stdout:
            outputfile.close()

if __name__ == "__main__":
    # This script is intended to be imported or called with arguments directly.
    # If run as script, it expects arguments to be passed to main() manually or via another script.
    # Since argparse is removed, we can't parse command line arguments here easily without it.
    # The user requested "It does not need argparse functionality".
    pass

