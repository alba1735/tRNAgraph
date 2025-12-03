#!/usr/bin/env python3

import pysam
import sys
import os.path
from collections import defaultdict
import itertools
from multiprocessing import Process, Queue, Pool
from . import toolsTG

def getdupes(namelist):
    allset = set()
    for currname in namelist:
        if currname in allset:
            yield currname
        else:
            allset.add(currname)

def enddict():
    return defaultdict(int)

class featurecount:
    def __init__(self, samplename, bamfile, trnas = list(), trnaloci = list(), emblgenes = list(), otherfeats = list()):
        self.samplename = samplename
        self.bamfile = bamfile
        self.trnas = trnas
        self.trnaloci = trnaloci
        self.emblgenes = emblgenes
        self.otherfeats = otherfeats
        
        self.counts = defaultdict(int)
        self.trnacounts = defaultdict(int)
        self.antitrnacount = defaultdict(int)
        self.trnawholecounts = defaultdict(int)
        self.trnafivecounts = defaultdict(int)
        self.trnathreecounts = defaultdict(int)
        self.trnalocuscounts = defaultdict(int)
        self.trnalocustrailercounts = defaultdict(int)
        self.partialtrnalocuscounts = defaultdict(int)
        self.fulltrnalocuscounts  = defaultdict(int)
        self.trnauniquecounts = defaultdict(int)
        self.aminocounts  = defaultdict(int)
        self.anticodoncounts =  defaultdict(int) 
        self.trnaendtypecounts = defaultdict(enddict)
        self.lengthsum = defaultdict(int) 
        self.lengthtotal = defaultdict(int) 
        
        self.gcpercent = defaultdict(int) 
        self.gctotal = defaultdict(int) 
        
        self.genetypes = dict()
        
        self.trnauniquewholecounts = defaultdict(int) 
        self.trnauniquefivecounts = defaultdict(int) 
        self.trnauniquethreecounts = defaultdict(int) 
        
    def setgenetype(self, genename, genetype):
        self.genetypes[genename] = genetype
    def addcount(self, genename):
       self.counts[genename] += 1
    def addantitrnacount(self, genename):
       self.antitrnacount[genename] += 1
       
    def addlocuscount(self, genename):
       self.trnalocuscounts[genename] += 1
    def addpartiallocuscount(self, genename):
       self.partialtrnalocuscounts[genename] += 1
    def addfulllocuscount(self, genename):
       self.fulltrnalocuscounts[genename] += 1
    def addlocustrailercount(self, genename):
       self.trnalocustrailercounts[genename] += 1
    def addtrnacount(self, genename):
       self.trnacounts[genename] += 1
    def adduniquecount(self, genename, fragtype):
       self.trnauniquecounts[genename] += 1
       if fragtype == "Whole":
           self.trnauniquewholecounts[genename] += 1
       elif fragtype == "Fiveprime":
           self.trnauniquefivecounts[genename] += 1
       elif fragtype == "Threeprime":
           self.trnauniquethreecounts[genename] += 1
       else:
           pass
    def addaminocount(self, amino):
       self.aminocounts[amino] += 1
    def addanticodoncount(self, anticodon):
       self.anticodoncounts[anticodon] += 1
    def addfragcount(self, featname, fragtype):    
        if fragtype == "Whole":
            self.trnawholecounts[featname] += 1
        elif fragtype == "Fiveprime":
            self.trnafivecounts[featname] += 1
        elif fragtype == "Threeprime":
            self.trnathreecounts[featname] += 1
    def addendcount(self, featname, endtype):    
        if endtype is not None:
            self.trnaendtypecounts[featname][endtype] += 1
            
    def addreadlength(self, genename, length):
       self.lengthsum[genename] += length
       self.lengthtotal[genename] += 1
       
    def addgc(self, genename, gc, length):
       self.gcpercent[genename] += gc
       self.gctotal[genename] += length
       
       
    def getgenecount(self, genename):
       return self.counts[genename]
    def getantitrnacount(self, genename):
       return self.antitrnacount[genename]
    def getlocuscount(self, genename):
       return self.trnalocuscounts[genename]
    def getpartiallocuscount(self, genename):
       return self.partialtrnalocuscounts[genename]
    def getfulllocuscount(self, genename):
       return self.fulltrnalocuscounts[genename] 
    def getlocustrailercount(self, genename):
       return self.trnalocustrailercounts[genename]
    def gettrnacount(self, genename):
       return self.trnacounts[genename]
       
       
    def getuniquecount(self, genename):
       return self.trnauniquecounts[genename]
    def getfiveuniquecount(self, genename):
       return self.trnauniquefivecounts[genename]
    def getthreeuniquecount(self, genename):
       return self.trnauniquethreecounts[genename]
    def getwholeuniquecount(self, genename):
       return self.trnauniquewholecounts[genename]
    def getotheruniquecount(self, genename):
       return self.trnauniquecounts[genename] - (self.trnauniquewholecounts[genename] + self.trnauniquethreecounts[genename] + self.trnauniquefivecounts[genename])
       
       
    def getaminocount(self, amino):
       return self.aminocounts[amino]
    def getanticodoncount(self, anticodon):
       return self.anticodoncounts[anticodon]

    def getfivecount(self, genename):
       return self.trnafivecounts[genename]
    def getthreecount(self, genename):
       return self.trnathreecounts[genename]
    def getwholecount(self, genename):
       return self.trnawholecounts[genename]
    def getendtypecount(self, genename):
       return self.trnaendtypecounts[genename]
    def getreadavglength(self, genename):
        if self.lengthtotal[genename] == 0:
            return 0
        else:
            return self.lengthsum[genename] / self.lengthtotal[genename]
    def getreadavggc(self, genename):
        if self.gctotal[genename] == 0:
            return 0
        else:
            return self.gcpercent[genename] / self.gctotal[genename]

allfragtypes = set(["Whole","Fiveprime","Threeprime","Other"])
fragnames = {"Whole":"wholecounts","Fiveprime":"fiveprime","Threeprime":"threeprime","Other":"other"}

def getchromdict(features):
    chromdict = defaultdict(list)
    for curr in features:
        chromdict[curr.chrom].append(curr)
    return chromdict

class counttypes:
    def __init__(self, samplename, bamfile, trnas = list(), trnaloci = list(), emblgenes = list(), otherfeats = list()):
        self.samplename = samplename
        self.bamfile = bamfile
        self.trnas = trnas
        self.trnaloci = trnaloci
        self.emblgenes = emblgenes
        self.otherfeats = otherfeats
        self.emblbiotypes = set()
        self.aminos = set()
        self.bedtypes = set()
        self.extraseqtypes = set()
        self.anticodons = set()
        
        self.embltypecounts =defaultdict(int)
        self.bedtypecounts =defaultdict(int)
        self.trnafragtypes =defaultdict(int)
        self.trnafragtypes =defaultdict(int)
        self.totalreads = 0
        self.trnareads = 0
        self.otherreads = 0
        self.readlengths = defaultdict(int)
        self.trnareadlengths = defaultdict(int)
        self.trnacounts = defaultdict(int)
        self.aminocounts = defaultdict(int)
        self.anticodoncounts = defaultdict(int)
        self.indelreads = defaultdict(int)
        self.trnaanticounts = defaultdict(int)
        self.trnaemblcounts = defaultdict(int)
        self.pretrnareadlengths = defaultdict(int)
        self.trnalocuscounts = defaultdict(int)
        self.partiallocuscounts = defaultdict(int)
        self.fulllocuscounts = defaultdict(int)
        self.extraseqcounts = defaultdict(int)
        self.trnaantilocuscounts  = defaultdict(int)
        self.mismatchcounts  = defaultdict(int)
        self.trnamismatchcounts  = defaultdict(int)
        self.aminouniqcounts = dict()
        self.anticodonuniqcounts = dict()
        self.trnatranscriptuniqcounts = dict()
        self.trnatranscriptcounts = dict()
        self.aminotranscriptcounts = dict()
        self.anticodontranscriptcounts = dict()
        for currfrag in allfragtypes:
            self.aminouniqcounts[currfrag] = defaultdict(int)
            self.anticodonuniqcounts[currfrag] = defaultdict(int)
            self.trnatranscriptuniqcounts[currfrag] = defaultdict(int)
            self.trnatranscriptcounts[currfrag] = defaultdict(int)
            self.aminotranscriptcounts[currfrag] = defaultdict(int)
            self.anticodontranscriptcounts[currfrag] = defaultdict(int)
    def addsamplecounts(self):
        self.totalreads += 1
    def addreadlengths(self, length):
        self.readlengths[length] += 1
    def addtrnareadlengths(self, length):
        self.trnareadlengths[length] += 1
    def addpretrnareadlengths(self, length):
        self.pretrnareadlengths[length] += 1
    def addpartiallocuscounts(self, currbed):
        self.partiallocuscounts[currbed] += 1     
    def addtrnasamplecounts(self):
        self.trnareads += 1

    def addtrnaantilocuscounts(self, currbed):
        self.trnaantilocuscounts[currbed] += 1
    def addtrnalocuscounts(self, currbed):
        self.trnalocuscounts[currbed] += 1
    def addfulllocuscounts(self, currbed):
        self.fulllocuscounts[currbed] += 1
    def addaminocounts(self, curramino, fragtype = None, unique = None):
        self.aminos.add(curramino)
        self.aminocounts[curramino] += 1
        if fragtype is not None:
            self.aminotranscriptcounts[fragtype][curramino] += 1
        if unique is None or unique:

            if fragtype is not None:
                self.aminouniqcounts[fragtype][curramino] += 1
    def addanticodoncounts(self, curranticodon, fragtype = None, unique = None):
        self.anticodons.add(curranticodon)
        self.anticodoncounts[curranticodon] += 1
        if fragtype is not None:
            self.anticodontranscriptcounts[fragtype][curranticodon] += 1
        if unique is None or unique:

            if fragtype is not None:
                self.anticodonuniqcounts[fragtype][curranticodon] += 1
                
    def addtrnacounts(self, currbed, currtrna, fragtype = None, unique = None):
        self.trnacounts[currbed] += 1
        if fragtype is not None:
            self.trnatranscriptcounts[fragtype][currtrna] += 1
        if unique is None or unique:

            if fragtype is not None:
                self.trnatranscriptuniqcounts[fragtype][currtrna] += 1
                
                
                
    def addindelreads(self, curramino):
        self.indelreads[curramino] += 1
    def addmismatchcounts(self, mismatchcounts):
        self.mismatchcounts[mismatchcounts] += 1
    def addtrnamismatchcounts(self, mismatchcounts):
        self.trnamismatchcounts[mismatchcounts] += 1
    def addotherreads(self):
        self.otherreads += 1
    def addtrnaantisense(self, currbed):
        self.trnaanticounts[currbed] += 1
    def addemblcounts(self, currtype):
        self.emblbiotypes.add(currtype)
        self.embltypecounts[currtype] += 1     
    def addbedcounts(self, genetype):
        self.bedtypes.add(genetype)
        self.bedtypecounts[genetype] += 1
    def addextracounts(self, genetype):
        self.extraseqtypes.add(genetype)
        self.extraseqcounts[genetype] += 1
    def anticodonallfrags(self, curranticodon, uniq = False):
        if uniq:
            return sum(self.anticodonuniqcounts[fragtype][curranticodon] for fragtype in allfragtypes)
        else:
            return sum(self.anticodontranscriptcounts[fragtype][curranticodon] for fragtype in allfragtypes)
    def trnaallfrags(self, currtrna, uniq = False):
        if uniq:
            return sum(self.trnatranscriptuniqcounts[fragtype][currtrna] for fragtype in allfragtypes)
        else:
            return sum(self.trnatranscriptcounts[fragtype][currtrna] for fragtype in allfragtypes)
    def aminoallfrags(self, curramino, uniq = False):
        if uniq:
            return sum(self.aminouniqcounts[fragtype][curramino] for fragtype in allfragtypes)
        else:
            return sum(self.aminotranscriptcounts[fragtype][curramino] for fragtype in allfragtypes)

def getbamcounts(bamfile, samplename,trnainfo, trnaloci, trnalist,featurelist = dict(),otherseqdict = dict(), embllist = list(), bedfiles = list(),nomultimap = False, allowindels = True, maxmismatches = None):
    samplecounts = featurecount(samplename, bamfile, trnas = trnalist, trnaloci = trnaloci, emblgenes = embllist, otherfeats = featurelist)
    fullpretrnathreshold = 2
    minpretrnaextend = 5
    #minimum mapq
    #nomultimap = False
    minmapq = 0
    if nomultimap:
        minmapq = 2
    #minimum number of reads for a feature to be reported
    minreads = 5
    #print >>sys.stderr, embllist
    
    genetypes = dict()
    currbam = bamfile
    
    for currfile in bedfiles:
        bedfeatures = list(toolsTG.readfeatures(currfile, removepseudo = False))
        for curr in bedfeatures:
            genetypes[curr.name] = os.path.basename(currfile)
            
        featurelist[currfile] = bedfeatures
    
    #print >>sys.stderr, currsample
    #doing this thing here why I only index the bamfile if the if the index file isn't there or is older than the map file
    try:
        if not os.path.isfile(currbam+".bai") or os.path.getmtime(currbam+".bai") < os.path.getmtime(currbam):
            pysam.index(""+currbam)
        bamfile = pysam.Samfile(""+currbam, "rb" )  
    except IOError as xxx_todo_changeme1:
        ( strerror) = xxx_todo_changeme1
        print(strerror, file=sys.stderr)
        sys.exit(1)
    except pysam.utils.SamtoolsError:
        print("Can not index "+currbam, file=sys.stderr)
        print("Exiting...", file=sys.stderr)
        sys.exit(1)
        
    
    for currfile in featurelist.keys():
        for currfeat in featurelist[currfile]:
            #try catch is to account for weird chromosomes and the like that aren't in the genome
            #means that if I can't find a feature, I record no counts for it rather than bailing
            try:
                for currread in toolsTG.getbamrange(bamfile, currfeat, singleonly = nomultimap, maxmismatches = maxmismatches,allowindels = allowindels):
                    if currfeat.coverage(currread) > 10:
                        samplecounts.addcount(currfeat.name)
                        samplecounts.addreadlength(currfeat.name, currread.length())
                        #samplecounts.addgc(currfeat.name, currread.getgc(), currread.length())
                        samplecounts.setgenetype(currfeat.name,os.path.basename(currfile))
            except ValueError:
                pass

    #extra sequences built during database creation (experimental)
    for currtype in otherseqdict.keys():
        for currfeat in otherseqdict[currtype]:
            for currread in toolsTG.getbamrange(bamfile, currfeat, singleonly = nomultimap, maxmismatches = maxmismatches,allowindels = allowindels):
                #print >>sys.stderr, currfeat.name
                samplecounts.addcount(currfeat.name)
                samplecounts.addreadlength(currfeat.name, currread.length())
                #samplecounts.addgc(currfeat.name, currread.getgc(), currread.length())
                samplecounts.setgenetype(currfeat.name,currtype)
    for genename, featset in itertools.groupby(embllist,lambda x: x.data["genename"]):
        #print >>sys.stderr, "**"
        #pass 
        try:
            allreads = set()
            for currfeat in list(featset):
                
                for currread in toolsTG.getbamrangeshort(bamfile, currfeat, singleonly = nomultimap, maxmismatches = maxmismatches,allowindels = allowindels, skiptags = True):
                    #print >>sys.stderr, "**"+currread.name 
                    #continue
                    
                    if currfeat.coverage(currread) > 10:
                        
                        samplecounts.addcount(genename)
                        samplecounts.addreadlength(currfeat.name, currread.length())
                        #samplecounts.addgc(currfeat.name, currread.getgc(), currread.length())
                        #print >>sys.stderr, "**"+currread.name
                        samplecounts.setgenetype(genename,currfeat.data["biotype"])
                        #print >>sys.stderr, currfeat.bedstring()
        except ValueError:
            pass
    for currfeat in trnaloci:
        #print >>sys.stderr,  currfeat.bedstring()
        #print >>sys.stderr,  currfeat.getdownstream(30).bedstring()
        for currread in toolsTG.getbamrangeshort(bamfile, currfeat.addmargin(30), singleonly = nomultimap, maxmismatches = maxmismatches,allowindels = allowindels, skiptags = True):
            #gotta be more than 5 bases off one end to be a true pre-tRNA
            #might want to shove these to the real tRNA at some point, but they are for now just ignored
            
            if currfeat.coverage(currread) > 10 and (currread.start + minpretrnaextend <= currfeat.start or currread.end - minpretrnaextend >= currfeat.end):
                samplecounts.addlocuscount(currfeat.name)
                samplecounts.addreadlength(currfeat.name, currread.length())
                #samplecounts.addgc(currfeat.name, currread.getgc(), currread.length())
                if currread.start + fullpretrnathreshold <  currfeat.start and currread.end - fullpretrnathreshold + 3 >  currfeat.end:
                    samplecounts.addfulllocuscount(currfeat.name)
                else:
                    #partialtrnalocuscounts[currsample][currfeat.name] += 1
                    samplecounts.addpartiallocuscount(currfeat.name)
            elif currfeat.getdownstream(30).coverage(currread) > 10:  #need the elif otherwise fragments that include trailer get in there
                samplecounts.addlocuscount(currfeat.name)
                samplecounts.addlocustrailercount(currfeat.name)
            else:
                #print >>sys.stderr,  currfeat.getdownstream(30).coverage(currread)
                pass
    #print >>sys.stderr, samplename+" threadA "+str(time.time())        
    for currfeat in trnalist:
        #print >>sys.stderr, samplename+":"+currfeat.name
        
        featreads = 0
        for currread in toolsTG.getbam(bamfile, currfeat, singleonly = nomultimap, allowindels = allowindels):
            #samplecounts.addgc(currfeat.name, currread.getgc(), currread.length())
            if maxmismatches is not None and currread.getmismatches() > maxmismatches:
                continue
            samplecounts.addreadlength(currfeat.name, currread.length())

            featreads += 1
            if not currfeat.strand == currread.strand:
                samplecounts.addantitrnacount(currfeat.name)
                continue
            if not currfeat.coverage(currread) > 10:
                continue

            curramino = trnainfo.getamino(currfeat.name)
            curranticodon = trnainfo.getanticodon(currfeat.name)
            #samplecounts.addfragcount(currfeat.name, fragtype)
            samplecounts.addtrnacount(currfeat.name)
                
            fragtype = toolsTG.getfragtype(currfeat, currread)
            samplecounts.addfragcount(currfeat.name, fragtype)
            endtype = toolsTG.getendtype(currfeat, currread)
            #print >>sys.stderr, endtype
            samplecounts.addendcount(currfeat.name, endtype)
            if currread.isuniquetrnamapping():
                samplecounts.adduniquecount(currfeat.name, fragtype)
            if currread.isuniqueaminomapping():
                pass
            elif currread.isuniqueacmapping():
                samplecounts.addaminocount(curramino)
            else:
                samplecounts.addaminocount(curramino)
                samplecounts.addanticodoncount(curramino)
        #print >>sys.stderr, str(featreads)+"/"+str(samplecounts.gettrnacount(currfeat.name))
    #print >>sys.stderr, samplename+" thread "+str(time.time())
            
    return samplecounts

def counttypereads(bamfile, samplename,trnainfo, trnaloci, trnalist,maturenames,featurelist = dict(), otherseqlist = list(), embllist = list(),nomultimap = False, allowindels = True, maxmismatches = None, bamnofeature = False, countfrags = False):
    
    bedlist = list(featurelist.keys())
    readtypecounts = counttypes(samplename, bamfile, trnas = trnalist, trnaloci = trnaloci, emblgenes = embllist, otherfeats = bedlist)
    mitochrom = None
    fullpretrnathreshold = 2
    minpretrnaextend = 5
    ncrnaorder = defaultdict(int)
    currbam = bamfile
    dumpotherreads = True

    for i, curr in enumerate(reversed(list(["snoRNA","miRNA", "rRNA","snRNA","misc_RNA","lincRNA", "protein_coding"]))):
        ncrnaorder[curr] = i + 1

        
    try:
        #print >>sys.stderr, currbam
        
        
        if not os.path.isfile(currbam+".bai") or  os.path.getmtime(currbam+".bai") < os.path.getmtime(currbam):
            pysam.index(""+currbam)
        bamfile = pysam.Samfile(""+currbam, "rb" )
        if bamnofeature:
            outname = os.path.splitext(currbam)[0]+"_nofeat.bam"
            outbamnofeature =  pysam.Samfile( outname, "wb", template =  bamfile)
    except IOError as xxx_todo_changeme1:
        ( strerror) = xxx_todo_changeme1
        print(strerror, file=sys.stderr)
        sys.exit()
    #continue #point0
    #print >>sys.stderr, "**||"+currbam

    for i, currread in enumerate(toolsTG.getbam(bamfile, primaryonly = True)):

        isindel = False
        hasmiamatch  = False
        readlength = currread.getlength()
        gotread = False
        #continue #point1
        readtypecounts.addsamplecounts()
        readtypecounts.addmismatchcounts(currread.getmismatches())
        if currread.hasindel():
            readtypecounts.addindelreads(readlength)
            isindel = True
            #continue

        else:
            pass
            #continue
        readtypecounts.addreadlengths(readlength)
        #readlengths[currsample][readlength] += 1
        #continue #point2
        #if currread.name == "NB501427:473:H3YJ2BGXG:3:11512:1635:4586":
        #    print >>sys.stderr, "**foundread"
        for currbed in trnaloci:
            for currfeat in trnaloci[currbed].getbin(currread):
                expandfeat = currfeat.addmargin(30)
                #if currread.name == "NB501427:473:H3YJ2BGXG:3:11512:1635:4586"  and currfeat.name == "tRNA-Val-AAC-5-1":
                #    print >>sys.stderr, "**trnaread:" +str(currfeat.coverage(currread))
                #    print >>sys.stderr, "**trnaread:" +str(currfeat.coverage(currread))
                
                cov = currfeat.coverage(currread)
                if cov > 10: # and (currread.start + minpretrnaextend <= currfeat.start or currread.end - minpretrnaextend >= currfeat.end):
                    #if currread.name == "NB501427:473:H3YJ2BGXG:3:11512:1635:4586":
                    #    print >>sys.stderr, "**foundread:" +str(currfeat.name)
                    if currfeat.strand != currread.strand:
                        readtypecounts.addtrnaantilocuscounts(currbed)
                        break
                    if (currread.start + minpretrnaextend <= currfeat.start or currread.end - minpretrnaextend >= currfeat.end):
                        pass
                    readtypecounts.addpretrnareadlengths(readlength)
                    readtypecounts.addtrnalocuscounts(currbed)
                    if currread.start + fullpretrnathreshold <  currfeat.start and currread.end - fullpretrnathreshold + 3 >  currfeat.end:
                        #fulltrnalocuscounts[currsample][currbed] += 1
                        readtypecounts.addfulllocuscounts(currbed)
                    else:# currread.start + fullpretrnathreshold <  currfeat.start or currread.end - fullpretrnathreshold +3 >  currfeat.end:
                        #partialtrnalocuscounts[currsample][currbed] += 1
                        readtypecounts.addpartiallocuscounts(currbed)
                        #print >>sys.stderr, "***"
                    gotread = True
                    break
                
                if currfeat.getdownstream(30).coverage(currread) > 10:
                    readtypecounts.addtrnaantisense(currbed)
                    #readtypecounts.addpretrnareadlengths(readlength)
                    #print >>sys.stderr, currfeat.bedstring()
                    gotread = True
                    break
                elif expandfeat.antisense().coverage(currread) > 5:
                    #trnaantisense[currsample][currbed] += 1
                    readtypecounts.addtrnaantisense(currbed)
                    gotread = True
                    break
        #if currread.name == "NB501427:473:H3YJ2BGXG:3:11512:1635:4586":
        #    print >>sys.stderr, "**foundread2"            
        if gotread: 
            continue
        #continue #point3

        for currbed in trnalist:
            if currread.chrom in maturenames[currbed]:
                currfeat = maturenames[currbed][currread.chrom]
                if currread.strand == "+":
                    

                    fragtype = None
                    #aminos.add(trnainfo.getamino(currfeat.name))
                    fragtype = toolsTG.getfragtype(currfeat, currread)
                    if fragtype == "Whole":
                        
                        #trnawholecounts[currsample][currbed] += 1
                        pass
                    elif fragtype == "Fiveprime":
                        #trnafivecounts[currsample][currbed] += 1
                        pass
                    elif fragtype == "Threeprime":
                        pass
                        #trnathreecounts[currsample][currbed] += 1
                    elif fragtype == "Trailer":
                        #trnatrailercounts[currsample][currbed] += 1
                        pass
                    else:
                        fragtype = "Other"
                    readtypecounts.addtrnareadlengths(readlength)
                    readtypecounts.addtrnasamplecounts()
                    
                    readtypecounts.addaminocounts(trnainfo.getamino(currfeat.name), fragtype = fragtype, unique = currread.isuniqueaminomapping())
                    readtypecounts.addanticodoncounts(trnainfo.getanticodon(currfeat.name), fragtype = fragtype, unique = currread.isuniqueaminomapping())
                    readtypecounts.addtrnacounts(currbed, currfeat.name, fragtype = fragtype, unique = currread.isuniquetrnamapping())

                    readtypecounts.addtrnamismatchcounts(currread.getmismatches())
                    
                    gotread = True
                    break
                        #print >>sys.stderr, str(currread.start - currfeat.start)+"-"+str(currread.end - currfeat.start)  
                        #print >>sys.stderr, str(currfeat.start - currfeat.start)+"-"+str(currfeat.end - currfeat.start)
                        #print >>sys.stderr, "****"
                elif currfeat.antisense().coverage(currread) > 10:
                    readtypecounts.addtrnaantisense(currbed)
                    #trnaantisense[currsample][currbed] += 1
                    gotread = True
                    break
        if gotread: 
            continue
        #continue #point4
        if embllist is not None:
            currtype = None
            for currfeat in embllist.getbin(currread):
                if currfeat.coverage(currread) > 10: 
                    if currfeat.data["biotype"] == "processed_transcript":
                        #print >>sys.stderr, currfeat.bedstring()
                        
                        pass

                    if currtype is None or ncrnaorder[currfeat.data["biotype"]] > ncrnaorder[currtype]:
                        currtype= currfeat.data["biotype"]
                        #if mitochrom == currread.chrom:
                            #currtype = "mito"+currtype
                    
                    
                    
            if currtype is not None:
                readtypecounts.addemblcounts(currtype)
                #emblcounts[currsample][currtype] += 1
                #emblbiotypes.add(currtype)
                gotread = True
                    #print >>sys.stderr, currbam +":"+ currbed
        if gotread: 
            continue
        #continue #point5
        #print >>sys.stderr, "**||"
        
        for currbed in bedlist:

            #if currread.name == "SRR10038183.1660151":
            #    print >>sys.stderr, "||"+currbed
            #    
            #    print >>sys.stderr, list(featurelist[currbed].getbin(currread))
            #    print >>sys.stderr, list(featurelist[currbed].getfeatbin("12-qE-23911.2"))
            #    print >>sys.stderr, list(featurelist[currbed].getbinnums(currread))
            
            for currfeat in featurelist[currbed].getbin(currread):
                
                if currfeat.coverage(currread) > 10:
                    #print >>sys.stderr, currbam +":"+ currbed
                    readtypecounts.addbedcounts(currbed)
                    #counts[currsample][currbed] += 1
                    gotread = True
                    break
                    #print >>sys.stderr, currbam +":"+ currbed
        if gotread: 
            continue
        for currbed in otherseqlist:
            #print >>sys.stderr, len(list(otherseqlist[currbed].getbin(currread)))
            #print >>sys.stderr, "**"+ otherseqlist[currbed].keys()[0]
            #print >>sys.stderr, "||"+ currfeat.chrom
            for currfeat in otherseqlist[currbed][currread.chrom]:
                #print >>sys.stderr, "**"+currbed
                if currfeat.coverage(currread) > 10:
                    
                    readtypecounts.addextracounts(currbed)
                    gotread = True
                    break 
                    #print >>sys.stderr, currbam +":"+ currbed
        #print >>sys.stderr, "**||"
        if gotread: 
            continue
        readtypecounts.addotherreads()
        if not gotread and embllist is not None and mitochrom == currread.chrom:
            currtype = "Mitochondrial_other"
            readtypecounts.addemblcounts(currtype)
            #emblbiotypes.add(currtype)
        if not gotread and bamnofeature:
            outbamnofeature.write(currread.bamline)
    
    return readtypecounts

def counttypereadsqueue(countqueue,currsample, *args, **kwargs):
    countqueue.put([currsample,counttypereads(*args, **kwargs)])

def counttypereadspool(args):
    return counttypereads(*args[0], **args[1])

def printcountfile(countfile, samples,  samplecounts, trnalist, trnaloci, featurelist, embllist, otherseqdict = dict(),minreads = 5, includebase = False):
    print("\t".join(samples), file=countfile)
    trnanames = set()
    for currfeat in trnalist:
        #print >>sys.stderr, samplecounts
        if max(itertools.chain((samplecounts[currsample].gettrnacount(currfeat.name) for currsample in samples), [0])) < minreads:
            continue
        if includebase:
            print(currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].gettrnacount(currfeat.name)) for currsample in samples), file=countfile)
            print(currfeat.name+"_antisense\t"+"\t".join(str(samplecounts[currsample].getantitrnacount(currfeat.name)) for currsample in samples), file=countfile)

        else:
            print(currfeat.name+"_wholecounts\t"+"\t".join(str(samplecounts[currsample].getwholecount(currfeat.name)) for currsample in samples), file=countfile)
            print(currfeat.name+"_fiveprime\t"+"\t".join(str(samplecounts[currsample].getfivecount(currfeat.name) ) for currsample in samples), file=countfile)
            print(currfeat.name+"_threeprime\t"+"\t".join(str(samplecounts[currsample].getthreecount(currfeat.name)) for currsample in samples), file=countfile)
            print(currfeat.name+"_other\t"+"\t".join(str(samplecounts[currsample].gettrnacount(currfeat.name) - (samplecounts[currsample].getwholecount(currfeat.name) + samplecounts[currsample].getfivecount(currfeat.name) + samplecounts[currsample].getthreecount(currfeat.name))) for currsample in samples), file=countfile)
            

            print(currfeat.name+"_antisense\t"+"\t".join(str(samplecounts[currsample].getantitrnacount(currfeat.name)) for currsample in samples), file=countfile)

    for currfeat in trnaloci:
        if max(itertools.chain((samplecounts[currsample].getlocuscount(currfeat.name) for currsample in samples),[0])) < minreads:
            continue
        if includebase:
            print(currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].getlocuscount(currfeat.name)) for currsample in samples), file=countfile)
        else:
            print(currfeat.name+"_wholeprecounts\t"+"\t".join(str(samplecounts[currsample].getfulllocuscount(currfeat.name) ) for currsample in samples), file=countfile)
            print(currfeat.name+"_partialprecounts\t"+"\t".join(str(samplecounts[currsample].getpartiallocuscount(currfeat.name) ) for currsample in samples), file=countfile)
            print(currfeat.name+"_trailercounts\t"+"\t".join(str(samplecounts[currsample].getlocustrailercount(currfeat.name)) for currsample in samples), file=countfile)        
    for currbed in featurelist.keys():
        for currfeat in featurelist[currbed]:
            if currfeat.name in trnanames:
                continue
            trnanames.add(currfeat.name)
            if max(samplecounts[currsample].getgenecount(currfeat.name) for currsample in samples) > minreads:
                print(currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].getgenecount(currfeat.name)) for currsample in samples), file=countfile)
    for currtype in otherseqdict.keys():        
        for currfeat in otherseqdict[currtype] :
        
            trnanames.add(currfeat.name)
            if max(samplecounts[currsample].getgenecount(currfeat.name) for currsample in samples) > minreads:
                print(currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].getgenecount(currfeat.name)) for currsample in samples), file=countfile)
    for currfeat in embllist:
        
        genename = currfeat.data['genename']
        if genename in trnanames:
            continue
        trnanames.add(genename)
        
        if genename is None:
            print(currfeat.name, file=sys.stderr)
            sys.exit(1)
        #print >>sys.stderr, list(samplecounts[currsample].getgenecount(currfeat.name) for currsample in samples)
        if max(samplecounts[currsample].getgenecount(genename) for currsample in samples) > minreads:
            print(genename+"\t"+"\t".join(str(samplecounts[currsample].getgenecount(genename)) for currsample in samples), file=countfile)

def printtrnacountfile(trnacountfilename,samples,  samplecounts, trnalist , includebase = False, minreads = 5):
    trnacountfile = open(trnacountfilename, "w")
    print("\t".join(samples), file=trnacountfile)
    for currfeat in trnalist:
        #print >>sys.stderr, samplecounts
        if max(itertools.chain((samplecounts[currsample].gettrnacount(currfeat.name) for currsample in samples), [0])) < minreads:
            continue
        if includebase:
            print(currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].gettrnacount(currfeat.name)) for currsample in samples), file=trnacountfile)

        else:
            print(currfeat.name+"_wholecounts\t"+"\t".join(str(samplecounts[currsample].getwholecount(currfeat.name)) for currsample in samples), file=trnacountfile)
            print(currfeat.name+"_fiveprime\t"+"\t".join(str(samplecounts[currsample].getfivecount(currfeat.name) ) for currsample in samples), file=trnacountfile)
            print(currfeat.name+"_threeprime\t"+"\t".join(str(samplecounts[currsample].getthreecount(currfeat.name)) for currsample in samples), file=trnacountfile)
            print(currfeat.name+"_other\t"+"\t".join(str(samplecounts[currsample].gettrnacount(currfeat.name) - (samplecounts[currsample].getwholecount(currfeat.name) + samplecounts[currsample].getfivecount(currfeat.name) + samplecounts[currsample].getthreecount(currfeat.name))) for currsample in samples), file=trnacountfile)
            


def averagesamples(allcounts, genename,samples):
    
    return str(sum(allcounts[currsample].lengthsum[genename] for currsample in samples)/(.01+1.*sum(allcounts[currsample].lengthtotal[genename] for currsample in samples)))
    
def gcsamples(allcounts, genename,samples):
        return str(sum(allcounts[currsample].gcpercent[genename] for currsample in samples)/(.01+1.*sum(allcounts[currsample].gctotal[genename] for currsample in samples)))

def printgenetypes(genetypeout,samples, allcounts,trnalist, trnaloci, featurelist, embllist , otherseqdict = dict(),minreads = 5):
    trnanames = set()
    genetypes = dict()
    for currsample in samples:
        genetypes.update(allcounts[currsample].genetypes)
    for currbed in featurelist.keys():
        for currfeat in featurelist[currbed] :
            if currfeat.name in trnanames:
                continue
            trnanames.add(currfeat.name)
            if max(allcounts[currsample].counts[currfeat.name] for currsample in samples) > minreads:
                print(currfeat.name+"\t"+genetypes[currfeat.name]   +"\t"+currfeat.chrom+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
    
        
    for currfeat in trnaloci:
        print(currfeat.name+"_wholeprecounts"+"\t"+"trna_wholeprecounts" +"\t"+currfeat.chrom +"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
        print(currfeat.name+"_partialprecounts"+"\t"+"trna_partialprecounts"+"\t"+currfeat.chrom+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
        print(currfeat.name+"_trailercounts"+"\t"+"trna_trailercounts"+"\t"+currfeat.chrom+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
        print(currfeat.name+""+"\t"+"tRNA_locus"+"\t"+currfeat.chrom+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
    for currfeat in trnalist:
        print(currfeat.name+"_wholecounts"+"\t"+"trna_wholecounts"+"\t"+"tRNA"+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
        print(currfeat.name+"_fiveprime"+"\t"+"trna_fiveprime"+"\t"+"tRNA"+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
        print(currfeat.name+"_threeprime"+"\t"+"trna_threeprime"+"\t"+"tRNA"+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
        print(currfeat.name+"_other"+"\t"+"trna_other"+"\t"+"tRNA"+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
        print(currfeat.name+"_antisense"+"\t"+"trna_antisense"+"\t"+"tRNA"+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
        print(currfeat.name+""+"\t"+"tRNA"+"\t"+"tRNA"+"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)
    
    for currfeat in embllist:
        genename = currfeat.data['genename']
        if genename in trnanames:
            continue
        trnanames.add(genename)
        
        if genename is None:
            #print >>sys.stderr, currfeat.name
            continue
        if max(allcounts[currsample].counts[genename] for currsample in samples) > minreads:
            print(genename+"\t"+genetypes[genename]    +"\t"+currfeat.chrom +"\t"+averagesamples(allcounts, currfeat.name, samples), file=genetypeout)         

def getbamcountsthr(results,currsample, *args, **kwargs):
    results[currsample] = getbamcounts(*args, **kwargs)
def getbamcountsqueue(countqueue,currsample, *args, **kwargs):
    countqueue.put([currsample,getbamcounts(*args, **kwargs)])
    
def countreadspool(args):
    return getbamcounts(*args[0], **args[1])

trnaends = list(["CCA","CC","C",""])    
def printtrnaendfile(trnaendfilename,samples,  samplecounts, trnalist, trnaloci , minreads = 5):
    trnaendfile = open(trnaendfilename, "w")
    print("end\t"+"\t".join(currsample for currsample in samples), file=trnaendfile)

    for currfeat in trnalist:
        if max(samplecounts[currsample].gettrnacount(currfeat.name) for currsample in samples) < minreads:
            continue
        for currend in trnaends:
            endstring = currend
            if currend == "":
                endstring = "Trimmed"
            print(currfeat.name+"\t"+endstring+"\t"+"\t".join(str(samplecounts[currsample].getendtypecount(currfeat.name)[currend]) for currsample in samples), file=trnaendfile)
    trnaendfile.close()   

def compressargs( *args, **kwargs):
    return tuple([args, kwargs])

def countreads_main(**argdict):
    trnauniquefilename = None
    argdict = defaultdict(lambda: None, argdict)
    includebase = argdict["nofrag"]
    trnatable = argdict["trnatable"]
    removepseudo = argdict["removepseudo"]
    ensemblgtf = argdict["ensemblgtf"]
    if argdict["maxmismatches"] is not None:
        maxmismatches = int(argdict["maxmismatches"])
    else:
        maxmismatches = None
    cores = argdict["cores"]
    trnaendfilename = argdict["trnaends"]
    threadmode = True
    if cores == 1:
        threadmode = False
    otherseqs = toolsTG.extraseqfile(argdict["otherseqs"])
    
    if "bamdir" not in argdict:
        bamdir = "./"
    bamdir = argdict["bamdir"]
    sampledata = toolsTG.samplefile(argdict["samplefile"], bamdir = bamdir)
    
    bedfiles = list() 
    if "trnauniquecounts" in argdict:
        trnauniquefilename = argdict["trnauniquecounts"]
    if "bedfile"  in argdict:
        bedfiles = argdict["bedfile"]
    trnalocifiles = list()
    if "trnaloci"  in argdict:
        trnalocifiles = argdict["trnaloci"]
    maturetrnas = list()
    if "maturetrnas" in argdict:
        maturetrnas = argdict["maturetrnas"]
        
    #trnalocifiles = argdict["trnaloci"]
    #maturetrnas=argdict["maturetrnas"]
    genetypefile = argdict["genetypefile"]
    trnacountfilename = argdict["trnacounts"]

    
    trnainfo = toolsTG.transcriptfile(trnatable)

    #print >>sys.stderr, bedfiles
    
    
    samples = sampledata.getsamples()
    genetypes = dict()
    otherseqdict = dict()
    #Grabbing all the features to count
    #print >>sys.stderr, otherseqs
    try:
        featurelist = dict()
        trnaloci = list()
        for currfile in bedfiles:
            bedfeatures = list(toolsTG.readfeatures(currfile, removepseudo = removepseudo))
            for curr in bedfeatures:
                genetypes[curr.name] = os.path.basename(currfile)
                
            featurelist[currfile] = bedfeatures
        trnalist = list()
        for currfile in trnalocifiles:
            trnaloci.extend(list(toolsTG.readbed(currfile)))
        for currfile in maturetrnas:
            trnalist.extend(list(toolsTG.readbed(currfile)))
        if ensemblgtf is not None:    
            embllist = list(toolsTG.readgtf(ensemblgtf, filterpsuedo = removepseudo))
        else:
            embllist = list()
        for name, currfile in otherseqs.getseqbeds().items():
            otherseqdict[name] = list(toolsTG.readbed(currfile))

    except IOError as e:
        print(e, file=sys.stderr)
        sys.exit()
    allfeats = trnaloci+trnalist
    if len(set(curr.name for curr in allfeats)) < len(list(curr.name for curr in allfeats )):
        print("Duplicate names in feature list", file=sys.stderr)
    
    
    #featurelist = list(curr for curr in featurelist if curr.name == 'unknown20') 
    #alltrnas = list(curr.name for curr in featurelist)
    #print >>sys.stderr, "***"
    #setting up all the feature count dictionaries
                            
    allcounts = dict()
    #threadmode = False
    #print  list(curr.name for curr in trnalist)
    #print >>sys.stderr, "**||"
    #print >>sys.stderr, maxmismatches
    
    #sys.exit()
    if threadmode:

        countpool = Pool(processes=cores)
        arglist = list()
        for currsample in samples:
            currbam = sampledata.getbam(currsample)
            arglist.append(compressargs(currbam, currsample,trnainfo, trnaloci, trnalist, otherseqdict = otherseqdict, embllist = embllist, featurelist = featurelist, maxmismatches = maxmismatches))
        #arglist = list((tuple([currsample, sampledata.getbam(currsample)]) for currsample in samples))
        results = countpool.map(countreadspool, arglist)
        for i, curr in enumerate(samples):
            allcounts[curr] = results[i]

    else:

        for currsample in samples:
            
            currbam = sampledata.getbam(currsample)
            allcounts[currsample] = getbamcounts(currbam, currsample,trnainfo, trnaloci, trnalist, otherseqdict = otherseqdict,embllist = embllist, featurelist = featurelist, bedfiles = bedfiles, maxmismatches = maxmismatches)
            #getbamcountsthr(allcounts, allcounts)
            #threads[currsample] = threading.Thread(target=getbamcountsthr, args=(allcounts,currsample,currbam, currsample,trnainfo, trnaloci, trnalist), kwargs = {'embllist' : embllist, 'featurelist' : featurelist, 'maxmismatches' : maxmismatches})
            #threads[currsample].start()
    #print >>sys.stderr, "time:" +str(endtime-starttime)
    if "countfile" not in argdict or argdict["countfile"] == "stdout":
        countfile = sys.stdout
    else:
        countfile = open(argdict["countfile"], "w")
    printcountfile(countfile, samples, allcounts,trnalist, trnaloci, featurelist, embllist, otherseqdict = otherseqdict,includebase = includebase)
            

    if genetypefile is not None:
        genetypeout = open(genetypefile, "w")
        printgenetypes(genetypeout,samples, allcounts,trnalist, trnaloci, featurelist, embllist,otherseqdict = otherseqdict )
    #it's currently not used, but here is where I could count by amino acid or anticodon
    #if typefile:
    #    trnacountfile = open(trnacountfilename, "w")
    #    for curramino in trnainfo.allaminos():
    #            print("AminoTotal_"+curramino+"\t"+"\t".join(str(aminocounts[currsample][curramino]) for currsample in samples), file=typefile)
    #    for currac in trnainfo.allanticodons():
    #            print("AnticodonTotal_"+currac+"\t"+"\t".join(str(anticodoncounts[currsample][currac]) for currsample in samples), file=typefile)


            
    if trnacountfilename is not None:
        #trnauniquefile = open(trnauniquefilename, "w")
        printtrnacountfile(trnacountfilename, samples, allcounts,trnalist,includebase = includebase)
        printtrnaendfile(trnaendfilename,samples,  allcounts, trnalist, trnaloci)
       
        
    if trnauniquefilename is not None:
        
        #trnauniquefile = open(trnauniquefilename, "w")
        printtrnauniquecountcountfile(trnauniquefilename,samples,  allcounts, trnalist, trnaloci,fragsep = includebase )
        #print >>trnauniquefile, "\t".join(currsample for currsample in samples)
        #for currfeat in trnalist:
        #    if max(trnacounts[currsample][currfeat.name] for currsample in samples) < minreads:
        #        continue
        #    print  >>trnauniquefile, currfeat.name+"\t"+"\t".join(str(trnacounts[currsample][currfeat.name]) for currsample in samples)
        #trnauniquefile.close()          
        pass

def printtrnauniquecountcountfile(trnauniquefile,samples,  samplecounts, trnalist, trnaloci , minreads = 5, fragsep = True):
    trnauniquefile = open(trnauniquefile, "w")
    print("\t".join(currsample for currsample in samples), file=trnauniquefile)
    for currfeat in trnalist:
        if max(samplecounts[currsample].getuniquecount(currfeat.name) for currsample in samples) < minreads:
            continue
        if fragsep:
            print(currfeat.name+"_fiveprime\t"+"\t".join(str(samplecounts[currsample].getfiveuniquecount(currfeat.name)) for currsample in samples), file=trnauniquefile)
            print(currfeat.name+"_threeprime\t"+"\t".join(str(samplecounts[currsample].getthreeuniquecount(currfeat.name)) for currsample in samples), file=trnauniquefile)
            print(currfeat.name+"_whole\t"+"\t".join(str(samplecounts[currsample].getwholeuniquecount(currfeat.name)) for currsample in samples), file=trnauniquefile)
            print(currfeat.name+"_other\t"+"\t".join(str(samplecounts[currsample].getotheruniquecount(currfeat.name)) for currsample in samples), file=trnauniquefile)
        else:
            print(currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].getuniquecount(currfeat.name)) for currsample in samples), file=trnauniquefile)
    trnauniquefile.close() 

def printrealcounts(countfile,samples, sampledata,allcounts,trnalist, trnaloci, bedtypes, emblbiotypes,extraseqtypes = set()):

    biotypefirst = ['snoRNA','snRNA','scaRNA','sRNA','miRNA']         
    biotypelast = ['Mt_rRNA','Mt_tRNA','rRNA']
    otherbiotypes = list(set(emblbiotypes) - (set(biotypefirst) | set(biotypelast)))
    biotypeorder = biotypefirst + otherbiotypes + biotypelast
        
    print("\t".join(samples), file=countfile)
    
    print("other"+"\t"+"\t".join(str(allcounts[currsample].otherreads) for currsample in samples), file=countfile)
    for currbed in bedtypes:
        print(os.path.basename(currbed)+"\t"+"\t".join(str(allcounts[currsample].bedtypecounts[currbed]) for currsample in samples), file=countfile)
    for currname in extraseqtypes:
         print(currname+"_seq\t"+"\t".join(str(allcounts[currsample].extraseqcounts[currname]) for currsample in samples), file=countfile)

    for currbiotype in reversed(biotypeorder):
        print(currbiotype+"\t"+"\t".join(str(allcounts[currsample].embltypecounts[currbiotype]) for currsample in samples), file=countfile)
        
    for currbed in trnaloci:
 
        print("pretRNA\t"+"\t".join(str(allcounts[currsample].trnalocuscounts[currbed]) for currsample in samples), file=countfile)
    for currbed in trnalist:     
        print("tRNA_antisense\t"+"\t".join(str(allcounts[currsample].trnaanticounts[currbed]) for currsample in samples), file=countfile)
        print("tRNA\t"+"\t".join(str(allcounts[currsample].trnacounts[currbed]) for currsample in samples), file=countfile)
        

def printaminocounts(trnaaminofilename, sampledata,trnainfo,allcounts, sizefactor, uniquemode = False, fragmode = True):
    trnaaminofile = open(trnaaminofilename, "w")
    
    aminos = trnainfo.allaminos()
    repmode = False
    if repmode:
        replicates = list(sampledata.allreplicates())
        
        print("\t".join(replicates), file=trnaaminofile)
        for curramino in aminos:
            print(curramino+"\t"+"\t".join(str(sum(allcounts[currsample].aminocounts[curramino]/sizefactor[currsample] for currsample in sampledata.getrepsamples(currrep))) for currrep in replicates), file=trnaaminofile)
    else:
        
        allsamples = list(sampledata.getsamples())
        
        
        print("\t".join(allsamples), file=trnaaminofile)
        for curramino in aminos:
            if uniquemode:
                if fragmode:
                    for fragtype in allfragtypes:
                        print(curramino+"_"+fragnames[fragtype]+"\t"+"\t".join(str(allcounts[currsample].aminouniqcounts[fragtype][curramino]) for currsample in allsamples), file=trnaaminofile)
                else:
                    print(curramino+"\t"+"\t".join(str(allcounts[currsample].aminoallfrags(curramino)) for currsample in allsamples), file=trnaaminofile)
            else:
                if fragmode:
                    for fragtype in allfragtypes:
                        print(curramino+"_"+fragnames[fragtype]+"\t"+"\t".join(str(allcounts[currsample].aminotranscriptcounts[fragtype][curramino]) for currsample in allsamples), file=trnaaminofile)
                else:
                    print(curramino+"\t"+"\t".join(str(allcounts[currsample].aminoallfrags(curramino, uniq = False)) for currsample in allsamples), file=trnaaminofile)
                    
def printanticodoncounts(trnaanticodonfilename, sampledata,trnainfo,allcounts, sizefactor, uniquemode = False, fragmode = True):
    anticodons = trnainfo.allanticodons()
    trnaanticodonfile = open(trnaanticodonfilename, "w")
    repmode = False    
    if repmode:
        replicates = list(sampledata.allreplicates())
        
        print("\t"+"\t".join(replicates), file=trnaanticodonfile)
        for curranticodon in anticodons:
            print(curranticodon+"\t"+"\t".join(str(sum(allcounts[currsample].anticodoncounts[curranticodon]/sizefactor[currsample] for currsample in sampledata.getrepsamples(currrep))) for currrep in replicates), file=trnaanticodonfile)
    else:
        allsamples = list(sampledata.getsamples())
        
        
        print("\t".join(allsamples), file=trnaanticodonfile)
        for curranticodon in anticodons:
            if uniquemode:
                if fragmode:
                    for fragtype in allfragtypes:
                        print(curranticodon+"_"+fragnames[fragtype]+"\t"+"\t".join(str(allcounts[currsample].anticodonuniqcounts[fragtype][curranticodon]) for currsample in allsamples), file=trnaanticodonfile)
                else:
                    print(curranticodon+"\t"+"\t".join(str(allcounts[currsample].anticodonallfrags(curranticodon)) for currsample in allsamples), file=trnaanticodonfile)
            else:
                if fragmode:
                    for fragtype in allfragtypes:
                        print(curranticodon+"_"+fragnames[fragtype]+"\t"+"\t".join(str(allcounts[currsample].anticodontranscriptcounts[fragtype][curranticodon]) for currsample in allsamples), file=trnaanticodonfile)
                else:
                    print(curranticodon+"\t"+"\t".join(str(allcounts[currsample].anticodonallfrags(curranticodon, uniq = False)) for currsample in allsamples), file=trnaanticodonfile)
                    
                
                
def printtrnacounts(trnacountfilename, sampledata,trnainfo,allcounts, sizefactor, uniquemode = False, fragmode = True):
    alltrnas = trnainfo.gettranscripts()
    trnacountfile = open(trnacountfilename, "w")
    repmode = False
    if repmode:
        replicates = list(sampledata.allreplicates())
        
        print("\t"+"\t".join(replicates), file=trnacountfile)
        # Assuming we want to iterate over trnas here if repmode was True
        for currtrna in alltrnas:
             print(currtrna+"\t"+"\t".join(str(sum(allcounts[currsample].trnatranscriptuniqcounts[currtrna]/sizefactor[currsample] for currsample in sampledata.getrepsamples(currrep))) for currrep in replicates), file=trnacountfile)
    else:
        allsamples = list(sampledata.getsamples())
        
        
        print("\t".join(allsamples), file=trnacountfile)
        for currtrna in alltrnas:
            if uniquemode:
                if fragmode:
                    for fragtype in allfragtypes:
                        print(currtrna+"_"+fragnames[fragtype]+"\t"+"\t".join(str(allcounts[currsample].trnatranscriptuniqcounts[fragtype][currtrna]) for currsample in allsamples), file=trnacountfile)
                else:
                    print(currtrna+"\t"+"\t".join(str(allcounts[currsample].trnaallfrags(currtrna)) for currsample in allsamples), file=trnacountfile)

            else:
                if fragmode:
                    for fragtype in allfragtypes:
                        print(currtrna+"_"+fragnames[fragtype]+"\t"+"\t".join(str(allcounts[currsample].trnatranscriptcounts[fragtype][currtrna]) for currsample in allsamples), file=trnacountfile)
                else:
                    print(currtrna+"\t"+"\t".join(str(allcounts[currsample].trnaallfrags(currtrna, uniq = False)) for currsample in allsamples), file=trnacountfile)


def printmismatchcounts(trnamismatchname, sampledata,trnainfo,allcounts, sizefactor):
    trnamismatchfile = open(trnamismatchname, "w")
    repmode = False
    mismatchcounts = list(range(10))
    if repmode:# not tested yet
        replicates = list(sampledata.allreplicates())
        
        print("count\ttype\t"+"\t".join(replicates), file=trnamismatchfile)
        for currmismatch in mismatchcounts:
            
            print(str(currmismatch)+"\ttrna\t"+"\t".join(str(sum(allcounts[currsample].trnamismatchcounts[currmismatch] for currsample in sampledata.getrepsamples(currrep))) for currrep in replicates), file=trnamismatchfile)
            print(str(currmismatch)+"\tnontrna\t"+"\t".join(str(sum(allcounts[currsample].mismatchcounts[currmismatch] for currsample in sampledata.getrepsamples(currrep)) - sum(allcounts[currsample].trnamismatchcounts[currmismatch] for currsample in sampledata.getrepsamples(currrep))) for currrep in replicates), file=trnamismatchfile)

    else:
        allsamples = list(sampledata.getsamples())
        
        
        print("count\ttype\t"+"\t".join(allsamples), file=trnamismatchfile)
        for currmismatch in mismatchcounts:
            print(str(currmismatch)+"\ttrna\t"+"\t".join(str(allcounts[currsample].trnamismatchcounts[currmismatch]/sizefactor[currsample]) for currsample in allsamples), file=trnamismatchfile)
            print(str(currmismatch)+"\tnontrna\t"+"\t".join(str(allcounts[currsample].mismatchcounts[currmismatch]/sizefactor[currsample] - allcounts[currsample].trnamismatchcounts[currmismatch]/sizefactor[currsample]) for currsample in allsamples), file=trnamismatchfile)



def printtrnanormfile(trnanormfile, samples, allcounts):         
    trnasamplecounts = {sample: sum(allcounts[sample].trnacounts.values()) for sample in samples}
    trnanormfile = open(trnanormfile, "w")
    mean = 1.*sum(trnasamplecounts.values())/len(list(trnasamplecounts.values()))
    print("\t".join(samples), file=trnanormfile)
    print("\t".join(str(trnasamplecounts[currsample]/mean) for currsample in samples), file=trnanormfile)

def printallreadsnormfile(allreadsnormfile, samples, allcounts):
    totalsamplecounts = {sample: sum(allcounts[sample].counts.values()) + sum(allcounts[sample].trnacounts.values()) for sample in samples}
    allreadsnormfile = open(allreadsnormfile, "w")
    mean = 1.*sum(totalsamplecounts.values())/len(list(totalsamplecounts.values()))
    print("\t".join(samples), file=allreadsnormfile)
    print("\t".join(str(totalsamplecounts[currsample]/mean) for currsample in samples), file=allreadsnormfile)

def printlengthfile(readlengthfile, samples,allcounts):
    readlengthfile = open(readlengthfile, "w")
    print("Length\tSample\tother\ttrnas\tpretrnas", file=readlengthfile)
    for currsample in samples:
        for curr in range(0,max(allcounts[currsample].readlengths.keys())+1):
            othercount = allcounts[currsample].trnareadlengths[curr] + allcounts[currsample].pretrnareadlengths[curr]
            print(str(curr)+"\t"+currsample+"\t"+str(allcounts[currsample].readlengths[curr] - othercount)+"\t"+str(allcounts[currsample].trnareadlengths[curr]) +"\t"+str(allcounts[currsample].pretrnareadlengths[curr]), file=readlengthfile)

def printtypefile(countfile,samples, sampledata,allcounts,trnalist, trnaloci, bedtypes, emblbiotypes, sizefactor,extraseqtypes = set(),countfrags = False, combinereps = True):

    def sumsamples(countdict,sampledata, repname, currfeat = None, sizefactors = defaultdict(lambda: 1)):
        return sum(countdict[currsample]/sizefactors[currsample] for currsample in sampledata.getrepsamples(repname))
  
     
    
    if combinereps:
        replicates = list(sampledata.allreplicates())
        print("\t".join(replicates), file=countfile)
        print("other"+"\t"+"\t".join(str(sumsamples({s: allcounts[s].otherreads for s in samples},sampledata,currrep, sizefactors = sizefactor)) for currrep in replicates), file=countfile)
        for currbed in bedtypes:
            print(os.path.basename(currbed)+"\t"+"\t".join(str(sumsamples({s: allcounts[s].bedtypecounts[currbed] for s in samples},sampledata,currrep, sizefactors = sizefactor)) for currrep in replicates), file=countfile)
        for currname in extraseqtypes:
             print(currname+"_seq\t"+"\t".join(str(sumsamples({s: allcounts[s].extraseqcounts[currname] for s in samples},sampledata,currrep, sizefactors = sizefactor)) for currrep in replicates), file=countfile)

        biotypefirst = ['snoRNA','snRNA','scaRNA','sRNA','miRNA']         
        biotypelast = ['Mt_rRNA','Mt_tRNA','rRNA']
        otherbiotypes = list(set(emblbiotypes) - (set(biotypefirst) | set(biotypelast)))
        biotypeorder = biotypefirst + otherbiotypes + biotypelast
        for currbiotype in reversed(biotypeorder):
            print(currbiotype+"\t"+"\t".join(str(sumsamples({s: allcounts[s].embltypecounts[currbiotype] for s in samples},sampledata,currrep, sizefactors = sizefactor)) for currrep in replicates), file=countfile)
            
        for currbed in trnaloci:
     
            print("pretRNA_antisense\t"+"\t".join(str(sumsamples({s: allcounts[s].trnaantilocuscounts[currbed] for s in samples},sampledata,currrep, sizefactors = sizefactor)) for currrep in replicates), file=countfile)
            print("pretRNA\t"+"\t".join(str(sumsamples({s: allcounts[s].trnalocuscounts[currbed] for s in samples},sampledata,currrep, sizefactors = sizefactor)) for currrep in replicates), file=countfile)
        for currbed in trnalist:     
            print("tRNA_antisense\t"+"\t".join(str(sumsamples({s: allcounts[s].trnaanticounts[currbed] for s in samples},sampledata,currrep, sizefactors = sizefactor)) for currrep in replicates), file=countfile)
            print("tRNA\t"+"\t".join(str(sumsamples({s: allcounts[s].trnacounts[currbed] for s in samples},sampledata,currrep, sizefactors = sizefactor)) for currrep in replicates), file=countfile)
            
        
    else:
        pass

def main(**argdict):
    argdict = defaultdict(lambda: None, argdict)
    countfrags = argdict["countfrags"]
    ensemblgtf = argdict["ensemblgtf"]
    bamnofeature = argdict["bamnofeature"]
    trnatable = argdict["trnatable"]
    uniquename = argdict["uniquename"]
    fraguniq = argdict["fraguniq"]
    otherseqs = toolsTG.extraseqfile(argdict["otherseqs"])
    
    if "bamdir" not in argdict:
        bamdir = "./"
    bamdir = argdict["bamdir"]
    sampledata = toolsTG.samplefile(argdict["samplefile"], bamdir = bamdir)
    cores = argdict["cores"]
    threadmode = True
    if cores == 1:
        threadmode = False
    
    sizefactor = defaultdict(lambda: 1)
    if argdict["sizefactors"]:
        sizefactor = toolsTG.getsizefactors(argdict["sizefactors"]) 
        for currsample in sampledata.getsamples():
            if currsample not in sizefactor:
                print("Size factor file "+argdict["sizefactors"]+" missing "+currsample, file=sys.stderr)
                sys.exit(1)
        
    bedfiles = list()
    
    if argdict["bedfile"]  is not None:
        bedfiles = argdict["bedfile"]
    
    locifiles = argdict["trnaloci"]
    maturetrnafiles = argdict["maturetrnas"]     
    trnaaminofilename = argdict["trnaaminofile"]
    trnaanticodonfilename = argdict["trnaanticodonfile"]
    readlengthfile = argdict["readlengthfile"]
    mismatchfilename = argdict["mismatchfile"]    
    
    if argdict["realcountfile"] == "stdout":
        realcountfile = sys.stdout
    else:
        realcountfile = open(argdict["realcountfile"],"w")
    
    
    if argdict["countfile"] == "stdout":
        countfile = sys.stdout
    else:
        countfile = open(argdict["countfile"],"w")
    
    maturenames = dict()
    otherseqlist = dict()
    
    trnainfo = toolsTG.transcriptfile(trnatable)
    
    samples = list(sampledata.getsamples())
    
    try:
        featurelist = dict()
        trnaloci = dict()
        trnalist = dict()
        otherseqlist = dict()
        for currfile in bedfiles:
            featurelist[currfile] = toolsTG.RangeBin(toolsTG.readfeatures(currfile))
        
        for currfile in locifiles:
            trnaloci[currfile] = toolsTG.RangeBin(toolsTG.readbed(currfile), binfactor = 10000)
        for currfile in maturetrnafiles:
            matlist = list(toolsTG.readbed(currfile))
            trnalist[currfile] = list(matlist)
            maturenames[currfile] = {curr.name:curr for curr in matlist}
        if ensemblgtf is not None:    
            embllist = toolsTG.RangeBin(toolsTG.readgtf(ensemblgtf, filtertypes = set()))
        else:
            embllist = None
        
        for currname, currfile in otherseqs.getseqbeds().items():
            otherseqlist[currname] = getchromdict(toolsTG.readbed(currfile))
    except IOError as e:
        print(e, file=sys.stderr)
        sys.exit()
    
    maxmismatches = None
    allcounts = dict()
    poolmode = True
    
    if threadmode:
        countqueue = Queue()
        threads = dict()
        if poolmode:
            countpool = Pool(processes=cores)
            arglist = list()
            for currsample in samples:
                currbam = sampledata.getbam(currsample)
                arglist.append(compressargs(currbam,currsample, trnainfo, trnaloci, trnalist,maturenames, otherseqlist = otherseqlist, embllist = embllist, featurelist = featurelist, maxmismatches = maxmismatches, bamnofeature = bamnofeature))
            results = countpool.map(counttypereadspool, arglist)
            for i, curr in enumerate(samples):
                allcounts[curr] = results[i]
        else:
            for currsample in samples:
                currbam = sampledata.getbam(currsample)
                threads[currsample] = Process(target=counttypereadsqueue,args = (countqueue,currsample,currbam, currsample,trnainfo, trnaloci, trnalist,maturenames), kwargs = { "embllist" : embllist, "featurelist" : featurelist, "maxmismatches" : maxmismatches, "bamnofeature" : bamnofeature})
                threads[currsample].start()
            for sample in threads.keys():
                currsample, counts = countqueue.get()
                allcounts[currsample] = counts
            
    else:
        for currsample in samples:
            currbam = sampledata.getbam(currsample)
            allcounts[currsample] = counttypereads(currbam, currsample,trnainfo, trnaloci, trnalist,maturenames, otherseqlist = otherseqlist, embllist = embllist, featurelist = featurelist, maxmismatches = maxmismatches, bamnofeature = bamnofeature)
        
    emblbiotypes  = set(itertools.chain.from_iterable(curr.emblbiotypes for curr in list(allcounts.values())))        
    bedtypes  = set(itertools.chain.from_iterable(curr.bedtypes for curr in list(allcounts.values()))) 
    extraseqtypes  = set(itertools.chain.from_iterable(curr.extraseqtypes for curr in list(allcounts.values())))   
    
    printtypefile(countfile, samples, sampledata,allcounts,trnalist, trnaloci, bedtypes, emblbiotypes,sizefactor, countfrags = countfrags , extraseqtypes = extraseqtypes)
    printrealcounts(realcountfile, samples, sampledata,allcounts,trnalist, trnaloci, bedtypes, emblbiotypes , extraseqtypes = extraseqtypes)

    if readlengthfile is not None:
        printlengthfile(readlengthfile, samples, allcounts)

    if trnaaminofilename is not None:
        printaminocounts(trnaaminofilename, sampledata,trnainfo, allcounts, sizefactor, fragmode = False)
    if trnaanticodonfilename is not None:
        printanticodoncounts(trnaanticodonfilename, sampledata,trnainfo, allcounts, sizefactor, fragmode = False)
    if mismatchfilename is not None:
        printmismatchcounts(mismatchfilename, sampledata,trnainfo, allcounts, sizefactor)
    if uniquename is not None:
        printaminocounts(uniquename+"-aminos.txt", sampledata,trnainfo, allcounts, sizefactor, uniquemode = True, fragmode =fraguniq)
        printanticodoncounts(uniquename+"-anticodons.txt", sampledata,trnainfo, allcounts, sizefactor, uniquemode = True, fragmode =fraguniq)
        printtrnacounts(uniquename+"-trnas.txt", sampledata,trnainfo, allcounts, sizefactor, uniquemode = True, fragmode = fraguniq)
