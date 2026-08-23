#!/usr/bin/env python3

import sys
import os
import subprocess
import tempfile
import time
import logging
import pysam
from pathlib import Path
from typing import List, Optional, Union, Any
from multiprocessing import Pool
from shutil import which

# Try to import from local toolsTG, otherwise assume it's in the path
try:
    from .toolsTG import samplefile, getsizefactors, transcriptfile, getnamedict, readbed, revcom, sort_temp_dir
except ImportError:
    try:
        from toolsTG import samplefile, getsizefactors, transcriptfile, getnamedict, readbed, revcom, sort_temp_dir
    except ImportError:
        # Fallback for when running from tRNAgraph directory
        sys.path.append(os.path.dirname(__file__))
        from toolsTG import samplefile, getsizefactors, transcriptfile, getnamedict, readbed, revcom, sort_temp_dir

def cigarlength(cigar):
    # 0=M, 1=I. Counts length on query sequence.
    return sum(curr[1] for curr in cigar if curr[0] in {0, 1})

def reverseintrons(introns, length):
    for currint in reversed(introns):
        yield (length - currint[0] + 1, currint[1])

class TrackHubBuilder:
    def __init__(self, 
                 genomedatabase: str, 
                 samplefilename: str, 
                 expname: str, 
                 threads: int = 8):
        self.dbname = genomedatabase
        self.samplefilename = samplefilename
        self.expname = expname
        self.threads = threads
        self.logger = logging.getLogger(__name__)

        self.sampledata = samplefile(self.samplefilename)
        self.trackdir = Path(self.expname) / "trackhub"
        
        size_factors_file = f"{self.expname}/{self.expname}-SizeFactors.txt"
        if os.path.exists(size_factors_file):
            self.sizefactors = getsizefactors(size_factors_file)
        else:
            self.sizefactors = {}

    def get_location(self, program: str, allowfail: bool = False) -> Optional[str]:
        progloc = which(program)
        if progloc is None and not allowfail:
            self.logger.error(f"Could not find {program} in path")
            self.logger.error("Aborting")
            sys.exit(1)
        return progloc

    def convertbam(self, inputbam: str, outputbam: str, force: bool = False, logfile: Optional[Any] = sys.stderr) -> None:
        if not os.path.isfile(outputbam) or force:
            tempprefix = f"{Path(inputbam).stem}_{os.getpid()}"
            temp_dir = sort_temp_dir(outputbam)
            
            self.logger.info(f"Converting {inputbam} to genome coordinates...")
            
            # Prepare sorting process
            sort_cmd = f"samtools sort -T {temp_dir}/convert{tempprefix} - -o {outputbam}"
            sort_process = subprocess.Popen(sort_cmd, shell=True, stdin=subprocess.PIPE)
            
            try:
                # Open input BAM
                bamfile = pysam.AlignmentFile(inputbam, "rb")
                
                # Open output BAM stream to sort process
                outfile = pysam.AlignmentFile(sort_process.stdin, "wb", template=bamfile)
                
                # Load tRNA info
                trnainfo = transcriptfile(f"{self.dbname}-trnatable.txt")
                trnaloci = getnamedict(readbed(f"{self.dbname}-trnaloci.bed", includeintrons=True))
                trnatranscripts = trnainfo.gettranscripts()
                
                npad = 20
                bam_cins = 1
                bam_cdel = 2
                
                for currmap in bamfile:
                    chromname = bamfile.getrname(currmap.tid)
                    readquals = currmap.query_qualities
                    readseq = currmap.query_sequence
                    readtags = currmap.get_tags()
                    origcigar = currmap.cigartuples
                    
                    if chromname in trnatranscripts:
                        origstart = currmap.reference_start
                        
                        # Iterate over loci for this transcript
                        for currlocusname in trnainfo.transcriptdict[chromname]:
                            currlocus = trnaloci[currlocusname]
                            
                            # We need a new read object to write
                            new_read = pysam.AlignedSegment(bamfile.header)
                            new_read.query_name = currmap.query_name
                            new_read.query_sequence = currmap.query_sequence
                            new_read.query_qualities = currmap.query_qualities
                            new_read.flag = currmap.flag
                            new_read.reference_id = bamfile.get_tid(currlocus.chrom)
                            new_read.mapping_quality = currmap.mapping_quality
                            new_read.tags = readtags
                            new_read.next_reference_id = currmap.next_reference_id
                            new_read.next_reference_start = currmap.next_reference_start
                            new_read.template_length = currmap.template_length
                            
                            introns = list()
                            if currlocus.data and "blockcount" in currlocus.data and currlocus.data["blockcount"] > 1:
                                lastblock = 0
                                for i in range(currlocus.data["blockcount"]):
                                    currblocksize = int(currlocus.data["blocksizes"][i])
                                    currblockstart = int(currlocus.data["blockstarts"][i])
                                    if lastblock != 0:
                                        introns.append([lastblock, currblockstart - lastblock])
                                    lastblock += currblocksize
                                    
                            if currlocus.strand == '+':
                                new_read.reference_start = origstart - npad + currlocus.start
                                new_read.query_sequence = readseq
                                new_read.query_qualities = readquals
                                new_read.is_reverse = False
                                
                                currpoint = origstart - npad
                                newcigar = list()
                                
                                # Adjust start for introns
                                for i in range(len(introns)):
                                    if currpoint >= introns[i][0]:
                                        new_read.reference_start += introns[i][1]
                                        currpoint += introns[i][1]
                                        
                                for currcigar in origcigar:
                                    if currcigar[0] in {0, bam_cins}: # M or I
                                        foundintron = False
                                        for intronstart, intronlength in introns:
                                            if currpoint < intronstart < currpoint + currcigar[1]:
                                                firseglength = intronstart - currpoint
                                                secseglength = currpoint + currcigar[1] - intronstart
                                                newcigar.append((currcigar[0], firseglength))
                                                newcigar.append((bam_cdel, intronlength))
                                                newcigar.append((currcigar[0], secseglength))
                                                foundintron = True
                                        if not foundintron:
                                            newcigar.append(currcigar)
                                        if currcigar[0] in {0, bam_cdel}: # M or D (D consumes ref)
                                            currpoint += currcigar[1]
                                    else:
                                        newcigar.append(currcigar)
                                        
                                new_read.cigartuples = newcigar
                                
                            else: # Minus strand
                                new_read.reference_start = currlocus.end - (origstart - npad + cigarlength(origcigar) + sum(curr[1] for curr in introns))
                                new_read.query_sequence = revcom(readseq)
                                new_read.query_qualities = list(reversed(readquals))
                                new_read.is_reverse = True
                                
                                currpoint = origstart - npad
                                newcigar = list()
                                
                                for intronstart, intronlength in reverseintrons(introns, currlocus.length()):
                                    if currpoint + cigarlength(origcigar) < intronstart:
                                        new_read.reference_start += intronlength
                                        currpoint -= intronlength
                                        
                                for currcigar in reversed(origcigar):
                                    if currcigar[0] in {0, bam_cins}:
                                        foundintron = False
                                        for intronstart, intronlength in introns:
                                            if currpoint < intronstart < currpoint + currcigar[1]:
                                                firseglength = intronstart - currpoint + 1
                                                secseglength = currpoint + currcigar[1] - intronstart - 1
                                                newcigar.append((currcigar[0], secseglength))
                                                newcigar.append((bam_cdel, intronlength))
                                                newcigar.append((currcigar[0], firseglength))
                                                foundintron = True
                                        if not foundintron:
                                            newcigar.append(currcigar)
                                        if currcigar[0] in {0, bam_cdel}:
                                            currpoint += currcigar[1]
                                    else:
                                        newcigar.append(currcigar)
                                        
                                new_read.cigartuples = newcigar
                                
                            outfile.write(new_read)
                    else:
                        outfile.write(currmap)
                        
                outfile.close()
                bamfile.close()
                sort_process.wait()
                
                if sort_process.returncode != 0:
                    self.logger.error("Failure to convert bam to genome space")
                    sys.exit(1)
                
                subprocess.check_call(f"samtools index {outputbam}", shell=True, stderr=logfile)
                
            except Exception as e:
                self.logger.error(f"Error converting BAM: {e}")
                if sort_process.poll() is None:
                    sort_process.kill()
                sys.exit(1)

        else:
            self.logger.warning(f"Skipping {outputbam} (already exists)")

    def samtoolsmerge(self, bamfiles: List[str], outbam: str, force: bool = False) -> None:
        samtoolsloc = self.get_location("samtools")
        if not os.path.isfile(outbam) or force:
            samcommand = [samtoolsloc, "merge", "-f", outbam] + bamfiles
            
            self.logger.info(" ".join(samcommand))
            result = subprocess.run(samcommand, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, check=False)
            if result.stdout:
                self.logger.info(result.stdout)

    def createmultiwigtrackdb(self, trackfile: Any, shortlabel: str = "", longlabel: str = "", suffix: str = '', startpriority: float = 3.0, stacked: bool = False) -> None:
        trackcolors = ['0,217,47', '47,142,248', '220,21,235', '264,115,6', '95,238,230']
        currpriority = startpriority
        
        trackfile.write(f"track {self.expname}{suffix}tracks            \n")
        trackfile.write("superTrack on show\n")
        trackfile.write(f"shortLabel {self.expname} {shortlabel}\n")
        trackfile.write(f"longLabel Data from {self.expname} {longlabel}\n")
        trackfile.write("visibility full\n")
        trackfile.write("\n")
        
        for currrep in self.sampledata.allreplicates():
            for currstrand in ['Plus', 'Minus']:
                trackfile.write(f"\ttrack {currrep}{suffix}_{currstrand}tracks\n")
                trackfile.write("\tcontainer multiWig\n")
                trackfile.write(f"\tshortLabel {currrep}{suffix} {currstrand} Strand\n")
                trackfile.write(f"\tlongLabel Data from {self.expname} {currrep}{suffix} {currstrand} Strand\n")
                trackfile.write("\ttype bigWig\n")
                trackfile.write(f"\tparent {self.expname}{suffix}tracks on\n")
                trackfile.write("\tdragAndDrop on\n")
                
                if stacked:
                    trackfile.write("\taggregate solidOverlay\n")
                else:
                    trackfile.write("\taggregate transparentOverlay\n")
                    
                trackfile.write("\tshowSubtrackColorOnUi on\n")
                trackfile.write("\tautoScale on\n")
                trackfile.write("\talwaysZero on\n")
                trackfile.write(f"\tpriority {currpriority + 0.1}  \n")
                trackfile.write("\tmaxHeightPixels 256:100:32\n")
                trackfile.write("\tvisibility full\n")
                trackfile.write("\n")
                
                currpriority += 0.2
                repsamples = self.sampledata.getrepsamples(currrep)
                for i, currsample in enumerate(repsamples):
                    trackfile.write(f"\t\ttrack {currsample}{suffix}_{currstrand}track\n")
                    trackfile.write("\t\ttype bigWig\n")
                    trackfile.write(f"\t\tparent {currrep}{suffix}_{currstrand}tracks\n")
                    trackfile.write(f"\t\tshortLabel {currsample}{suffix} {currstrand} Strand\n")
                    trackfile.write(f"\t\tlongLabel Data from {self.expname} {currsample}{suffix} {currstrand} Strand\n")
                    trackfile.write(f"\t\tcolor {trackcolors[i % len(trackcolors)]}\n")
                    trackfile.write(f"\t\tbigDataUrl {currsample}{suffix}.{currstrand}.bw\n")
                    trackfile.write("\t\tvisibility full\n")
                    trackfile.write("\n")

    def createtrackdb(self, allreps: List[str]) -> None:
        with open(f"{self.expname}/trackhub/trackdb.txt", "w") as trackdb:
            currpriority = 2.3
            
            trackdb.write(f"track {self.expname}tracks            \n")
            trackdb.write("compositeTrack on                     \n")
            trackdb.write(f"shortLabel Data from {self.expname}   \n")
            trackdb.write(f"longLabel Data from {self.expname}    \n")
            trackdb.write("type bigWig                           \n")
            trackdb.write("dragAndDrop on                        \n")
            trackdb.write("autoScale on                          \n")
            trackdb.write("alwaysZero on                         \n")
            trackdb.write("maxHeightPixels 100:32:8              \n")
            trackdb.write("\n")
            
            for currrep in allreps:
                trackdb.write(f"track {currrep}plus                                             \n")
                trackdb.write(f"parent {self.expname}tracks                                   \n")
                trackdb.write(f"bigDataUrl {currrep}.Plus.bw                                 \n")
                trackdb.write(f"shortLabel Plus {currrep}                                    \n")
                trackdb.write(f"longLabel Plus strand coverage {currrep} all mapped reads    \n")
                trackdb.write("color 220,148,44                                             \n")
                trackdb.write("type bigWig                                                  \n")
                trackdb.write(f"priority {currpriority}                                       \n")
                trackdb.write("\n")
                
                trackdb.write(f"track {currrep}Minus                                             \n")
                trackdb.write(f"parent {self.expname}tracks                                    \n")
                trackdb.write(f"bigDataUrl {currrep}.Minus.bw                                 \n")
                trackdb.write(f"shortLabel Minus {currrep}                                    \n")
                trackdb.write(f"longLabel Minus strand coverage {currrep} all mapped reads    \n")
                trackdb.write("color 112,73,18                                       \n")
                trackdb.write("type bigWig                                                  \n")
                trackdb.write(f"priority {currpriority + 0.1}                                  \n")
                trackdb.write("\n\n\n")

                currpriority += 0.2

    def makebigwigs(self, bamfile: str, repname: str, faifile: str, directory: Union[str, Path], suffix: str = '', scalefactor: float = 1) -> None:
        self.logger.info(faifile)
        
        # Pipeline: samtools view -> genomeCoverageBed -> sort -> temp_bedgraph
        
        cmd_pipeline = f"samtools view -b -F 0x10 {bamfile} | genomeCoverageBed -scale {1./scalefactor} -bg -ibam stdin -g {faifile} | sort -k1,1 -k2,2n"
        
        # Plus strand
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmp_plus:
            tmp_plus_name = tmp_plus.name
            
        self.logger.info(f"Generating Plus strand BigWig for {repname}...")
        # Run pipeline to generate bedgraph
        full_cmd_plus = f"{cmd_pipeline} > {tmp_plus_name}"
        subprocess.run(full_cmd_plus, shell=True, check=True)
        
        # Run bedGraphToBigWig
        out_plus = Path(directory) / f"{repname}{suffix}.Plus.bw"
        subprocess.run(f"bedGraphToBigWig {tmp_plus_name} {faifile} {out_plus}", shell=True, check=True)
        os.remove(tmp_plus_name)

        # Minus strand (samtools view -f 0x10)
        cmd_pipeline_minus = f"samtools view -b -f 0x10 {bamfile} | genomeCoverageBed -scale {1./scalefactor} -bg -ibam stdin -g {faifile} | sort -k1,1 -k2,2n"
        
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmp_minus:
            tmp_minus_name = tmp_minus.name
            
        self.logger.info(f"Generating Minus strand BigWig for {repname}...")
        full_cmd_minus = f"{cmd_pipeline_minus} > {tmp_minus_name}"
        subprocess.run(full_cmd_minus, shell=True, check=True)
        
        out_minus = Path(directory) / f"{repname}{suffix}.Minus.bw"
        subprocess.run(f"bedGraphToBigWig {tmp_minus_name} {faifile} {out_minus}", shell=True, check=True)
        os.remove(tmp_minus_name)

    def maketracks(self, currbam: str, genomebam: str, currsample: str) -> str:
        # sizefactors is accessed from self.sizefactors
        scalefactor = self.sizefactors.get(currsample, 1)
        if scalefactor == 0:
             scalefactor = 1 # Avoid division by zero
             
        self.convertbam(currbam, genomebam, force=True)
        self.makebigwigs(genomebam, currsample, f"{self.dbname}-tRNAgenome.fa.fai", self.trackdir, scalefactor=scalefactor)
        return genomebam

    def run(self):
        if not self.trackdir.exists():
            self.trackdir.mkdir(parents=True)
            
        allsamples = self.sampledata.getsamples()
        
        # Index genome fasta
        subprocess.run(f"samtools faidx {self.dbname}-tRNAgenome.fa", shell=True, check=True)
        
        trackargs = []
        for currsample in allsamples:
            currbam = self.sampledata.getbam(currsample)
            genomebam = f"{currsample}-genome.bam"
            trackargs.append((currbam, genomebam, currsample))
            
        starttime = time.time()
        
        # Use multiprocessing
        if self.threads > 1:
            with Pool(processes=self.threads) as pool:
                results = pool.starmap(self.maketracks, trackargs)
                for i, res in enumerate(results):
                    self.logger.info(f"{res}:{time.time() - starttime}")
        else:
            for args in trackargs:
                res = self.maketracks(*args)
                self.logger.info(f"{res}:{time.time() - starttime}")

        with open(self.trackdir / "trackdb.txt", "w") as trackfile:
            self.createmultiwigtrackdb(trackfile, shortlabel="all", longlabel="all")

        # Merge replicates
        for currrep in self.sampledata.allreplicates():
            repsamples = self.sampledata.getrepsamples(currrep)
            merge_bam = f"{currrep}-mergegenome.bam"
            input_bams = [f"{curr}-genome.bam" for curr in repsamples]
            
            self.samtoolsmerge(input_bams, merge_bam, force=True)
            pysam.index(merge_bam)
            self.makebigwigs(merge_bam, currrep, f"{self.dbname}-tRNAgenome.fa.fai", self.trackdir)
