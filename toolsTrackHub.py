#!/usr/bin/env python3

import sys
import os
import subprocess
import tempfile
import time
import pysam
from pathlib import Path
from typing import List, Optional, Union, Any
from multiprocessing import Pool
from shutil import which

# Try to import from local toolsTG, otherwise assume it's in the path
try:
    from .toolsTG import samplefile, getsizefactors
except ImportError:
    try:
        from toolsTG import samplefile, getsizefactors
    except ImportError:
        # Fallback for when running from tRNAgraph directory
        sys.path.append(os.path.dirname(__file__))
        from toolsTG import samplefile, getsizefactors

class TrackHubBuilder:
    def __init__(self, 
                 genomedatabase: str, 
                 samplefilename: str, 
                 expname: str, 
                 scriptdir: Optional[str] = None,
                 threads: int = 8):
        self.dbname = genomedatabase
        self.samplefilename = samplefilename
        self.expname = expname
        self.threads = threads
        
        # Determine script directory (tRAX directory)
        if scriptdir:
            self.scriptdir = Path(scriptdir)
        else:
            # Assume tRAX is sibling to tRNAgraph or in parent
            # Current file is in tRNAgraph/toolsTrackHub.py
            # tRAX scripts are in ../tRAX/
            current_file = Path(__file__).resolve()
            self.scriptdir = current_file.parent.parent / 'tRAX'
            
        if not self.scriptdir.exists():
             # Fallback to current directory or check if we are in tRAX
             if (Path.cwd() / 'convertbam.py').exists():
                 self.scriptdir = Path.cwd()
             else:
                 print(f"Warning: Could not find tRAX script directory at {self.scriptdir}", file=sys.stderr)

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
            print(f"Could not find {program} in path", file=sys.stderr)
            print("Aborting", file=sys.stderr)
            sys.exit(1)
        return progloc

    def convertbam(self, inputbam: str, outputbam: str, force: bool = False, logfile: Optional[Any] = sys.stderr) -> None:
        if not os.path.isfile(outputbam) or force:
            tempprefix = f"{Path(inputbam).stem}_{os.getpid()}"
            temp_dir = tempfile.gettempdir()
            
            convert_script = self.scriptdir / 'convertbam.py'
            
            cmd = f"{convert_script} {inputbam} {self.dbname} | samtools sort -T {temp_dir}/convert{tempprefix} - -o {outputbam}"
            
            print(cmd, file=sys.stderr)
            
            if logfile:
                logfile.flush()
                
            process = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, universal_newlines=True)
            _, errinfo = process.communicate()
            
            if logfile and errinfo:
                print(errinfo, file=logfile)
                logfile.flush()
                
            if process.returncode:
                print("Failure to convert bam to genome space", file=sys.stderr)
                print("check logfile", file=sys.stderr)
                sys.exit(1)

    def samtoolsmerge(self, bamfiles: List[str], outbam: str, force: bool = False) -> None:
        samtoolsloc = self.get_location("samtools")
        if not os.path.isfile(outbam) or force:
            samcommand = [samtoolsloc, "merge", "-f", outbam] + bamfiles
            
            print(" ".join(samcommand), file=sys.stderr)
            result = subprocess.run(samcommand, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, check=False)
            if result.stdout:
                print(result.stdout, file=sys.stderr)

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

    def makebigwigs(self, bamfile: str, repname: str, faifile: str, directory: Union[str, Path], filterloci: bool = False, suffix: str = '', scalefactor: float = 1) -> None:
        filtercommand = ''
        if filterloci:
            filter_script = self.scriptdir / 'filterunique.py'
            filtercommand = f"{filter_script} --uniqloci | "
        
        print(faifile, file=sys.stderr)
        
        # Pipeline: samtools view -> (filter) -> genomeCoverageBed -> sort -> temp_bedgraph
        
        cmd_pipeline = f"samtools view -b -F 0x10 {bamfile} | {filtercommand} genomeCoverageBed -scale {1./scalefactor} -bg -ibam stdin -g {faifile} | sort -k1,1 -k2,2n"
        
        # Plus strand
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmp_plus:
            tmp_plus_name = tmp_plus.name
            
        print(f"Generating Plus strand BigWig for {repname}...", file=sys.stderr)
        # Run pipeline to generate bedgraph
        full_cmd_plus = f"{cmd_pipeline} > {tmp_plus_name}"
        subprocess.run(full_cmd_plus, shell=True, check=True)
        
        # Run bedGraphToBigWig
        out_plus = Path(directory) / f"{repname}{suffix}.Plus.bw"
        subprocess.run(f"bedGraphToBigWig {tmp_plus_name} {faifile} {out_plus}", shell=True, check=True)
        os.remove(tmp_plus_name)

        # Minus strand (samtools view -f 0x10)
        cmd_pipeline_minus = f"samtools view -b -f 0x10 {bamfile} | {filtercommand} genomeCoverageBed -scale {1./scalefactor} -bg -ibam stdin -g {faifile} | sort -k1,1 -k2,2n"
        
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmp_minus:
            tmp_minus_name = tmp_minus.name
            
        print(f"Generating Minus strand BigWig for {repname}...", file=sys.stderr)
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
                    print(f"{res}:{time.time() - starttime}", file=sys.stderr)
        else:
            for args in trackargs:
                res = self.maketracks(*args)
                print(f"{res}:{time.time() - starttime}", file=sys.stderr)

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
