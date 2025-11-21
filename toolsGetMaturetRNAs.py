#!/usr/bin/env python3

import sys
import re
import subprocess
from collections import defaultdict
from contextlib import ExitStack
from pathlib import Path
from typing import List, Optional, Dict

import toolsTG

def get_mature_trnas(
    trnascan: Optional[List[str]] = None,
    genome: Optional[str] = None,
    maturetrnafa: Optional[str] = None,
    bedfile: Optional[str] = None,
    maturetrnatable: Optional[str] = None,
    locibed: Optional[str] = None,
    trnaalignment: Optional[str] = None,
    gtrnafa: Optional[str] = None,
    cmmodel: Optional[str] = None,
    prokmode: bool = False,
    namemap: Optional[str] = None,
    addtrna: Optional[str] = None
) -> None:
    """
    Process and generate mature tRNA sequences and annotations.

    Args:
        trnascan: List of tRNAscan-SE output files.
        genome: Path to the genome FASTA file.
        maturetrnafa: Output path for mature tRNA FASTA.
        bedfile: Output path for tRNA BED file.
        maturetrnatable: Output path for tRNA table.
        locibed: Output path for loci BED file.
        trnaalignment: Output path for tRNA alignment (Stockholm).
        gtrnafa: Path to GtRNAdb FASTA file for name mapping.
        cmmodel: Path to covariance model for alignment.
        prokmode: Boolean flag for prokaryotic mode (affects CCA addition).
        namemap: Path to name mapping file.
        addtrna: Path to additional tRNA FASTA file.
    """
    
    trnascan_files = trnascan or []
    
    if not genome:
        raise ValueError("Genome file must be specified.")
        
    genome_path = Path(genome).expanduser()
    genome_str = str(genome_path)

    with ExitStack() as stack:
        maturetrnafa_file = stack.enter_context(open(maturetrnafa, "w")) if maturetrnafa else sys.stdout
        trnabed_file = stack.enter_context(open(bedfile, "w")) if bedfile else None
        trnatable_file = stack.enter_context(open(maturetrnatable, "w")) if maturetrnatable else None
        locibed_file = stack.enter_context(open(locibed, "w")) if locibed else None

        gtrnatrans: Optional[Dict[str, str]] = None
        
        if namemap:
            gtrnatrans = {}
            with open(namemap, "r") as f:
                for currline in f:
                    fields = currline.rstrip().split()
                    if "tRNAscan-SE_id" == fields[0] or len(fields) < 2 or len(fields[1].split("-")) < 4:
                        continue
                    shortname = fields[0].split("-")[0]
                    gtrnatrans[shortname] = fields[1]
        elif gtrnafa:
            trnanamere = re.compile(r"^\>\S+_((?:\w+\-)?tRNA\-\w+\-[\w\?]+\-\d+\-\d+)\s+\((?:tRNAscan\-SE\s+ID:\s+)?(\S+)\)")
            gtrnatrans = {}
            with open(gtrnafa, "r") as f:
                for currline in f:
                    trnamatch = trnanamere.match(currline)
                    if trnamatch:
                        gtrnatrans[trnamatch.group(2)] = trnamatch.group(1)
            
            if not gtrnatrans:
                print("Could not extract names from gtrnadb fasta file", file=sys.stderr)
                sys.exit(1)

        alltrnas: List[toolsTG.tRNAtranscript] = []
        trnascantrnas: List[toolsTG.tRNAlocus] = []
        trnadbtrnas: List[toolsTG.tRNAtranscript] = []
        
        for currfile in trnascan_files:
            if gtrnatrans:
                trnadbtrnas.extend(toolsTG.readtRNAdb(currfile, genome_str, gtrnatrans))
            else:
                trnascantrnas.extend(toolsTG.readtRNAscan(currfile, genome_str))
                
        extratrnas: List[toolsTG.tRNAtranscript] = []
        if addtrna:
            for currname, currseq in toolsTG.read_multi_fasta(addtrna):
                namefields = currname.split("-")
                if len(namefields) != 4:
                    print(f"additional tRNA {currname} from {addtrna} does not use a valid tRNA name", file=sys.stderr)
                    sys.exit(1)
                
                curramino = namefields[1]
                currac = namefields[2]
                if curramino == "Und":
                    curramino = "Undet"
                extratrnas.append(toolsTG.tRNAtranscript(
                    currseq.replace("U", "T"), None, curramino, currac, [], None, name=currname, artificialtrna=True
                ))
                
        alltrnas = list(toolsTG.getuniquetRNAs(trnascantrnas)) + trnadbtrnas + extratrnas
        
        margin = 20
        anticodoncount: Dict[str, int] = defaultdict(int)
        
        for currtrans in alltrnas:
            if currtrans.name is None:
                name = f'tRNA-{currtrans.amino}{currtrans.anticodon}{anticodoncount[currtrans.anticodon] + 1}'
                currtrans.name = name
                anticodoncount[currtrans.anticodon] += 1
                
        if not alltrnas:
            print("No trna sequences", file=sys.stderr)
            sys.exit(1)
            
        if trnaalignment and cmmodel:
            # Use toolsTG.tempmultifasta with context manager
            with toolsTG.tempmultifasta(((currtrans.name, currtrans.getmatureseq(addcca=not prokmode)) for currtrans in alltrnas)) as seqfile:
                cmcommand = ['cmalign', "-o", trnaalignment, "--nonbanded", "--notrunc", "-g", cmmodel, seqfile.name]
                
                try:
                    subprocess.run(cmcommand, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
                except subprocess.CalledProcessError as e:
                    print(f"Command failed: {' '.join(cmcommand)}", file=sys.stderr)
                    print(f"Error: {e.stderr}", file=sys.stderr)
                    print("Failure to align tRNAs", file=sys.stderr)
                    sys.exit(1)

        locibed_lines = []

        for currtrans in alltrnas:
            name = currtrans.name
            # Write to maturetrnafa
            print(f">{name}", file=maturetrnafa_file)
            matureseq = currtrans.getmatureseq(addcca=not prokmode)
            print(f"{'N' * margin}{matureseq}{'N' * margin}", file=maturetrnafa_file)

            locinames = ",".join(currlocus.name for currlocus in sorted(currtrans.loci, key=lambda x: x.name)) if currtrans.loci else "NA"
            if not locinames:
                locinames = "NA"

            if trnatable_file:
                print(f"{name}\t{locinames}\t{currtrans.amino}\t{currtrans.anticodon}", file=trnatable_file)
            
            if trnabed_file:
                transcriptrange = toolsTG.GenomeRange("genome", name, margin, margin + len(matureseq), strand="+", name=name)
                print(transcriptrange.bedstring(), file=trnabed_file)
            
            if locibed_file:
                itemrgb = "0"
                for currlocus in currtrans.loci:
                    trnalength = currlocus.loc.end - currlocus.loc.start
                    if currlocus.intron is None:
                        line = f"{currlocus.loc.bedstring()}\t{currlocus.loc.start}\t{currlocus.loc.end}\t{itemrgb}\t1\t{trnalength}\t0"
                    else:
                        blockcounts = 2
                        blocksizes = f"{currlocus.intron[0] + 1},{trnalength - currlocus.intron[1] - 1}"
                        blockstarts = f"0,{currlocus.intron[1] + 1}"
                        line = f"{currlocus.loc.bedstring()}\t{currlocus.loc.start}\t{currlocus.loc.end}\t{itemrgb}\t{blockcounts}\t{blocksizes}\t{blockstarts}"
                    locibed_lines.append((currlocus.name, line))

        if locibed_file:
            # Sort by tRNA name
            locibed_lines.sort(key=lambda x: x[0])
            for _, line in locibed_lines:
                print(line, file=locibed_file)

if __name__ == '__main__':
    pass