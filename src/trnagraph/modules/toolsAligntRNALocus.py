#!/usr/bin/env python3

import sys
import subprocess
import logging
from pathlib import Path
from typing import Union

from . import toolsTG

logger = logging.getLogger(__name__)

def align_trna_locus(
    stkfile: Union[str, Path],
    genomefile: Union[str, Path],
    trnaloci: Union[str, Path],
    cmmodel: Union[str, Path]
) -> None:
    """
    Aligns tRNA loci sequences to a covariance model using cmalign.

    Args:
        stkfile: Output Stockholm file path.
        genomefile: Input genome FASTA file path.
        trnaloci: Input BED file path containing tRNA loci.
        cmmodel: Covariance model file path.
    """
    genome_path = Path(genomefile).expanduser()
    
    # Read bed file and get sequences
    # toolsTG.readbed expects strings for file paths
    trnaloci_list = list(toolsTG.readbed(str(trnaloci), orgdb="genome", seqfile=str(genome_path)))
    
    # toolsTG.getseqdict expects faifiles dict
    lociseqs = toolsTG.getseqdict(trnaloci_list, faifiles={"genome": f"{genome_path}.fai"})
    
    # Create temporary fasta file and run cmalign
    # Using context manager ensures the temporary file is properly closed/deleted
    with toolsTG.tempmultifasta(lociseqs.items()) as seqfile:
        cmcommand = [
            'cmalign', 
            "-o", str(stkfile), 
            "--nonbanded", 
            "--notrunc", 
            "-g", 
            str(cmmodel), 
            seqfile.name
        ]
        
        try:
            subprocess.run(
                cmcommand, 
                check=True, 
                stdout=subprocess.DEVNULL, 
                stderr=subprocess.PIPE,
                text=True
            )
        except subprocess.CalledProcessError as e:
            logger.error(f"Command failed: {' '.join(cmcommand)}")
            logger.error(f"Error: {e.stderr}")
            logger.error("Failure to align tRNAs")
            sys.exit(1)

if __name__ == "__main__":
    pass