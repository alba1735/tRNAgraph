#!/usr/bin/env python3

import sys
import os
import shutil
import time
import subprocess
import multiprocessing
import logging
from collections import defaultdict
from typing import Tuple

# Custom external modules (Assuming these exist in the PYTHONPATH)
from . import toolsGetMaturetRNAs
from . import toolsAligntRNALocus
from . import toolsTG


def _covariance_models(orgmode: str) -> Tuple[str, str, bool]:
    '''
    Infernal covariance models for an organism mode: (mature model, tRNA model, prokaryotic mode).

    Prokaryotic mode is a property of the domain, not a separate switch -- bacteria and archaea
    lack the eukaryotic 3' CCA addition, so their models pair with prok_mode=True. `--forcecca`
    overrides that downstream; this function only reports the domain default.

    Model paths come from toolsTG.assets_dir() so they resolve inside the installed package under
    both `pip install -e .` and a plain `pip install .`.
    '''
    models = {
        'euk': ('trnamature-euk.cm', 'TRNAinf-euk.cm', False),
        'arch': ('trnamature-arch.cm', 'TRNAinf-arch.cm', True),
        'mito': ('TRNAMatureMitoinf.cm', 'TRNAinf.cm', False),
        'bact': ('trnamature-bact.cm', 'TRNAinf-bact.cm', True),
    }
    if orgmode not in models:
        # tRNADatabaseBuilder.__init__ already validates orgmode against POSITION_TABLES, so this
        # is unreachable in practice. It exists because the if/elif chain this replaced simply
        # left both model names unbound for an unrecognised mode, turning a typo into an
        # UnboundLocalError raised well away from its cause.
        raise ValueError(
            f"Unknown orgmode {orgmode!r} for covariance model selection. "
            f"Expected one of: {', '.join(sorted(models))}."
        )
    mature_name, trna_name, prok_mode = models[orgmode]
    cm_dir = os.path.join(toolsTG.assets_dir(), 'cm')
    return os.path.join(cm_dir, mature_name), os.path.join(cm_dir, trna_name), prok_mode


class tRNADatabaseBuilder:
    '''
    Class to generate a fasta file and indices containing mature tRNA sequences.
    '''
    def __init__(self, args):
        self.args = args
        self.logger = logging.getLogger(__name__)
        self.db_name = os.path.basename(args.output)
        
        # Handle directory pathing
        toolsTG.builder(args.output)
        self.db_directory = os.path.dirname(args.output)
        if self.db_directory and not self.db_directory.endswith('/'):
            self.db_directory += '/'
        elif self.db_directory == "":
            self.db_directory = "" # Keep empty if in current dir

        self.genome = args.genome
        self.trnaout = args.trnaout
        self.trnafa = args.trnafa
        self.namemap = args.namemap
        self.addtrna = args.addtrna
        self.addseqs = args.addseqs
        self.forcecca = args.forcecca
        self.orgmode = args.orgmode or "euk"
        from . import toolsGetCoverage
        if self.orgmode not in toolsGetCoverage.POSITION_TABLES:
            # An unrecognised mode used to fall through to eukaryotic positions,
            # so a typo produced a plausible-looking but wrongly numbered database.
            raise ValueError(
                f"Unknown organism mode {self.orgmode!r}. Expected one of: "
                + ", ".join(sorted(toolsGetCoverage.POSITION_TABLES))
            )
        self.threads = args.threads
        
        # Point to the package root (src/trnagraph) instead of modules/
        self.script_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + "/"
        
        # Initialize Sprinzl positions
        self._init_positions()

    def _init_positions(self):
        '''Sprinzl numbering positions per organism mode.

        Sourced from toolsGetCoverage.POSITION_TABLES rather than redefined here,
        so the database build and the coverage step can never disagree about what
        a position means.
        '''
        from . import toolsGetCoverage
        self.pos_maps = {
            mode: list(table)
            for mode, table in toolsGetCoverage.POSITION_TABLES.items()
        }

    def check_dependencies(self):
        '''Ensure required external tools are available'''
        for prog in ["samtools", "bowtie2-build"]:
            if shutil.which(prog) is None:
                self.logger.error(f"Error: Could not find '{prog}' in path.")
                sys.exit(1)

    def _run_shell(self, command, fail_quit=False):
        '''Wrapper for running shell commands'''
        try:
            result = subprocess.run(
                command, 
                shell=True, 
                check=True, 
                text=True, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT
            )
            if result.stdout:
                self.logger.info(result.stdout)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Command failed: {command}")
            if e.stdout:
                self.logger.error(e.stdout)
            if fail_quit:
                self.logger.error("Aborting program...")
                sys.exit(1)
            return e.returncode
        return 0

    def get_git_hash(self):
        '''Retrieve git version info matching trnagraph style'''
        version_str = "Unknown"
        hash_str = "Unknown"

        # Try to get version from package metadata
        try:
            from importlib.metadata import version
            version_str = version("tRNAgraph")
        except Exception:
            pass

        # Try to get git info
        try:
            # Assuming script_dir is src/trnagraph/, repo root is two levels up
            repo_root = os.path.dirname(os.path.dirname(self.script_dir.rstrip('/')))
            git_dir = os.path.join(repo_root, '.git')
            
            if os.path.exists(git_dir):
                # Get hash
                hash_str = subprocess.check_output(['git', '--git-dir='+git_dir, 'rev-parse', 'HEAD'], text=True).strip()
                
                # If version is still unknown, try git describe
                if version_str == "Unknown":
                     version_str = subprocess.check_output(['git', '--git-dir='+git_dir, 'describe'], text=True).strip()
        except Exception:
            pass
            
        return version_str, hash_str

    def get_trna_nums(self, trna_align, margin=0, mode="transcript"):
        '''
        Gets the tRNA numbers by the Sprinzl numbering system.
        '''
        trna_num = []
        curr_count = 0
        enum = 1
        gap_num = 1
        intron_num = 1
        
        positions = self.pos_maps.get(self.orgmode, self.pos_maps["euk"])

        for i in range(margin):
            trna_num.append('head' + str(margin - i))

        for i, struct in enumerate(trna_align.consensus):
            if curr_count >= len(positions):
                trna_num.append('gap' + str(gap_num))
                gap_num += 1
                curr_count += 1
            elif struct in set("+=*"):
                # Special case to account for differences between loci/transcripts
                if curr_count == 0 and struct == '=' and self.orgmode != "bact":
                    curr_count = 1
                    gap_num = 1
                
                if curr_count < len(positions) and str(positions[curr_count]).startswith('e'):
                    trna_num.append('e' + str(enum))
                    enum += 1
                    curr_count += 1
                    gap_num = 1
                elif curr_count < len(positions) and positions[curr_count] == '-':
                    trna_num.append(str(curr_count) + '.gap' + str(gap_num))
                    gap_num += 1
                    curr_count += 1
                else:
                    trna_num.append(str(positions[curr_count]))
                    curr_count += 1
                    gap_num = 1
            else:
                # Intron handling
                if curr_count < len(positions) and positions[curr_count] == 38:
                    trna_num.append('intron' + str(intron_num))
                    intron_num += 1
                else:
                    trna_num.append(str(curr_count) + '.gap' + str(gap_num))
                    gap_num += 1
        
        for i in range(margin):
            trna_num.append('tail' + str(i + 1))
            
        return trna_num

    def process_additional_seqs(self):
        '''Handle additional sequences if provided'''
        if not self.addseqs:
            return ""

        seq_fasta_name = self.db_directory + self.db_name + "-additionals.fa"
        seq_files = {}
        
        with open(seq_fasta_name, "w") as seq_fasta, \
             open(self.db_directory + self.db_name + "-otherseqs.txt", "w") as other_seqs:
            
            seq_counts = defaultdict(int)
            
            for curr_line in open(self.addseqs):
                fields = curr_line.split()
                if len(fields) < 2:
                    continue
                    
                seqs_name = fields[0]
                seqs_file = fields[1]
                
                print(f"{seqs_name}\t{self.db_name}-{seqs_name}_seq.fa\t{self.db_name}-{seqs_name}_seq.bed", file=other_seqs)
                seq_files[seqs_name] = seqs_file

                with open(self.db_directory + self.db_name + "-" + seqs_name + "_seq.bed", "w") as seq_bed:
                    for name, seq in toolsTG.read_multi_fasta(seqs_file):
                        display_name = name
                        if name not in seq_counts:
                            print(f">{name}", file=seq_fasta)
                        else:
                            display_name = f"{name}.{seq_counts[name]}"
                            print(f">{display_name}", file=seq_fasta)
                        
                        seq_counts[name] += 1
                        
                        # Add N buffers
                        padded_seq = (20 * "N") + seq.upper() + (20 * "N")
                        print(padded_seq, file=seq_fasta)
                        
                        # Write BED
                        print(f"{name}\t20\t{len(seq) + 20}\t{display_name}\t1000\t+", file=seq_bed)
                        
        return seq_fasta_name

    def main(self):
        runtime = time.time()
        loc_time = time.localtime(runtime)
        
        self.check_dependencies()

        # Index genome if needed
        if not os.path.isfile(self.genome + ".fai"):
            self._run_shell("samtools faidx " + self.genome)

        # Determine Covariance Models
        mature_model, trna_model, prok_mode = _covariance_models(self.orgmode)

        if self.forcecca:
            prok_mode = False

        # File paths
        transcript_stk = self.db_directory + self.db_name + "-trnaalign.stk"
        locus_stk = self.db_directory + self.db_name + "-trnaloci.stk"
        mature_trna_bed = self.db_directory + self.db_name + "-maturetRNAs.bed"
        trna_table = self.db_directory + self.db_name + "-trnatable.txt"
        loci_bed = self.db_directory + self.db_name + "-trnaloci.bed"
        mature_trna_fa = self.db_directory + self.db_name + "-maturetRNAs.fa"

        # Run external modules
        toolsGetMaturetRNAs.get_mature_trnas(
            trnascan=[self.trnaout], 
            genome=self.genome,
            gtrnafa=self.trnafa,
            addtrna=self.addtrna, 
            namemap=self.namemap, 
            bedfile=mature_trna_bed,
            maturetrnatable=trna_table,
            trnaalignment=transcript_stk,
            locibed=loci_bed,
            maturetrnafa=mature_trna_fa,
            cmmodel=mature_model, 
            prokmode=prok_mode
        )

        toolsAligntRNALocus.align_trna_locus(
            genomefile=self.genome,
            stkfile=locus_stk,
            trnaloci=loci_bed, 
            cmmodel=trna_model
        )

        # Process Alignments
        # Assuming toolsTG.read_rna_stk returns an iterator/list
        with open(transcript_stk, "r") as f:
            trans_align = list(toolsTG.read_rna_stk(f))[0]
        
        trans_nums = self.get_trna_nums(trans_align, margin=0, mode="transcript")
        
        with open(self.db_directory + self.db_name + "-alignnum.txt", "w") as align_num_file:
            for i, _ in enumerate(trans_align.consensus):
                print(f"{i}\t{trans_align.consensus[i]}\t{trans_align.structure[i]}\t{trans_nums[i]}", file=align_num_file)

        with open(locus_stk, "r") as f:
            locial_ign = list(toolsTG.read_rna_stk(f))[0]
            
        loci_nums = self.get_trna_nums(locial_ign, margin=0, mode="locus")
        
        with open(self.db_directory + self.db_name + "-locusnum.txt", "w") as loci_num_file:
            for i, _ in enumerate(locial_ign.consensus):
                print(f"{i}\t{locial_ign.consensus[i]}\t{locial_ign.structure[i]}\t{loci_nums[i]}", file=loci_num_file)

        # Process additional sequences
        seq_fasta_name = self.process_additional_seqs()
        
        # Concatenate genome files
        final_fasta = self.db_directory + self.db_name + "-tRNAgenome.fa"
        files_to_cat = [mature_trna_fa, self.genome]
        if seq_fasta_name:
            files_to_cat.append(seq_fasta_name)
            
        self._run_shell(f"cat {' '.join(files_to_cat)} > {final_fasta}", fail_quit=True)

        # Generate DB Info
        git_version, _ = self.get_git_hash()
        
        with open(self.db_directory + self.db_name + "-dbinfo.txt", "w") as db_info:
            print(f"time\t{runtime}({loc_time.tm_mon}/{loc_time.tm_mday}/{loc_time.tm_year})", file=db_info)
            print(f"creation\t{' '.join(sys.argv)}", file=db_info)
            print(f"genome\t{self.genome}", file=db_info)
            print(f"trnascanout\t{self.trnaout}", file=db_info)
            print(f"orgmode\t{self.orgmode}", file=db_info)
            print(f"forcecca\t{self.forcecca}", file=db_info)
            print(f"git version\t{git_version}", file=db_info)
            if self.addseqs:
                # Note: Re-parsing logic for addseqs names might be needed here if specific logging is required,
                # currently simplified to match original intent
                print(f"additional transcripts file\t{self.addseqs}", file=db_info)

        # Build Bowtie Index
        index_option = "--large-index"
        index_base = self.db_directory + self.db_name + "-tRNAgenome"

        self.logger.info(f"Building Bowtie2 index using {self.threads} threads...")
        self._run_shell(
            f"bowtie2-build {final_fasta} {index_base} {index_option} -p {self.threads}",
            fail_quit=True
        )
        self.logger.info("Database creation complete.")


if __name__ == '__main__':
    pass