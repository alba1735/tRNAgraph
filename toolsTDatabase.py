#!/usr/bin/env python3

import sys
import os
import shutil
import time
import subprocess
import multiprocessing
from collections import defaultdict

# Custom external modules (Assuming these exist in the PYTHONPATH)
import toolsGetMaturetRNAs
import toolsAligntRNALocus
import toolsTG

class tRNADatabaseBuilder:
    '''
    Class to generate a fasta file and indices containing mature tRNA sequences.
    '''
    def __init__(self, args):
        self.args = args
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
        self.threads = args.threads
        
        self.script_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
        
        # Initialize Sprinzl positions
        self._init_positions()

    def _init_positions(self):
        '''Define Sprinzl numbering positions for different organism types'''
        common_start = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45]
        common_end = [46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76]
        e_vars = ['e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14']
        
        self.pos_maps = {
            "euk": [-1] + common_start[:17] + ['17a'] + common_start[17:20] + ['20a','20b'] + common_start[20:] + e_vars + ['e15','e16','e17','e18','e19'] + common_end,
            "arch": [-1] + common_start[:17] + ['17a'] + common_start[17:20] + ['20a','20b'] + common_start[20:] + e_vars + common_end,
            "bact": common_start + e_vars + ['e15','e16','e17'] + common_end,
            "mito": [-1] + common_start + e_vars + ['e15','e16','e17'] + common_end
        }

    def check_dependencies(self):
        '''Ensure required external tools are available'''
        for prog in ["samtools", "bowtie2-build"]:
            if shutil.which(prog) is None:
                print(f"Error: Could not find '{prog}' in path.", file=sys.stderr)
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
                print(result.stdout, end='')
        except subprocess.CalledProcessError as e:
            print(f"Command failed: {command}", file=sys.stderr)
            if e.stdout:
                print(e.stdout, file=sys.stderr, end='')
            if fail_quit:
                print("Aborting program...", file=sys.stderr)
                sys.exit(1)
            return e.returncode
        return 0

    def get_git_hash(self):
        '''Retrieve git version info matching trnagraph style'''
        try:
            git_version = subprocess.check_output(['git', '--git-dir='+self.script_dir+'.git', 'describe'], text=True).strip()
            git_hash = subprocess.check_output(['git', '--git-dir='+self.script_dir+'.git', 'rev-parse', 'HEAD'], text=True).strip()
            return git_version, git_hash
        except Exception:
            return "Unknown", "Unknown"

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
        prok_mode = False
        if self.orgmode == "euk":
            mature_model = self.script_dir + 'assets/cm/trnamature-euk.cm'
            trna_model = self.script_dir + 'assets/cm/TRNAinf-euk.cm'
        elif self.orgmode == "arch":
            mature_model = self.script_dir + 'assets/cm/trnamature-arch.cm'
            trna_model = self.script_dir + 'assets/cm/TRNAinf-arch.cm'
            prok_mode = True
        elif self.orgmode == "mito":
            mature_model = self.script_dir + 'assets/cm/TRNAMatureMitoinf.cm'
            trna_model = self.script_dir + 'assets/cm/TRNAinf.cm'
        elif self.orgmode == "bact":
            mature_model = self.script_dir + 'assets/cm/trnamature-bact.cm'
            trna_model = self.script_dir + 'assets/cm/TRNAinf-bact.cm'
            prok_mode = True
            
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

        print(f"Building Bowtie2 index using {self.threads} threads...")
        self._run_shell(
            f"bowtie2-build {final_fasta} {index_base} {index_option} -p {self.threads}",
            fail_quit=True
        )
        print("Database creation complete.")


if __name__ == '__main__':
    pass