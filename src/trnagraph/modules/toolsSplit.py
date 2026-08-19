#!/usr/bin/env python3

import os
import sys
import shutil
import subprocess
import multiprocessing
import re
import itertools
import pysam

class BamSplitter:
    '''
    Splits BAM files by read length. Internal helper used by `analyze build
    --readlengthsplit`/`analyze addsplit` (adataBuild.py) -- not exposed as its own CLI
    command. Callers treat the `u<N>`/`o<N>` output directories this writes under `bamdir`
    as scratch: by default they're deleted once merged into the AnnData object, unless
    `--savesplitbams` is passed.
    '''
    def __init__(self, args):
        self.args = args
        self.manifest = os.path.abspath(args.input)
        self.cutoff = args.readlengthsplit
        # If bamdir is not provided, assume current directory
        self.bamdir = os.path.abspath(args.bamdir) if args.bamdir else os.getcwd()
        self.overwrite = getattr(args, 'overwritebams', False)
        
        # Validate inputs
        if not os.path.isfile(self.manifest):
            print(f"Error: Manifest file '{self.manifest}' not found.")
            sys.exit(1)
            
        if not os.path.isdir(self.bamdir):
            print(f"Error: BAM directory '{self.bamdir}' not found.")
            sys.exit(1)
            
        # Check for samtools
        if shutil.which('samtools') is None:
            print("Error: samtools not found in PATH. Please install samtools.")
            sys.exit(1)

        self.samples = self._parse_manifest()

    def _parse_manifest(self):
        '''
        Parses tab-delimited manifest file to get sample names.
        Assumes first column is Sample ID.
        '''
        samples = []
        try:
            with open(self.manifest, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#') or line.lower().startswith('fastq'):
                        continue
                    parts = line.split()
                    # New format: fastq, sample, group. Sample is 2nd column.
                    if len(parts) >= 2:
                        samples.append(parts[1])
                        
        except Exception as e:
            print(f"Error parsing manifest: {e}")
            sys.exit(1)
        return samples

    def _split_bam(self, sample):
        '''
        Splits a single BAM file.
        '''
        bam_filename = f"{sample}.bam"
        bam_path = os.path.join(self.bamdir, bam_filename)
        
        if not os.path.isfile(bam_path):
            print(f"Warning: File {bam_path} not found. Skipping.")
            return

        # Define output directories
        dir_under = os.path.join(self.bamdir, f"u{self.cutoff}")
        dir_over = os.path.join(self.bamdir, f"o{self.cutoff}")
        
        os.makedirs(dir_under, exist_ok=True)
        os.makedirs(dir_over, exist_ok=True)
        
        out_under = os.path.join(dir_under, f"{sample}.bam")
        out_over = os.path.join(dir_over, f"{sample}.bam")
        
        # Check if files exist and overwrite is False
        if not self.overwrite and os.path.isfile(out_under) and os.path.isfile(out_over):
            print(f"Skipping {sample} (split files exist)")
            return

        print(f"Splitting {bam_filename} (Cutoff: {self.cutoff})...")
        
        # Commands
        # Reads UNDER (<) the cutoff
        cmd_under = [
            "samtools", "view", "-h", "-b", 
            "-e", f"length(seq) < {self.cutoff}", 
            "-o", out_under, bam_path
        ]
        
        # Reads OVER OR EQUAL (>=) to the cutoff
        cmd_over = [
            "samtools", "view", "-h", "-b", 
            "-e", f"length(seq) >= {self.cutoff}", 
            "-o", out_over, bam_path
        ]
        
        # Run in parallel (using subprocess.Popen)
        try:
            p1 = subprocess.Popen(cmd_under)
            p2 = subprocess.Popen(cmd_over)
            
            p1.wait()
            p2.wait()
            
            if p1.returncode != 0:
                print(f"Error splitting (under) for {sample}")
            if p2.returncode != 0:
                print(f"Error splitting (over) for {sample}")
            else:
                print(f"  -> Created {out_under}")
                print(f"  -> Created {out_over}")
                
        except Exception as e:
            print(f"Error executing samtools for {sample}: {e}")

    def _write_mapinfo_files(self, output_dir, basename):
        '''
        Writes mapinfo.txt and trnamapinfo.txt with zero values.
        These files are created during split since bowtie2 is not run, 
        but the files are needed for downstream pipeline steps.
        '''
        results_dir = os.path.join(output_dir, "results")
        graphs_dir = os.path.join(output_dir, "graphs")
        
        os.makedirs(results_dir, exist_ok=True)
        os.makedirs(graphs_dir, exist_ok=True)
        
        # Write mapinfo.txt with zeros
        mapinfo_path = os.path.join(results_dir, f"{basename}-mapinfo.txt")
        with open(mapinfo_path, 'w') as f:
            # Header row with sample names
            f.write("\t".join(self.samples) + "\n")
            # Data rows with zeros
            f.write("unmap\t" + "\t".join(["0"] * len(self.samples)) + "\n")
            f.write("single\t" + "\t".join(["0"] * len(self.samples)) + "\n")
            f.write("multi\t" + "\t".join(["0"] * len(self.samples)) + "\n")
        print(f"  -> Created {mapinfo_path}")
        
        # Write trnamapinfo.txt with zeros
        trnamapinfo_path = os.path.join(results_dir, f"{basename}-trnamapinfo.txt")
        with open(trnamapinfo_path, 'w') as f:
            # Header row with sample names
            f.write("\t".join(self.samples) + "\n")
            # Data rows with zeros
            f.write("multtrans\t" + "\t".join(["0"] * len(self.samples)) + "\n")
            f.write("multac\t" + "\t".join(["0"] * len(self.samples)) + "\n")
            f.write("multamino\t" + "\t".join(["0"] * len(self.samples)) + "\n")
            f.write("trna\t" + "\t".join(["0"] * len(self.samples)) + "\n")
            f.write("singlenon\t" + "\t".join(["0"] * len(self.samples)) + "\n")
            f.write("multiplenon\t" + "\t".join(["0"] * len(self.samples)) + "\n")
            f.write("total\t" + "\t".join(["0"] * len(self.samples)) + "\n")
        print(f"  -> Created {trnamapinfo_path}")

    def process(self):
        print(f"Processing BAM files in {self.bamdir}...")
        
        # If user provided threads, we can use them.
        threads = self.args.threads if hasattr(self.args, 'threads') and self.args.threads > 0 else 1
        
        if threads > 1:
            with multiprocessing.Pool(threads) as pool:
                pool.map(self._split_bam, self.samples)
        else:
            for sample in self.samples:
                self._split_bam(sample)
        
        parent_dir = os.path.dirname(self.bamdir)  
        exp_basename = os.path.basename(parent_dir)  
        exp_dir = os.path.dirname(parent_dir)  
        
        base_output_dir = exp_basename
        
        under_output_dir = base_output_dir
        over_output_dir = base_output_dir
        
        under_results_name = f"results_u{self.cutoff}"
        over_results_name = f"results_o{self.cutoff}"
        under_graphs_name = f"graphs_u{self.cutoff}"
        over_graphs_name = f"graphs_o{self.cutoff}"
        
        print(f"Creating output directories for u{self.cutoff}...")
        under_exp_dir = base_output_dir
        os.makedirs(os.path.join(under_exp_dir, under_results_name), exist_ok=True)
        os.makedirs(os.path.join(under_exp_dir, under_graphs_name), exist_ok=True)
        
        self._write_mapinfo_for_split(under_exp_dir, under_results_name, exp_basename, f"u{self.cutoff}")
        
        print(f"Creating output directories for o{self.cutoff}...")
        os.makedirs(os.path.join(over_output_dir, over_results_name), exist_ok=True)
        os.makedirs(os.path.join(over_output_dir, over_graphs_name), exist_ok=True)
        
        self._write_mapinfo_for_split(over_output_dir, over_results_name, exp_basename, f"o{self.cutoff}")
                
        print("Done.")

    def _analyze_bam(self, bam_path):
        '''
        Analyze a BAM file to extract mapping statistics.
        Returns a dict with mapinfo and trnamapinfo values.
        
        The BAM files have already been processed by toolsMap.process_mappings(),
        which adds tags like YR (multi-tRNA), YA (multi-anticodon), YM (multi-amino)
        to tRNA reads. Non-tRNA reads don't have these tags.
        '''
        stats = {
            # mapinfo values
            'unmap': 0,
            'single': 0,
            'multi': 0,
            # trnamapinfo values
            'multtrans': 0,  # YR > 1
            'multac': 0,     # YA > 1
            'multamino': 0,  # YM > 1
            'trna': 0,       # All tRNA reads
            'singlenon': 0,  # Unique non-tRNA
            'multiplenon': 0, # Multi-mapped non-tRNA
        }
        
        if not os.path.isfile(bam_path):
            print(f"Warning: BAM file not found: {bam_path}")
            return stats
        
        try:
            bamfile = pysam.AlignmentFile(bam_path, "rb")
        except Exception as e:
            print(f"Warning: Could not open BAM file {bam_path}: {e}")
            return stats
        
        # Group reads by query name to count unique reads (not alignments)
        seen_reads = set()
        
        for read in bamfile.fetch(until_eof=True):
            # Skip if we've already counted this read
            if read.query_name in seen_reads:
                continue
            seen_reads.add(read.query_name)
            
            # Check if unmapped
            if read.is_unmapped:
                stats['unmap'] += 1
                continue
            
            # Check for tRNA tags (added by process_mappings)
            try:
                yr = read.get_tag('YR')  # Number of tRNA transcripts
                ya = read.get_tag('YA')  # Number of anticodons
                ym = read.get_tag('YM')  # Number of aminos
                
                # This is a tRNA read
                stats['trna'] += 1
                
                if yr > 1:
                    stats['multtrans'] += 1
                    stats['multi'] += 1
                else:
                    stats['single'] += 1
                    
                if ya > 1:
                    stats['multac'] += 1
                if ym > 1:
                    stats['multamino'] += 1
                    
            except KeyError:
                # No tRNA tags - this is a non-tRNA read
                # Check if it's a multi-mapped read by checking secondary flag or NH tag
                try:
                    nh = read.get_tag('NH')
                    if nh > 1:
                        stats['multiplenon'] += 1
                        stats['multi'] += 1
                    else:
                        stats['singlenon'] += 1
                        stats['single'] += 1
                except KeyError:
                    # No NH tag, check if secondary
                    if read.is_secondary or read.is_supplementary:
                        stats['multiplenon'] += 1
                        stats['multi'] += 1
                    else:
                        stats['singlenon'] += 1
                        stats['single'] += 1
        
        bamfile.close()
        
        # Calculate total
        stats['total'] = stats['trna'] + stats['singlenon'] + stats['multiplenon']
        
        return stats

    def _write_mapinfo_for_split(self, exp_dir, results_dir_name, exp_basename, suffix):
        '''
        Writes mapinfo.txt and trnamapinfo.txt with calculated values from BAM files.
        '''
        results_dir = os.path.join(exp_dir, results_dir_name)
        
        # Determine which BAM directory to use based on suffix
        split_bamdir = os.path.join(self.bamdir, suffix)
        
        # Collect stats for all samples
        print(f"  Analyzing BAM files in {split_bamdir}...")
        all_stats = {}
        for sample in self.samples:
            bam_path = os.path.join(split_bamdir, f"{sample}.bam")
            all_stats[sample] = self._analyze_bam(bam_path)
        
        # Write mapinfo.txt
        mapinfo_path = os.path.join(results_dir, f"{exp_basename}-mapinfo.txt")
        with open(mapinfo_path, 'w') as f:
            f.write("\t".join(self.samples) + "\n")
            f.write("unmap\t" + "\t".join([str(all_stats[s]['unmap']) for s in self.samples]) + "\n")
            f.write("single\t" + "\t".join([str(all_stats[s]['single']) for s in self.samples]) + "\n")
            f.write("multi\t" + "\t".join([str(all_stats[s]['multi']) for s in self.samples]) + "\n")
        print(f"  -> Created {mapinfo_path}")
        
        # Write trnamapinfo.txt
        trnamapinfo_path = os.path.join(results_dir, f"{exp_basename}-trnamapinfo.txt")
        with open(trnamapinfo_path, 'w') as f:
            f.write("\t".join(self.samples) + "\n")
            f.write("multtrans\t" + "\t".join([str(all_stats[s]['multtrans']) for s in self.samples]) + "\n")
            f.write("multac\t" + "\t".join([str(all_stats[s]['multac']) for s in self.samples]) + "\n")
            f.write("multamino\t" + "\t".join([str(all_stats[s]['multamino']) for s in self.samples]) + "\n")
            f.write("trna\t" + "\t".join([str(all_stats[s]['trna']) for s in self.samples]) + "\n")
            f.write("singlenon\t" + "\t".join([str(all_stats[s]['singlenon']) for s in self.samples]) + "\n")
            f.write("multiplenon\t" + "\t".join([str(all_stats[s]['multiplenon']) for s in self.samples]) + "\n")
            f.write("total\t" + "\t".join([str(all_stats[s]['total']) for s in self.samples]) + "\n")
        print(f"  -> Created {trnamapinfo_path}")
