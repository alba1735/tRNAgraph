#!/usr/bin/env python3

import os
import sys
import shutil
import subprocess
import multiprocessing
import re

class BamSplitter:
    '''
    Class to split BAM files based on read length.
    '''
    def __init__(self, args):
        self.args = args
        self.manifest = os.path.abspath(args.input)
        self.cutoff = args.cutoff
        # If bamdir is not provided, assume current directory
        self.bamdir = os.path.abspath(args.bamdir) if args.bamdir else os.getcwd()
        
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
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split()
                    if parts:
                        samples.append(parts[0])
        except Exception as e:
            print(f"Error parsing manifest: {e}")
            sys.exit(1)
        return samples

    def _create_manifests(self):
        '''
        Generates new manifest files with suffixes.
        '''
        base_dir = os.path.dirname(self.manifest)
        filename = os.path.basename(self.manifest)
        name, ext = os.path.splitext(filename)
        
        # Handle .txt extension if present (or others)
        suffix_under = f"_u{self.cutoff}"
        suffix_over = f"_o{self.cutoff}"
        
        out_under = os.path.join(base_dir, f"{name}{suffix_under}{ext}")
        out_over = os.path.join(base_dir, f"{name}{suffix_over}{ext}")
        
        print(f"Creating new metadata files:\n  {out_under}\n  {out_over}")
        
        try:
            with open(self.manifest, 'r') as fin, \
                 open(out_under, 'w') as funder, \
                 open(out_over, 'w') as fover:
                
                for line in fin:
                    if not line.strip() or line.startswith('#'):
                        funder.write(line)
                        fover.write(line)
                        continue
                    
                    # Regex to match the first word
                    pattern = re.compile(r'^([^\s]+)')
                    
                    line_under = pattern.sub(r'\1' + suffix_under, line)
                    line_over = pattern.sub(r'\1' + suffix_over, line)
                    
                    funder.write(line_under)
                    fover.write(line_over)
                    
        except Exception as e:
            print(f"Error creating metadata files: {e}")
            sys.exit(1)

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
        
        suffix_under = f"_u{self.cutoff}"
        suffix_over = f"_o{self.cutoff}"
        
        out_under = os.path.join(dir_under, f"{sample}{suffix_under}.bam")
        out_over = os.path.join(dir_over, f"{sample}{suffix_over}.bam")
        
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

    def process(self):
        print("Generating new metadata files...")
        self._create_manifests()
        
        print(f"Processing BAM files in {self.bamdir}...")
        
        # If user provided threads, we can use them.
        threads = self.args.threads if hasattr(self.args, 'threads') and self.args.threads > 0 else 1
        
        if threads > 1:
            with multiprocessing.Pool(threads) as pool:
                pool.map(self._split_bam, self.samples)
        else:
            for sample in self.samples:
                self._split_bam(sample)
                
        print("Done.")
