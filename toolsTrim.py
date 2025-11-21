#!/usr/bin/env python3

import os
import sys
import json
import subprocess
import pandas as pd
import multiprocessing
import plotsTrimmingStats

class FastpTrimmer:
    '''
    Class to handle adapter trimming, merging, and UMI extraction using fastp.
    '''
    def __init__(self, args):
        self.args = args
        self.manifest = args.manifest
        
        # Parse Manifest
        self.samples = self._parse_manifest()

        # Thread management: 
        # fastp uses threads internally (-w). 
        # We also want to run samples in parallel if possible.
        # Heuristic: Divide total available threads by samples, but keep fastp threads reasonable (<=16).
        total_cores = self.args.threads
        if total_cores == 0:
            total_cores = multiprocessing.cpu_count()
            
        self.jobs = min(len(self.samples), total_cores)
        # fastp threads per job. Ensure at least 1.
        self.fastp_threads = max(1, int(total_cores / self.jobs))
        
        # Cap fastp threads to prevent diminishing returns/overhead if running many single jobs
        if self.fastp_threads > 16: 
            self.fastp_threads = 16

    def _parse_manifest(self):
        '''
        Parses tab-delimited manifest file.
        Expected format: OutputPrefix <tab> R1_Path [<tab> R2_Path]
        '''
        samples = {}
        try:
            with open(self.manifest, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split()
                    output_prefix = parts[0]
                    r1 = parts[1]
                    r2 = parts[2] if len(parts) > 2 else None
                    
                    if not os.path.isfile(r1):
                        raise FileNotFoundError(f"Read 1 file not found: {r1}")
                    if r2 and not os.path.isfile(r2):
                        raise FileNotFoundError(f"Read 2 file not found: {r2}")
                        
                    samples[output_prefix] = {'r1': r1, 'r2': r2}
        except Exception as e:
            print(f"Error parsing manifest file: {e}")
            sys.exit(1)
        return samples

    def _construct_command(self, output_prefix, files):
        '''
        Constructs the fastp command line arguments.
        '''
        r1 = files['r1']
        r2 = files['r2']
        
        # Determine output directory and base name
        output_dir = os.path.dirname(output_prefix)
        
        # Ensure output directory exists
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
        
        # Base input/output
        cmd = [
            'fastp',
            '--in1', r1,
            '--thread', str(self.fastp_threads),
            '--length_required', str(self.args.length),
            '--json', f"{output_prefix}.json",
            '--html', f"{output_prefix}.html"
        ]

        # Adapter settings (fastp auto-detects, but overrides help)
        if self.args.adapter1:
            cmd.extend(['--adapter_sequence', self.args.adapter1])
        if self.args.adapter2:
            cmd.extend(['--adapter_sequence_r2', self.args.adapter2])

        # UMI Processing
        if self.args.umilength > 0:
            cmd.append('--umi')
            cmd.extend(['--umi_len', str(self.args.umilength)])
            
            # UMI Location logic
            if self.args.umi3:
                # If UMI is at 3' end, fastp doesn't have a direct 'read1_tail' locator in same way
                # but usually 3' adapters handle this. If strictly UMI:
                # Assuming standard UMI-tools logic where it's part of the read sequence
                # fastp default is 5' (--umi_loc=read1). 
                # For 3' UMI, we might rely on adapter trimming or specific location flags if supported by version.
                # Standard fastp handles 5' well. For 3', we warn or assume user configured adapter trimming to expose it.
                # However, common protocols (McSeq) often put UMI at 5'. 
                # If strictly 3', use --umi_loc read1/read2 combined with skip if needed.
                pass 
            else:
                # Default is 5' end of read1
                cmd.extend(['--umi_loc', 'read1'])

        # Paired-end specific settings
        if r2:
            cmd.extend(['--in2', r2])
            
            # Merging
            # Output merged file
            merged_out = f"{output_prefix}_merged.fastq.gz"
            cmd.extend(['--merge', '--merged_out', merged_out])
            
            # Unmerged outputs (optional, but good for debug)
            # If not specified, fastp discards them or puts them in out1/out2.
            # We usually want the merged ones for tRNA-seq.
            # To keep unmerged:
            out1 = f"{output_prefix}_unmerged_R1.fastq.gz"
            out2 = f"{output_prefix}_unmerged_R2.fastq.gz"
            cmd.extend(['--out1', out1, '--out2', out2])
            
            # Adapter detection is usually better enabled for PE
            cmd.append('--detect_adapter_for_pe')

        else:
            # Single-end output
            out1 = f"{output_prefix}_trimmed.fastq.gz"
            cmd.extend(['--out1', out1])

        return cmd

    def _run_process(self, output_prefix, files):
        '''
        Worker function for multiprocessing
        '''
        cmd = self._construct_command(output_prefix, files)
        cmd_str = ' '.join(cmd)
        
        if self.args.verbose:
            print(f"[{output_prefix}] Starting: {cmd_str}")
        
        try:
            # Check if fastp is installed
            process = subprocess.run(cmd, capture_output=True, text=True)
            if process.returncode != 0:
                return (output_prefix, False, process.stderr)
            return (output_prefix, True, "Success")
        except FileNotFoundError:
            return (output_prefix, False, "fastp executable not found in PATH.")

    def process(self):
        '''
        Main execution block
        '''
        print(f"Starting trimming on {len(self.samples)} samples.")
        print(f"Configuration: {self.jobs} concurrent jobs, {self.fastp_threads} threads per job.")
        
        tasks = [(name, files) for name, files in self.samples.items()]
        
        with multiprocessing.Pool(self.jobs) as pool:
            results = pool.starmap(self._run_process, tasks)

        # Check results
        failed = []
        for name, success, msg in results:
            if not success:
                print(f"ERROR: Sample {name} failed.\nMessage: {msg}")
                failed.append(name)
            else:
                print(f"Finished: {name}")
        
        if failed:
            print(f"\nWarning: {len(failed)} samples failed processing.")
        
        self._generate_summary()

    def _generate_summary(self):
        '''
        Parses fastp JSON reports and creates a summary DataFrame/CSV
        similar to the old logic but more robust.
        '''
        summary_data = []
        
        print("Generating summary report...")
        for output_prefix in self.samples:
            json_path = f"{output_prefix}.json"
            if not os.path.exists(json_path):
                continue
                
            with open(json_path, 'r') as f:
                data = json.load(f)
                
                # General stats
                row = {'Sample': os.path.basename(output_prefix)}
                summary = data.get('summary', {})
                before = summary.get('before_filtering', {})
                after = summary.get('after_filtering', {})
                
                row['Raw_Reads'] = before.get('total_reads', 0)
                row['Clean_Reads'] = after.get('total_reads', 0)
                row['Reads_Passed_Filter'] = row['Clean_Reads']
                
                # Filtering stats
                filter_stats = data.get('filtering_result', {})
                row['Reads_Too_Short'] = filter_stats.get('too_short_reads', 0)
                row['Reads_Low_Quality'] = filter_stats.get('low_quality_reads', 0)
                
                # Adapter stats
                adapter_stats = data.get('adapter_cutting', {})
                row['Reads_With_Adapter'] = adapter_stats.get('adapter_trimmed_reads', 0)
                
                # Merging stats (if paired)
                if 'merged_and_filtered' in data:
                    merge_stats = data['merged_and_filtered']
                    row['Merged_Reads'] = merge_stats.get('total_reads', 0)
                    row['Unmerged_Reads'] = row['Clean_Reads'] - row['Merged_Reads'] 
                else:
                    row['Merged_Reads'] = 0
                    row['Unmerged_Reads'] = row['Clean_Reads']

                summary_data.append(row)

        if summary_data:
            df = pd.DataFrame(summary_data)
            
            # Determine output directory from the first sample
            first_output_prefix = list(self.samples.keys())[0]
            output_dir = os.path.dirname(first_output_prefix)
            
            # Create manifest update file (Sample -> Final Output)
            manifest_out = os.path.join(output_dir, f"{self.args.runname}_trim_manifest_updated.txt")
            with open(manifest_out, 'w') as f:
                for output_prefix in self.samples:
                    # Logic to determine the 'primary' output file
                    if self.samples[output_prefix]['r2']:
                         # If paired, primary for tRNA is usually the merged file
                         outfile = f"{output_prefix}_merged.fastq.gz"
                    else:
                         outfile = f"{output_prefix}_trimmed.fastq.gz"
                    f.write(f"{os.path.basename(output_prefix)}\t{outfile}\n")

            # Save Stats
            stats_out = os.path.join(output_dir, f"{self.args.runname}_trim_stats.csv")
            df.to_csv(stats_out, index=False)
            print(f"Summary statistics written to {stats_out}")
            print(f"Updated manifest written to {manifest_out}")
            
            # Generate Plot
            plot_out = os.path.join(output_dir, f"{self.args.runname}_trim_feature_types.pdf")
            print("Generating feature types plot...")
            plotsTrimmingStats.visualizer(stats_out, plot_out).plot()
        else:
            print("No JSON reports found to summarize.")

if __name__ == "__main__":
    pass