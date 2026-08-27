#!/usr/bin/env python3

import os
import sys
import json
import logging
import subprocess
import pandas as pd
import multiprocessing
from pydantic import ValidationError
from . import plotsTrimmingStats, toolsTG

# Where trimmed output lands when a manifest's OutputPrefix carries no directory of its own.
DEFAULT_TRIM_DIR = 'processed/trimmed'

class FastpTrimmer:
    '''
    Class to handle adapter trimming, merging, and UMI extraction using fastp.
    '''
    def __init__(self, args):
        self.args = args
        self.manifest = args.input
        # A stable, module-level logger -- NOT configured with its own handlers here. Per
        # Python's logging convention, only the application entry point (cli.py's
        # configure_logging()/handle_output()) attaches FileHandler/StreamHandler; this module
        # just logs and relies on propagation. Using getLogger(__name__) (module-scoped, not
        # per-instance) also means repeated FastpTrimmer construction reuses the same cached
        # logger object rather than registering a new one every time -- loggers are never
        # garbage-collected by the logging module, so an instance-keyed name would otherwise
        # leak a logger (and any handler it held) for the life of the process.
        self.logger = logging.getLogger(__name__)

        # Load the style file if specified. Only the colors block is used here, namespaced
        # under 'trimtype' since there is no adata/obs column at trim time -- but it is the
        # same file `analyze graph` takes, so one file styles the whole pipeline.
        self.colormap = None
        self.style = None
        if getattr(args, 'style', None):
            self.style = toolsTG.load_style_file(args.style, self.logger)
            self.colormap = self.style.colors_for('trimtype')

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

                    # A bare OutputPrefix (no directory component) writes into DEFAULT_TRIM_DIR.
                    # Resolving that here -- once, before anything is keyed on it -- is what
                    # keeps self.samples pointing at the real output location. Deriving it
                    # again downstream is what previously left _generate_summary() hunting for
                    # fastp's JSON report in the working directory and finding nothing.
                    if not os.path.dirname(output_prefix):
                        output_prefix = os.path.join(DEFAULT_TRIM_DIR, output_prefix)

                    if not os.path.isfile(r1):
                        raise FileNotFoundError(f"Read 1 file not found: {r1}")
                    if r2 and not os.path.isfile(r2):
                        raise FileNotFoundError(f"Read 2 file not found: {r2}")
                        
                    samples[output_prefix] = {'r1': r1, 'r2': r2}
        except Exception as e:
            self.logger.error(f"Error parsing manifest file: {e}")
            sys.exit(1)
        return samples

    @staticmethod
    def _primary_output(output_prefix, files):
        '''
        The single definition of what fastp's primary output for a sample is called: paired-end
        input is merged into one file, single-end input is simply trimmed. Both
        _construct_command() (which tells fastp where to write) and _generate_summary() (which
        records the path in trim_metadata.tsv) read the rule from here instead of each spelling
        it out -- the two spellings had drifted, and the metadata template named a
        `_merged_trimmed.fastq.gz` that fastp never writes.
        '''
        if files['r2']:
            return f"{output_prefix}_merged.fastq.gz"
        return f"{output_prefix}_trimmed.fastq.gz"

    def _construct_command(self, output_prefix, files):
        '''
        Constructs the fastp command line arguments.
        '''
        r1 = files['r1']
        r2 = files['r2']
        
        # _parse_manifest() has already resolved a bare prefix to DEFAULT_TRIM_DIR, so the
        # prefix arrives fully qualified here -- only its directory still needs creating.
        output_dir = os.path.dirname(output_prefix)
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
        # fastp's own --umi_loc only supports head-anchored extraction (index1/index2/read1/
        # read2/per_index/per_read) -- there is no tail-anchored option in fastp itself, so a
        # --umi3 request can't be expressed as a fastp flag at all. For that case, fastp runs
        # with no UMI flags whatsoever (plain adapter trimming/merging), and a umi_tools
        # extract --3prime post-processing step (mirroring tRAX's own approach) handles the
        # actual 3' extraction afterward -- see _run_process().
        if self.args.umilength > 0 and not self.args.umi3:
            cmd.append('--umi')
            cmd.extend(['--umi_len', str(self.args.umilength)])
            # Default is 5' end of read1
            cmd.extend(['--umi_loc', 'read1'])
            # fastp's own default delimiter is ':', but `umi_tools` (both the --umi3 extraction
            # path below and `preprocess map --dedup`) uses '_', and a colon is ambiguous anyway
            # because Illumina read names are already colon-delimited. Pinning it here keeps one
            # read-name shape across both of tRNAgraph's UMI paths.
            cmd.extend(['--umi_delim', '_'])

        # Paired-end specific settings
        if r2:
            cmd.extend(['--in2', r2])

            # Merging
            # Output merged file
            merged_out = self._primary_output(output_prefix, files)
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
            primary_output = merged_out

        else:
            # Single-end output
            out1 = self._primary_output(output_prefix, files)
            cmd.extend(['--out1', out1])
            primary_output = out1

        return cmd, primary_output

    def _run_process(self, output_prefix, files):
        '''
        Worker function for multiprocessing
        '''
        cmd, primary_output = self._construct_command(output_prefix, files)
        cmd_str = ' '.join(cmd)

        if self.args.verbose:
            self.logger.info(f"[{output_prefix}] Starting: {cmd_str}")

        try:
            # Check if fastp is installed
            process = subprocess.run(cmd, capture_output=True, text=True)
        except FileNotFoundError:
            return (output_prefix, False, "fastp executable not found in PATH.")

        if process.returncode != 0:
            # Return fastp's own stderr regardless of success/failure -- it's discarded here
            # (rather than in process(), which runs sequentially in the parent after this
            # worker returns) since capturing/logging it must not happen inside a
            # multiprocessing.Pool worker, where concurrent writes to the same --log file
            # could interleave.
            return (output_prefix, False, process.stderr)

        stderr = process.stderr
        if self.args.umilength > 0 and self.args.umi3:
            umi_ok, umi_stderr = self._extract_three_prime_umi(output_prefix, primary_output)
            stderr += umi_stderr
            if not umi_ok:
                return (output_prefix, False, stderr)

        return (output_prefix, True, stderr)

    def _extract_three_prime_umi(self, output_prefix, fastq_path):
        '''
        Post-processing step for --umi3: fastp has no tail-anchored --umi_loc, so
        _construct_command() runs fastp with no UMI flags at all, and this runs
        `umi_tools extract --3prime` on fastp's own trimmed/merged output afterward,
        replacing it in place -- mirroring tRAX's original trimadapters.py, which ran the
        equivalent `umi_tools extract` (with `--3prime` for this same case) as a
        post-processing step after SeqPrep/cutadapt, over the same output slot.
        '''
        umi_log = f"{output_prefix}_umilog.txt"
        # umi_tools decides whether to gzip its output from the --stdout filename's own
        # extension, not from the input's -- the temp name must still end in .gz or the
        # written bytes are plain text despite os.replace() giving them a .fastq.gz name.
        extracted = fastq_path.replace('.fastq.gz', '_umi_extracted.fastq.gz')
        cmd = [
            'umi_tools', 'extract', '--3prime',
            '--bc-pattern', 'N' * self.args.umilength,
            '--stdin', fastq_path,
            '--stdout', extracted,
            '-L', umi_log,
        ]
        try:
            process = subprocess.run(cmd, capture_output=True, text=True)
        except FileNotFoundError:
            return False, "umi_tools executable not found in PATH."
        if process.returncode != 0:
            return False, process.stderr
        os.replace(extracted, fastq_path)
        return True, process.stderr

    def _run_process_unpacked(self, task):
        '''
        pool.imap_unordered() calls its worker function with a single argument, unlike
        starmap's automatic tuple-unpacking -- this just adapts _run_process() to that.
        '''
        name, files = task
        return self._run_process(name, files)

    def process(self):
        '''
        Main execution block
        '''
        self.logger.info(f"Starting trimming on {len(self.samples)} samples.")
        self.logger.info(f"Configuration: {self.jobs} concurrent jobs, {self.fastp_threads} threads per job.")

        tasks = [(name, files) for name, files in self.samples.items()]

        # imap_unordered (rather than starmap) so progress_iterator gets a real per-sample
        # completion signal to report on -- starmap blocks and hands back every result at once,
        # with no way to know how far through a many-hour run it is. Safe to consume out of
        # task order: the loop below reads `name` out of each result tuple rather than relying
        # on positional alignment with `tasks`.
        with multiprocessing.Pool(self.jobs) as pool:
            results = list(toolsTG.progress_iterator(
                pool.imap_unordered(self._run_process_unpacked, tasks),
                total=len(tasks), desc="Trimming samples", logger=self.logger,
                quiet=getattr(self.args, 'quiet', False),
            ))

        # Check results. fastp's own stderr (per-sample summary/diagnostic output) is logged
        # here -- back in the parent process, sequentially, after the pool returns -- for every
        # sample regardless of success/failure, so it's persisted to the run's log file even on
        # success instead of discarding it entirely.
        failed = []
        for name, success, stderr in results:
            if not success:
                self.logger.error(f"Sample {name} failed.\nMessage: {stderr}")
                failed.append(name)
            else:
                self.logger.info(f"Finished: {name}")
                if stderr:
                    self.logger.info(f"[{name}] fastp stderr:\n{stderr}")

        if failed:
            self.logger.warning(f"{len(failed)} samples failed processing.")

        self._generate_summary()

    def _generate_summary(self):
        '''
        Parses fastp JSON reports and creates a summary DataFrame/CSV
        similar to the old logic but more robust.
        '''
        summary_data = []
        # Prefixes fastp actually produced a report for. A sample whose run failed has no
        # trimmed output on disk, so it must not appear in the metadata template -- that would
        # hand `analyze build` a path to a file that isn't there.
        summarized = []

        self.logger.info("Generating summary report...")
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
                
                # Merging stats -- fastp only attempts merging for paired-end input, so a
                # single-end sample has no "unmerged" concept at all (every filter-passing
                # read is simply trimmed, never a merge candidate in the first place).
                if 'merged_and_filtered' in data:
                    merge_stats = data['merged_and_filtered']
                    row['Merged_Reads'] = merge_stats.get('total_reads', 0)
                    row['Unmerged_Reads'] = row['Clean_Reads'] - row['Merged_Reads']
                    row['Trimmed_Reads'] = 0
                else:
                    row['Merged_Reads'] = 0
                    row['Unmerged_Reads'] = 0
                    row['Trimmed_Reads'] = row['Clean_Reads']

                summary_data.append(row)
                summarized.append(output_prefix)

        if summary_data:
            df = pd.DataFrame(summary_data)
            
            # One run produces one summary, placed alongside the first sample that
            # produced output. Every prefix is directory-qualified by now, so a manifest
            # writing into several directories still gets a single template -- say which
            # directory it landed in rather than silently picking the first.
            output_dirs = list(dict.fromkeys(os.path.dirname(p) for p in summarized))
            output_dir = output_dirs[0]

            # Create manifest update file (Sample -> Final Output)
            manifest_out = os.path.join(output_dir, "trim_metadata.tsv")
            if len(output_dirs) > 1:
                self.logger.warning(
                    f"Manifest writes trimmed output into {len(output_dirs)} directories "
                    f"({', '.join(output_dirs)}); a single summary covering all of them "
                    f"is written to {output_dir}, including {manifest_out}."
                )
            with open(manifest_out, 'w') as f:
                # write header
                f.write("fastq\tsample\tgroup\n")
                for output_prefix in summarized:
                    outfile = self._primary_output(output_prefix, self.samples[output_prefix])
                    f.write(f"{outfile}\t{os.path.basename(output_prefix)}\t{os.path.basename(output_prefix)}\n")

            # Save Stats
            stats_out = os.path.join(output_dir, "trim_stats.csv")
            df.to_csv(stats_out, index=False)
            self.logger.info(f"Summary statistics written to {stats_out}")
            self.logger.info(f"Updated manifest written to {manifest_out}")

            # Generate Plot
            plot_out = os.path.join(output_dir, "trim_feature_types.pdf")
            self.logger.info("Generating feature types plot...")
            plotsTrimmingStats.visualizer(stats_out, plot_out, colormap=self.colormap, settings=toolsTG.resolve_plot_style(self.style, 'trimming')).plot()
        else:
            self.logger.warning("No JSON reports found to summarize.")

if __name__ == "__main__":
    pass