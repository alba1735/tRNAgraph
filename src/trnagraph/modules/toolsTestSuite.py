#!/usr/bin/env python3

import os
import re
import sys
import subprocess
import logging
import contextlib
import zipfile
import argparse


class _LiveBoxHandler(logging.Handler):
    '''
    Feeds demoPipeline._live_box()'s live display: captures every log message into a bounded
    tail (the "Recent activity" panel) and watches for a milestone format to drive the box's bar
    from real percentages whenever a wrapped subcommand emits one -- a subprocess's piped stdout
    is never a real terminal from its own point of view, regardless of whether THIS process is
    interactive, so toolsTG.py's progress_iterator/PhaseTracker always fall back to logging one of
    the two formats below rather than a live rich display.

    Two distinct formats exist, and this handler must tell them apart:

    - toolsTG.PhaseTracker's OUTER, phase-level format -- the literal word "phase" immediately
      before the fraction is the distinguishing marker, present in both of its two variants:
      the final "<desc> phase N/Total (P%) complete: <label>" (e.g. "Build phase 2/6 (33%)
      complete: Analyzing counts") logged once a phase's `with phase():` block exits, and the
      intermediate "<desc> phase N/Total (P%): <label>" (no "complete") that `advance()` logs
      periodically for a phase with a large weight (e.g. one per coverage plot, throttled to
      ~10%-of-that-phase's-weight steps).
    - toolsTG.progress_iterator's bare, per-item format ("N/Total (P%) complete", e.g.
      "Counting reads: 7/10 (70%) complete"), used both standalone (trim/map, which have no
      phase concept) and for the INNER per-sample loop within a phase-tracked command's
      "Counting Reads"/"Counting Read Types" phases.

    `phase_only=True` (set by whichever `_live_box()` call wraps a phase-tracked command, e.g.
    "Building AnnData object...") makes bare per-item milestones never drive the bar at all --
    they're still captured into the tail, just not used to move it. This has to be an upfront
    flag (known from which step this is), not something inferred reactively from whichever
    milestone type happens to arrive first: a phase-tracked command's FIRST phase can itself wrap
    an inner per-item loop (e.g. "Counting Reads"), whose own milestones would otherwise fire
    *before* that phase's own completion line does -- reactively waiting for "the first phase
    signal seen so far" would let the bar climb through the inner loop's 10/20/.../100% and then
    visibly drop back down to that phase's real (much lower) outer percentage once its
    completion line finally arrives, a jarring glitch. `phase_only=False` (the default, used by
    every other step -- trim/map/etc, which have no phase concept at all) leaves the original
    behavior fully intact: bare per-item milestones always drive the bar, exactly as before.

    This exists to fix a real regression: without `phase_only`, an inner per-sample counting
    milestone reaching "10/10 (100%) complete" almost immediately during a build used to pin the
    box at 100% through everything that ran afterward (DESeq2 fitting, coverage generation, VST,
    writing the h5ad), since the old bare-only regex matched both formats identically and simply
    took whichever line came last.
    '''
    _PHASE_MILESTONE_RE = re.compile(r'\bphase (\d+)/(\d+) \(\d+%\)(?: complete)?:\s*(.*)$')
    _ITEM_MILESTONE_RE = re.compile(r'(\d+)/(\d+) \(\d+%\) complete')

    def __init__(self, tail, on_milestone, phase_only: bool = False):
        super().__init__()
        self.tail = tail
        self.on_milestone = on_milestone
        self.phase_only = phase_only

    def emit(self, record: logging.LogRecord) -> None:
        message = self.format(record)
        self.tail.append(message)

        phase_match = self._PHASE_MILESTONE_RE.search(message)
        if phase_match:
            self.on_milestone(int(phase_match.group(1)), int(phase_match.group(2)), phase_match.group(3))
            return

        if self.phase_only:
            return

        item_match = self._ITEM_MILESTONE_RE.search(message)
        if item_match:
            self.on_milestone(int(item_match.group(1)), int(item_match.group(2)))


_VARIANT_LABEL_RE = re.compile(r'^\[(Under|Over) (\d+)\]')


def _friendly_variant_title(label: str):
    '''
    Maps a PhaseTracker phase label's split-variant bracket prefix (e.g. "[Under 60] Counting
    Reads", from AnalysisPipeline's variant_label="Under 60") to a friendlier box title -- "Under"
    means the fragment (u<N>) variant, "Over" the full-length (o<N>) one, matching roadmap.md's
    Phase 4 terminology. Returns None for a non-variant label (the main/full build), so the box's
    title is left as whatever _live_box() was originally given (e.g. "Building AnnData object...").
    '''
    match = _VARIANT_LABEL_RE.match(label)
    if not match:
        return None
    direction, cutoff = match.group(1), match.group(2)
    kind = 'Fragments' if direction == 'Under' else 'Full-length'
    return f'Building {direction.lower()} {cutoff}bp split ({kind})...'


class demoPipeline:
    """
    A pipeline class to run the tRNAgraph test suite.
    
    This class automates the download and processing of sample data (Vibrio cholerae)
    to verify the functionality of the tRNAgraph pipeline. It mirrors the steps
    found in the tRAX tutorial.
    """
    def __init__(self, args: argparse.Namespace) -> None:
        """
        Initialize the demoPipeline.

        Args:
            args (argparse.Namespace): Parsed command-line arguments.
        """
        self.args = args
        self.repo_root = os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))
        self.trnagraph_path = "trnagraph"
        
        # Set up working directory
        if self.args.directory:
             work_dir = os.path.abspath(self.args.directory)
        else:
             work_dir = os.path.join(self.repo_root, "test_vibrChol1")

        os.makedirs(work_dir, exist_ok=True)
        os.chdir(work_dir)

        # Configure logging. Deliberately self-contained (its own handlers, not
        # logging.basicConfig() on the root logger) with propagate=False: this logger's name
        # (__name__ == 'trnagraph.modules.toolsTestSuite') makes it a natural child of the
        # 'trnagraph' logger that cli.py's configure_logging() centrally configures for --log/
        # --quiet, and without propagate=False every message here would also flow up into
        # whatever that attached, double-printing everything (this class logs everything through
        # self.logger, no separate print() calls). toolsTestSuite.log is a fixed, dedicated
        # demo-pipeline log, always written regardless of --log/--quiet -- but this class's OWN
        # --quiet still controls its own separate console StreamHandler below, so
        # _run_command()'s live-streamed subprocess output is actually visible on screen while a
        # step runs, not just captured in the log file.
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.logger.propagate = False
        for handler in self.logger.handlers[:]:
            handler.close()
            self.logger.removeHandler(handler)
        file_handler = logging.FileHandler('toolsTestSuite.log', mode='w')
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        self.logger.addHandler(file_handler)
        # Kept as its own attribute (not just "whatever StreamHandler is on self.logger") so
        # _run_command's rich Live branch can detach exactly this one handler for the duration
        # of a live session -- logging.FileHandler is itself a StreamHandler subclass, so a
        # generic isinstance(h, logging.StreamHandler) check would also match file_handler.
        self.console_handler = None
        if not getattr(self.args, 'quiet', False):
            self.console_handler = logging.StreamHandler(sys.stdout)
            self.console_handler.setFormatter(logging.Formatter('%(message)s'))
            self.logger.addHandler(self.console_handler)

        # Clean up if requested
        if self.args.all:
            self._cleanup_workspace()

        # Copy assets
        self.assets_dir = os.path.join(self.repo_root, "src/trnagraph/assets")
        os.makedirs("config", exist_ok=True)
        with self._live_box("Copying assets..."):
            self._run_command(f"cp --update {self.assets_dir}/*.txt config/.")
            self._run_command(f"cp --update {self.assets_dir}/*.json config/.")

    @contextlib.contextmanager
    def _live_box(self, description: str, phase_only: bool = False):
        """
        Wrap one whole step (e.g. a single download_*()/trim_fastq()/etc. call, including every
        self.logger.info() and _run_command() subprocess it runs) in a single live rich display:
        a spinner, switching to a real progress bar the moment a wrapped step reports a
        progress_iterator()-style milestone line, plus a scrolling "Recent activity" tail of
        whatever gets logged during the step -- so the terminal shows one clean, auto-refreshing
        box per step instead of scrolling raw log lines. Falls back to plain logging (no display
        at all) when this isn't a real interactive terminal, or --quiet/--log is set.

        `phase_only=True` for a step whose wrapped command is phase-tracked (toolsTG.PhaseTracker,
        e.g. `analyze build`) -- see _LiveBoxHandler for why this must be an upfront flag set by
        the caller (who knows what it's wrapping) rather than inferred from whichever milestone
        format happens to arrive first.
        """
        use_live = sys.stdout.isatty() and self.console_handler is not None and not getattr(self.args, 'log', None)
        if not use_live:
            self.logger.info(description)
            yield
            return

        from collections import deque
        from rich.console import Group
        from rich.live import Live
        from rich.panel import Panel
        from rich.progress import BarColumn, Progress, TaskProgressColumn, TextColumn
        from rich.spinner import Spinner

        tail = deque(maxlen=10)
        spinner = Spinner('dots', text=description, style='green')
        progress = Progress(
            TextColumn('[bold green]{task.description}'),
            BarColumn(complete_style='green', finished_style='bright_green'),
            TaskProgressColumn(),
        )
        task_id = progress.add_task(description, total=None)
        state = {'has_milestone': False}

        def on_milestone(completed, total, label=None):
            # Only a split-variant phase label overrides the box's title (e.g. "Building under
            # 60bp split (Fragments)...") -- the main/full build's phases keep the original
            # `description` the box was opened with, since _friendly_variant_title() returns
            # None for those.
            update_kwargs = {'completed': completed, 'total': total}
            friendly_title = _friendly_variant_title(label) if label else None
            if friendly_title:
                update_kwargs['description'] = friendly_title
            progress.update(task_id, **update_kwargs)
            state['has_milestone'] = True

        def render():
            tail_text = '\n'.join(tail) or '...'
            header = progress if state['has_milestone'] else spinner
            return Group(header, Panel(tail_text, title='Recent activity', border_style='cyan', height=12))

        box_handler = _LiveBoxHandler(tail, on_milestone, phase_only=phase_only)
        box_handler.setFormatter(logging.Formatter('%(message)s'))

        # Swap the console handler out for the box's capturing handler for the whole step --
        # every self.logger.info() call made anywhere during this block (including inside any
        # _run_command() subprocess streaming) still reaches toolsTestSuite.log via
        # file_handler (untouched throughout), just routes into the box instead of the raw
        # terminal, which the Live display now owns exclusively (concurrent raw prints would
        # corrupt it). transient=True clears the box when the step finishes, so a full run
        # shows one clean box per step rather than stacking a permanent panel per step down the
        # terminal until it scrolls off the bottom. get_renderable (rather than a static
        # renderable) makes Live auto-refresh from current tail/bar state on its own, without
        # needing an explicit .update() call after every single log line.
        self.logger.removeHandler(self.console_handler)
        self.logger.addHandler(box_handler)
        try:
            with Live(get_renderable=render, refresh_per_second=4, transient=True):
                self.logger.info(description)
                yield
        finally:
            self.logger.removeHandler(box_handler)
            self.logger.addHandler(self.console_handler)

    def _run_command(self, command: str, description: str = "", check: bool = True) -> subprocess.CompletedProcess:
        """
        Helper method to run shell commands with logging.

        Args:
            command (str): The shell command to execute.
            description (str): A brief description for logging/printing.
            check (bool): Whether to raise an exception on failure.

        Returns:
            subprocess.CompletedProcess: The result of the command execution.
        """
        if description:
            self.logger.info(description)

        # Ensure the current python environment's bin directory is in the PATH
        env = os.environ.copy()
        python_bin_dir = os.path.dirname(sys.executable)
        if python_bin_dir not in env.get("PATH", ""):
            env["PATH"] = f"{python_bin_dir}{os.pathsep}{env.get('PATH', '')}"

        # Streamed rather than subprocess.run(capture_output=True): that fully buffers the
        # child process's entire output until it exits, then dumps it all at once -- meaning
        # nothing from the trim/map/build/graph subcommands this wraps was ever visible while a
        # step was actually running. stderr is merged into stdout (a single interleaved stream)
        # to avoid the complexity of reading two pipes concurrently without deadlocking; nothing
        # in this file inspects them separately. Each line just goes through self.logger --
        # whatever's currently attached (the box's capturing handler if a _live_box() is active,
        # or the plain console handler otherwise) decides where it actually ends up; this method
        # doesn't need to know or care which.
        process = subprocess.Popen(
            command, shell=True, env=env, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1,
        )
        output_lines = []
        for line in process.stdout:
            line = line.rstrip('\n')
            output_lines.append(line)
            self.logger.info(line)

        process.wait()
        full_output = '\n'.join(output_lines)

        if check and process.returncode != 0:
            self.logger.error(f"Command failed: {command}\nOutput:\n{full_output}")
            raise subprocess.CalledProcessError(process.returncode, command, output=full_output)
        return subprocess.CompletedProcess(command, process.returncode, stdout=full_output, stderr='')

    def _cleanup_workspace(self) -> None:
        """Removes generated files to ensure a clean run, keeping only the log file."""
        with self._live_box("Cleaning up workspace..."):
            # Remove all files and directories in the current working directory (the test
            # directory) except for the log file.
            self._run_command(
                'find . -maxdepth 1 -mindepth 1 -not -name "toolsTestSuite.log" -exec rm -rf {} +',
                "Removing all contents from test directory except log file..."
            )
            self.logger.info("Workspace cleaned.")

    def download_metadata(self) -> None:
        """Downloads metadata from SRA using pysradb."""
        with self._live_box("Downloading metadata from SRA..."):
            os.makedirs("raw/vibrChol1/fastq", exist_ok=True)

            # Unlike download_fastq/download_trna/download_genome, this had no skip-if-present
            # check at all -- every run re-hit pysradb's network call even when accessions.tsv
            # was already there from a prior run, making it the slowest step to redundantly repeat.
            if os.path.exists("raw/vibrChol1/fastq/accessions.tsv"):
                self.logger.info("accessions.tsv already exists, skipping metadata re-fetch.")
                self.logger.info("Done.")
                return

            cmd = (
                "pysradb metadata SRP254278 | "
                "grep -v -e 'Escherichia coli' -e 'trmK' -e 'miaA' -e 'ttcA' -e 'thiI' -e 'run_accession' | "
                "cut -f22 > raw/vibrChol1/fastq/accessions.tsv"
            )
            self._run_command(cmd, "Fetching metadata...")
            self.logger.info("Done.")

    def download_fastq(self) -> None:
        """Downloads FASTQ files for the accessions listed in accessions.tsv."""
        with self._live_box("Downloading fastq files..."):
            try:
                with open("raw/vibrChol1/fastq/accessions.tsv", "r") as f:
                    accessions = f.read().splitlines()

                for acc in accessions:
                    if os.path.exists(f"raw/vibrChol1/fastq/{acc}.fastq"):
                        self.logger.info(f"{acc}.fastq already exists, skipping download.")
                        continue

                    # Prefetch
                    self._run_command(
                        f"prefetch {acc} --output-file raw/vibrChol1/fastq/{acc}.sra",
                        f"Prefetching {acc}..."
                    )

                    # Fastq-dump
                    self._run_command(
                        f"fastq-dump -O raw/vibrChol1/fastq -X 100000 {acc}",
                        f"Dumping fastq for {acc}..."
                    )

                    # Cleanup SRA file
                    sra_file = f"raw/vibrChol1/fastq/{acc}.sra"
                    if os.path.exists(sra_file):
                        os.remove(sra_file)

            except FileNotFoundError:
                self.logger.error("accessions.tsv not found. Run download_metadata first.")
                raise

            self.logger.info("Done.")

    def download_trna(self) -> None:
        """Downloads and extracts Vibrio cholerae tRNA sequences."""
        with self._live_box("Downloading Vibrio cholerae tRNA sequences..."):
            os.makedirs("references/vibrChol1/trnas", exist_ok=True)
            self._run_command(f"cp --update {self.assets_dir}/*.gz .", "Copying assets...")

            if os.path.exists("references/vibrChol1/trnas/vibrChol1-tRNAs.fa"):
                self.logger.info("vibrChol1-tRNAs.fa already exists, skipping download.")
            else:
                # Using local tar.gz as per original logic - gtRNAdb has issues currently with downloads
                self._run_command(
                    "tar xzvf vibrChol1-tRNAs.tar.gz -C references/vibrChol1/trnas/",
                    "Extracting tRNA sequences..."
                )

            # cleanup
            self._run_command("rm -f vibrChol1-tRNAs.tar.gz", "Removing tRNA tar.gz...")

            self.logger.info("Done.")

    def download_genome(self) -> None:
        """Downloads Vibrio cholerae genome and GFF annotation, converting to GTF."""
        with self._live_box("Downloading Vibrio cholerae genome and GFF annotation..."):
            os.makedirs("references/vibrChol1/annotations", exist_ok=True)
            os.makedirs("references/vibrChol1/genomes", exist_ok=True)

            if os.path.exists("references/vibrChol1/annotations/GCF_000006745.1.gtf"):
                self.logger.info("Vibrio cholerae genome already exists, skipping download.")
                return

            # Download genome and GFF from FTP (API is flaky)
            ftp_base = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/745/GCF_000006745.1_ASM674v1"

            # Download fna.gz
            self._run_command(
                f"curl -s -o references/vibrChol1/genomes/genomic.fna.gz {ftp_base}/GCF_000006745.1_ASM674v1_genomic.fna.gz",
                "Downloading genome FASTA..."
            )
            self._run_command("gunzip -f references/vibrChol1/genomes/genomic.fna.gz", "Extracting genome FASTA...")

            # Download gff.gz
            self._run_command(
                f"curl -s -o references/vibrChol1/annotations/genomic.gff.gz {ftp_base}/GCF_000006745.1_ASM674v1_genomic.gff.gz",
                "Downloading genome GFF..."
            )
            self._run_command("gunzip -f references/vibrChol1/annotations/genomic.gff.gz", "Extracting genome GFF...")

            fna_path = "references/vibrChol1/genomes/genomic.fna"
            gff_path = "references/vibrChol1/annotations/genomic.gff"

            # Modify FASTA headers
            sed_cmd = (
                f'sed -i -e "/NC_002505.1/c\\>chrI" {fna_path} && '
                f'sed -i -e "/NC_002506.1/c\\>chrII" {fna_path}'
            )
            self._run_command(sed_cmd, "Modifying FASTA headers...")

            # Convert GFF to GTF
            gtf_path = "references/vibrChol1/annotations/genomic.gtf"
            self._run_command(f"gffread -E {gff_path} -T -o {gtf_path}", "Converting GFF to GTF...")

            # Modify GTF and filter safely
            filtered_gtf_path = "references/vibrChol1/annotations/genomic.filtered.gtf"
            final_gtf_cmd = (
                f"cat {gtf_path} | sed 's/NC_002505.1/chrI/g' | "
                f"sed 's/NC_002506.1/chrII/g' | grep -v '^#' > {filtered_gtf_path}"
            )
            self._run_command(final_gtf_cmd, "Finalizing GTF file...")
            # Atomically replace the original GTF with the filtered version
            os.replace(filtered_gtf_path, gtf_path)

            os.rename(fna_path, "references/vibrChol1/genomes/GCF_000006745.1_ASM674v1_genomic.fna")
            os.rename(gtf_path, "references/vibrChol1/annotations/GCF_000006745.1.gtf")

            # Cleanup GFF
            os.remove(gff_path)

            self.logger.info("Done.")

    def trim_fastq(self) -> None:
        """Trims adapters from FASTQ files using the tRNAgraph preprocess trim tool."""
        with self._live_box("Trimming fastq files with fastp..."):
            if os.path.exists("processed/trimmed/vibrChol1_trim_manifest_updated.txt"):
                self.logger.info("Trimmed fastq files already exist, skipping trimming.")
                return

            cmd = (
                f"{self.trnagraph_path} preprocess trim "
                "-i config/vibrChol1.manifest.txt -a1 ACTGTAGGCACCATCAATC --colormap config/colormap.json"
            )
            self._run_command(cmd, "Running trim command...")

            self.logger.info("Done.")

    def create_index(self) -> None:
        """Creates the Bowtie2 index for the tRNA database."""
        with self._live_box("Creating bowtie2 index..."):
            if os.path.exists("references/vibrChol1/trnadb/vibrChol1_db-tRNAgenome.1.bt2l"):
                self.logger.info("Bowtie2 index already exists, skipping creation.")
                return

            cmd = (
                f"{self.trnagraph_path} preprocess makedb "
                "-g references/vibrChol1/genomes/GCF_000006745.1_ASM674v1_genomic.fna "
                "-t references/vibrChol1/trnas/vibrChol1-tRNAs.out "
                "-r references/vibrChol1/trnas/vibrChol1-tRNAs.fa "
                "-m references/vibrChol1/trnas/vibrChol1-tRNAs_name_map.txt "
                "-s bact -o references/vibrChol1/trnadb/vibrChol1_db"
            )
            self._run_command(cmd, "Running makedb command...")

            self.logger.info("Done.")

    def map_reads(self) -> None:
        """Maps reads to the tRNA database."""
        with self._live_box("Mapping reads to tRNA genes..."):
            cmd = (
                f"{self.trnagraph_path} preprocess map "
                "-o vibrChol1 -d references/vibrChol1/trnadb/vibrChol1_db "
                "-i config/vibrChol1.metadata.txt "
                f"--bamdir processed/vibrChol1/bam"
            )
            self._run_command(cmd, "Running map command...")

            self.logger.info("Done.")

    def build_db(self) -> None:
        """Builds the AnnData object from the tRNAgraph output."""
        with self._live_box("Building AnnData object...", phase_only=True):
            extra_flags = ""
            if self.args.hubonly:
                extra_flags += " --hubonly"
            # Only add readlengthsplit if split_build is requested
            if getattr(self.args, 'split_build', False):
                extra_flags += " --readlengthsplit 60"

            cmd = (
                f"{self.trnagraph_path} analyze build "
                "-i config/vibrChol1.metadata.txt "
                "-d references/vibrChol1/trnadb/vibrChol1_db "
                "--gtf references/vibrChol1/annotations/GCF_000006745.1.gtf "
                "--pairs config/vibrChol1.pair.txt "
                "--bamdir processed/vibrChol1/bam "
                "--uniqueonly "
                "-o vibrChol1"
                f"{extra_flags}"
            )
            self._run_command(cmd, "Running build command...")

            self.logger.info("Done.")

    def _has_split_variant(self, h5ad_path: str, tag: str) -> bool:
        """Check whether split variant `tag` (e.g. 'u60') is present in uns['size_splits'] of an h5ad."""
        if not os.path.exists(h5ad_path):
            return False
        try:
            import anndata as ad
            return tag in ad.read_h5ad(h5ad_path).uns.get('size_splits', {})
        except Exception:
            return False

    def cluster_db(self) -> None:
        """Clusters the AnnData object. Split variants (added via --readlengthsplit at build
        time) now live inside the same vibrChol1.h5ad, so they're clustered in place via
        --variant rather than as separate h5ad files."""
        with self._live_box("Clustering AnnData object..."):
            cmd = (
                f"{self.trnagraph_path} analyze cluster "
                "-i vibrChol1/vibrChol1.h5ad -o vibrChol1/vibrChol1.h5ad --overwrite"
            )
            self._run_command(cmd, "Running cluster command...")

            if self._has_split_variant("vibrChol1/vibrChol1.h5ad", "u60"):
                cmd = (
                    f"{self.trnagraph_path} analyze cluster "
                    "-i vibrChol1/vibrChol1.h5ad -o vibrChol1/vibrChol1.h5ad --variant norm:u60 --overwrite"
                )
                self._run_command(cmd, "Running cluster command for under split...")

            if self._has_split_variant("vibrChol1/vibrChol1.h5ad", "o60"):
                cmd = (
                    f"{self.trnagraph_path} analyze cluster "
                    "-i vibrChol1/vibrChol1.h5ad -o vibrChol1/vibrChol1.h5ad --variant norm:o60 --overwrite"
                )
                self._run_command(cmd, "Running cluster command for over split...")

            self.logger.info("Done.")

    def graph_db(self) -> None:
        """Generates graphs from the AnnData object."""
        with self._live_box("Generating graphs..."):
            cmd = (
                f"{self.trnagraph_path} graph "
                "-i vibrChol1/vibrChol1.h5ad -o vibrChol1/graphs --colormap config/colormap.json"
            )
            self._run_command(cmd, "Running graph command...")

            self.logger.info("Done.")

    def graph_split_db(self) -> None:
        """Generates graphs for split variants. These now live inside the same vibrChol1.h5ad
        (rather than separate _u60.h5ad/_o60.h5ad files), selected in place via --variant."""
        with self._live_box("Generating graphs for split variants..."):
            cmd = (
                f"{self.trnagraph_path} graph "
                "-i vibrChol1/vibrChol1.h5ad -o vibrChol1/graphs_u60 --variant norm:u60 --colormap config/colormap.json"
            )
            self._run_command(cmd, "Running graph command for under split...")

            cmd = (
                f"{self.trnagraph_path} graph "
                "-i vibrChol1/vibrChol1.h5ad -o vibrChol1/graphs_o60 --variant norm:o60 --colormap config/colormap.json"
            )
            self._run_command(cmd, "Running graph command for over split...")

            self.logger.info("Done.")

    def main(self) -> None:
        """Main execution logic for the pipeline."""
        try:
            self.logger.info("Running tests...")
            
            specific_flags = [
                self.args.metadata, self.args.fastq, self.args.trna,
                self.args.genome, self.args.trim, self.args.makedb, self.args.map,
                self.args.build, getattr(self.args, 'split_build', False),
                self.args.cluster, self.args.merge, self.args.graph, getattr(self.args, 'split_graph', False),
                self.args.hubonly, self.args.maponly
            ]
            run_all = self.args.all or not any(specific_flags)
            # --skip-download skips metadata/fastq/trna/genome regardless of run_all, for a
            # faster dev loop against data that's already on disk (downloading is the slowest
            # part of a demo run and is often untouched between iterations). --all still forces
            # a full clean + redownload as before (--all wipes the workspace in __init__,
            # before this check ever runs, and is not itself affected by --skip-download).
            skip_download = getattr(self.args, 'skip_download', False)

            if not skip_download and (run_all or self.args.metadata):
                self.download_metadata()
            if not skip_download and (run_all or self.args.fastq):
                self.download_fastq()
            if not skip_download and (run_all or self.args.trna):
                self.download_trna()
            if not skip_download and (run_all or self.args.genome):
                self.download_genome()
            if run_all or self.args.trim:
                self.trim_fastq()
            if run_all or self.args.makedb:
                self.create_index()
            if run_all or self.args.map or self.args.maponly:
                self.map_reads()
            if (run_all and not self.args.maponly) or self.args.build or self.args.hubonly:
                self.build_db()
            # split_build runs build with readlengthsplit
            if getattr(self.args, 'split_build', False) or (run_all and not self.args.maponly):
                self.args.split_build = True  # Ensure flag is set
                self.build_db()
            if (run_all and not self.args.maponly) or self.args.cluster:
                self.cluster_db()
            if (run_all and not self.args.maponly) or self.args.graph:
                self.graph_db()
            # split_graph runs graphs for split h5ad files
            if getattr(self.args, 'split_graph', False) or (run_all and not self.args.maponly):
                self.graph_split_db()
                
            if self.args.cleanrun:
                self._cleanup_workspace()

            self.logger.info("All tests completed.")

        except Exception as e:
            self.logger.error(f"An error occurred during execution: {e}")

if __name__ == "__main__":
    pass