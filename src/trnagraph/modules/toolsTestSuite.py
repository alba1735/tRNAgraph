#!/usr/bin/env python3

import os
import sys
import subprocess
import logging
import zipfile
import argparse

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

        # Configure logging
        logging.basicConfig(
            filename='toolsTestSuite.log',
            filemode='w',
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)

        # Clean up if requested
        if self.args.all:
            self._cleanup_workspace()

        # Copy assets
        self.assets_dir = os.path.join(self.repo_root, "src/trnagraph/assets")
        os.makedirs("config", exist_ok=True)
        self._run_command(f"cp --update {self.assets_dir}/*.txt config/.", "Copying assets...")
        self._run_command(f"cp --update {self.assets_dir}/*.json config/.", "Copying assets...")

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
            print(description)
        
        # Ensure the current python environment's bin directory is in the PATH
        env = os.environ.copy()
        python_bin_dir = os.path.dirname(sys.executable)
        if python_bin_dir not in env.get("PATH", ""):
            env["PATH"] = f"{python_bin_dir}{os.pathsep}{env.get('PATH', '')}"

        try:
            result = subprocess.run(
                command,
                shell=True,
                check=check,
                capture_output=True,
                text=True,
                env=env
            )
            # Log output if verbose or on error (though check=True handles error)
            if result.stdout:
                self.logger.info(f"Stdout:\n{result.stdout}")
            if result.stderr:
                self.logger.info(f"Stderr:\n{result.stderr}")
            return result
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Command failed: {command}\nStdout: {e.stdout}\nStderr: {e.stderr}")
            raise e

    def _cleanup_workspace(self) -> None:
        """Removes generated files to ensure a clean run, keeping only the log file."""
        self.logger.info("Cleaning up workspace...")
        # Remove all files and directories in the current working directory (the test directory)
        # except for the log file.
        self._run_command(
            'find . -maxdepth 1 -mindepth 1 -not -name "toolsTestSuite.log" -exec rm -rf {} +',
            "Removing all contents from test directory except log file..."
        )
        self.logger.info("Workspace cleaned.")

    def download_metadata(self) -> None:
        """Downloads metadata from SRA using pysradb."""
        os.makedirs("raw/vibrChol1/fastq", exist_ok=True)
        self.logger.info("Downloading metadata from SRA...")
        print("Downloading metadata from SRA...")
        
        cmd = (
            "pysradb metadata SRP254278 | "
            "grep -v -e 'Escherichia coli' -e 'trmK' -e 'miaA' -e 'ttcA' -e 'thiI' -e 'run_accession' | "
            "cut -f22 > raw/vibrChol1/fastq/accessions.tsv"
        )
        self._run_command(cmd, "Fetching metadata...")
        self.logger.info("Done.")
        print("Done.")

    def download_fastq(self) -> None:
        """Downloads FASTQ files for the accessions listed in accessions.tsv."""
        self.logger.info("Downloading fastq files...")
        print("Downloading fastq files...")
        
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
        print("Done.")

    def download_trna(self) -> None:
        """Downloads and extracts Vibrio cholerae tRNA sequences."""
        os.makedirs("references/vibrChol1/trnas", exist_ok=True)
        self._run_command(f"cp --update {self.assets_dir}/*.gz .", "Copying assets...")
        self.logger.info("Downloading Vibrio cholerae tRNA sequences...")
        print("Downloading Vibrio cholerae tRNA sequences...")
        
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
        print("Done.")

    def download_genome(self) -> None:
        """Downloads Vibrio cholerae genome and GFF annotation, converting to GTF."""
        os.makedirs("references/vibrChol1/annotations", exist_ok=True)
        os.makedirs("references/vibrChol1/genomes", exist_ok=True)
        self.logger.info("Downloading Vibrio cholerae genome and GFF annotation...")
        print("Downloading Vibrio cholerae genome and GFF annotation...")
        
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
        print("Done.")

    def trim_fastq(self) -> None:
        """Trims adapters from FASTQ files using the tRNAgraph preprocess trim tool."""
        if os.path.exists("processed/trimmed/vibrChol1_trim_manifest_updated.txt"):
            self.logger.info("Trimmed fastq files already exist, skipping trimming.")
            print("Trimmed fastq files already exist, skipping trimming.")
            return

        self.logger.info("Trimming fastq files with fastp...")
        print("Trimming fastq files with fastp...")
        
        cmd = (
            f"{self.trnagraph_path} preprocess trim "
            "-i config/vibrChol1.manifest.txt -a1 ACTGTAGGCACCATCAATC --colormap config/colormap.json"
        )
        self._run_command(cmd, "Running trim command...")
        
        self.logger.info("Done.")
        print("Done.")

    def create_index(self) -> None:
        """Creates the Bowtie2 index for the tRNA database."""
        if os.path.exists("references/vibrChol1/trnadb/vibrChol1_db-tRNAgenome.1.bt2l"):
            self.logger.info("Bowtie2 index already exists, skipping creation.")
            print("Bowtie2 index already exists, skipping creation.")
            return

        self.logger.info("Creating bowtie2 index...")
        print("Creating bowtie2 index...")
        
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
        print("Done.")

    def map_reads(self) -> None:
        """Maps reads to the tRNA database."""
        self.logger.info("Mapping reads to tRNA genes...")
        print("Mapping reads to tRNA genes...")
        
        cmd = (
            f"{self.trnagraph_path} preprocess map "
            "-o vibrChol1 -d references/vibrChol1/trnadb/vibrChol1_db "
            "-i config/vibrChol1.metadata.txt "
            f"--bamdir processed/vibrChol1/bam"
        )
        self._run_command(cmd, "Running map command...")
        
        self.logger.info("Done.")
        print("Done.")

    def build_db(self) -> None:
        """Builds the AnnData object from the tRNAgraph output."""
        self.logger.info("Building AnnData object...")
        print("Building AnnData object...")
        
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
        print("Done.")

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
        self.logger.info("Clustering AnnData object...")
        print("Clustering AnnData object...")

        cmd = (
            f"{self.trnagraph_path} analyze cluster "
            "-i vibrChol1/vibrChol1.h5ad -o vibrChol1/vibrChol1.h5ad --overwrite"
        )
        self._run_command(cmd, "Running cluster command...")

        if self._has_split_variant("vibrChol1/vibrChol1.h5ad", "u60"):
            self.logger.info("Clustering AnnData object for under split...")
            print("Clustering AnnData object for under split...")

            cmd = (
                f"{self.trnagraph_path} analyze cluster "
                "-i vibrChol1/vibrChol1.h5ad -o vibrChol1/vibrChol1.h5ad --variant norm:u60 --overwrite"
            )
            self._run_command(cmd, "Running cluster command for under split...")

        if self._has_split_variant("vibrChol1/vibrChol1.h5ad", "o60"):
            self.logger.info("Clustering AnnData object for over split...")
            print("Clustering AnnData object for over split...")

            cmd = (
                f"{self.trnagraph_path} analyze cluster "
                "-i vibrChol1/vibrChol1.h5ad -o vibrChol1/vibrChol1.h5ad --variant norm:o60 --overwrite"
            )
            self._run_command(cmd, "Running cluster command for over split...")

        self.logger.info("Done.")
        print("Done.")

    def graph_db(self) -> None:
        """Generates graphs from the AnnData object."""
        self.logger.info("Generating graphs...")
        print("Generating graphs...")
        
        cmd = (
            f"{self.trnagraph_path} graph "
            "-i vibrChol1/vibrChol1.h5ad -o vibrChol1/graphs --colormap config/colormap.json"
        )
        self._run_command(cmd, "Running graph command...")
        
        self.logger.info("Done.")
        print("Done.")

    def graph_split_db(self) -> None:
        """Generates graphs for split variants. These now live inside the same vibrChol1.h5ad
        (rather than separate _u60.h5ad/_o60.h5ad files), selected in place via --variant."""
        self.logger.info("Generating graphs for under split...")
        print("Generating graphs for under split...")

        cmd = (
            f"{self.trnagraph_path} graph "
            "-i vibrChol1/vibrChol1.h5ad -o vibrChol1/graphs_u60 --variant norm:u60 --colormap config/colormap.json"
        )
        self._run_command(cmd, "Running graph command for under split...")

        self.logger.info("Generating graphs for over split...")
        print("Generating graphs for over split...")

        cmd = (
            f"{self.trnagraph_path} graph "
            "-i vibrChol1/vibrChol1.h5ad -o vibrChol1/graphs_o60 --variant norm:o60 --colormap config/colormap.json"
        )
        self._run_command(cmd, "Running graph command for over split...")

        self.logger.info("Done.")
        print("Done.")

    def main(self) -> None:
        """Main execution logic for the pipeline."""
        try:
            self.logger.info("Running tests...")
            print("Running tests...")
            
            specific_flags = [
                self.args.metadata, self.args.fastq, self.args.trna,
                self.args.genome, self.args.trim, self.args.makedb, self.args.map,
                self.args.build, getattr(self.args, 'split_build', False),
                self.args.cluster, self.args.merge, self.args.graph, getattr(self.args, 'split_graph', False),
                self.args.hubonly, self.args.maponly
            ]
            run_all = self.args.all or not any(specific_flags)

            if run_all or self.args.metadata:
                self.download_metadata()
            if run_all or self.args.fastq:
                self.download_fastq()
            if run_all or self.args.trna:
                self.download_trna()
            if run_all or self.args.genome:
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
                self.logger.info("Cleaning up test files...")
                print("Cleaning up test files...")
                self._cleanup_workspace()

            self.logger.info("All tests completed.")
            print("All tests completed.")

        except Exception as e:
            self.logger.error(f"An error occurred during execution: {e}")
            print(f"Error: {e}")

if __name__ == "__main__":
    pass