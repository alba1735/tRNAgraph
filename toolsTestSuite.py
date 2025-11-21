#!/usr/bin/env python3

import os
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
        self.repo_root = os.path.abspath(os.path.dirname(__file__))
        
        # Set up working directory
        os.chdir(self.repo_root)
        os.makedirs("tests/", exist_ok=True)
        os.chdir(os.path.join(self.repo_root, "tests/"))

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
        self._run_command("cp --update=none ../assets/* .", "Copying assets...")

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
        
        try:
            result = subprocess.run(
                command,
                shell=True,
                check=check,
                capture_output=True,
                text=True
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
        """Removes generated files to ensure a clean run."""
        dirs_to_remove = "fastq_raw fastq_trim trnadb vibrChol1 vibrChol1-tRNAs vibrChol1-tRAX"
        files_to_remove = "accessions.tsv mismatchcompare.txt vibrChol1-tRNAs.tar.gz vibrChol1.*.txt VC_*.bam VC_*.bai"
        
        self._run_command(f"rm -rf {dirs_to_remove}", "Cleaning directories...")
        self._run_command(f"rm -f {files_to_remove}", "Cleaning files...")

    def download_metadata(self) -> None:
        """Downloads metadata from SRA using pysradb."""
        self.logger.info("Downloading metadata from SRA...")
        print("Downloading metadata from SRA...")
        
        cmd = (
            "pysradb metadata SRP254278 | "
            "grep -v -e 'Escherichia coli' -e 'trmK' -e 'miaA' -e 'ttcA' -e 'thiI' -e 'run_accession' | "
            "cut -f22 > accessions.tsv"
        )
        self._run_command(cmd, "Fetching metadata...")
        self.logger.info("Done.")
        print("Done.")

    def download_fastq(self) -> None:
        """Downloads FASTQ files for the accessions listed in accessions.tsv."""
        self.logger.info("Downloading fastq files...")
        print("Downloading fastq files...")
        
        try:
            with open("accessions.tsv", "r") as f:
                accessions = f.read().splitlines()
            
            for acc in accessions:
                if os.path.exists(f"fastq_raw/{acc}.fastq"):
                    self.logger.info(f"{acc}.fastq already exists, skipping download.")
                    continue
                
                # Prefetch
                self._run_command(
                    f"prefetch {acc} --output-file fastq_raw/{acc}.sra",
                    f"Prefetching {acc}..."
                )
                
                # Fastq-dump
                self._run_command(
                    f"fastq-dump -O fastq_raw -X 100000 {acc}",
                    f"Dumping fastq for {acc}..."
                )
                
                # Cleanup SRA file
                sra_file = f"fastq_raw/{acc}.sra"
                if os.path.exists(sra_file):
                    os.remove(sra_file)
                    
        except FileNotFoundError:
            self.logger.error("accessions.tsv not found. Run download_metadata first.")
            raise

        self.logger.info("Done.")
        print("Done.")

    def download_trna(self) -> None:
        """Downloads and extracts Vibrio cholerae tRNA sequences."""
        os.makedirs("vibrChol1-tRNAs", exist_ok=True)
        self.logger.info("Downloading Vibrio cholerae tRNA sequences...")
        print("Downloading Vibrio cholerae tRNA sequences...")
        
        if os.path.exists("vibrChol1-tRNAs/vibrChol1-tRNAs.fa"):
            self.logger.info("vibrChol1-tRNAs.fa already exists, skipping download.")
        else:
            # Using local tar.gz as per original logic - gtRNAdb has issues currently with downloads
            self._run_command(
                "tar xzvf vibrChol1-tRNAs.tar.gz -C vibrChol1-tRNAs/",
                "Extracting tRNA sequences..."
            )
            
        self.logger.info("Done.")
        print("Done.")

    def download_genome(self) -> None:
        """Downloads Vibrio cholerae genome and GFF annotation, converting to GTF."""
        self.logger.info("Downloading Vibrio cholerae genome and GFF annotation...")
        print("Downloading Vibrio cholerae genome and GFF annotation...")
        
        os.makedirs("vibrChol1", exist_ok=True)
        if os.path.exists("vibrChol1/genes.gtf"):
            self.logger.info("Vibrio cholerae genome already exists, skipping download.")
            return

        # Download genome
        download_cmd = (
            'curl -s -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/'
            'GCF_000006745.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,'
            'CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000006745.1.zip" '
            '-H "Accept: application/zip"'
        )
        self._run_command(download_cmd, "Downloading genome zip...")

        # Extract
        with zipfile.ZipFile("GCF_000006745.1.zip", "r") as zip_ref:
            zip_ref.extractall("genomes")
        os.remove("GCF_000006745.1.zip")

        # Modify FASTA headers
        base_path = "genomes/ncbi_dataset/data/GCF_000006745.1"
        fna_path = f"{base_path}/GCF_000006745.1_ASM674v1_genomic.fna"
        
        sed_cmd = (
            f'sed -i -e "/NC_002505.1/c\\>chrI" {fna_path} && '
            f'sed -i -e "/NC_002506.1/c\\>chrII" {fna_path}'
        )
        self._run_command(sed_cmd, "Modifying FASTA headers...")

        # Convert GFF to GTF
        gff_path = f"{base_path}/genomic.gff"
        gtf_path = f"{base_path}/genomic.gtf"
        self._run_command(f"gffread -E {gff_path} -T -o {gtf_path}", "Converting GFF to GTF...")

        # Modify GTF and filter
        final_gtf_cmd = (
            f"cat {gtf_path} | sed 's/NC_002505.1/chrI/g' | "
            f"sed 's/NC_002506.1/chrII/g' | grep -v '^#' > {base_path}/genes.gtf"
        )
        self._run_command(final_gtf_cmd, "Finalizing GTF file...")

        # Move files
        os.rename(fna_path, "vibrChol1/GCF_000006745.1_ASM674v1_genomic.fna")
        os.rename(f"{base_path}/genes.gtf", "vibrChol1/genes.gtf")
        
        # Cleanup
        self._run_command("rm -rf genomes", "Removing temporary genomes directory...")
        
        self.logger.info("Done.")
        print("Done.")

    def trim_fastq(self) -> None:
        """Trims adapters from FASTQ files using the tRNAgraph preprocess trim tool."""
        if os.path.exists("fastq_trim/vibrChol1_trim_manifest_updated.txt"):
            self.logger.info("Trimmed fastq files already exist, skipping trimming.")
            print("Trimmed fastq files already exist, skipping trimming.")
            return

        self.logger.info("Trimming fastq files with fastp...")
        print("Trimming fastq files with fastp...")
        
        cmd = (
            "python3 ../trnagraph.py preprocess trim "
            "-r vibrChol1 -i vibrChol1.manifest.txt -a1 ACTGTAGGCACCATCAATC"
        )
        self._run_command(cmd, "Running trim command...")
        
        self.logger.info("Done.")
        print("Done.")

    def create_index(self) -> None:
        """Creates the Bowtie2 index for the tRNA database."""
        if os.path.exists("trnadb/vibrChol1_db-tRNAgenome.1.bt2l"):
            self.logger.info("Bowtie2 index already exists, skipping creation.")
            print("Bowtie2 index already exists, skipping creation.")
            return

        self.logger.info("Creating bowtie2 index...")
        print("Creating bowtie2 index...")
        
        cmd = (
            "python3 ../trnagraph.py preprocess makedb "
            "-g vibrChol1/GCF_000006745.1_ASM674v1_genomic.fna "
            "-t vibrChol1-tRNAs/vibrChol1-tRNAs.out "
            "-r vibrChol1-tRNAs/vibrChol1-tRNAs.fa "
            "-m vibrChol1-tRNAs/vibrChol1-tRNAs_name_map.txt "
            "-s bact -o trnadb/vibrChol1_db"
        )
        self._run_command(cmd, "Running makedb command...")
        
        self.logger.info("Done.")
        print("Done.")

    def map_reads(self) -> None:
        """Maps reads to the tRNA database."""
        self.logger.info("Mapping reads to tRNA genes...")
        print("Mapping reads to tRNA genes...")
        
        cmd = (
            "python3 ../trnagraph.py preprocess map "
            "-e vibrChol1-tRAX -d trnadb/vibrChol1_db "
            "-s vibrChol1.metadata.txt --pairs vibrChol1.pair.txt "
            "--gtf vibrChol1/genes.gtf --traxmode"
        )
        self._run_command(cmd, "Running map command...")
        
        self.logger.info("Done.")
        print("Done.")

    def main(self) -> None:
        """Main execution logic for the pipeline."""
        try:
            self.logger.info("Running tests...")
            print("Running tests...")
            
            specific_flags = [
                self.args.metadata, self.args.fastq, self.args.trna,
                self.args.genome, self.args.trim, self.args.makedb, self.args.map
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
            if run_all or self.args.map:
                self.map_reads()
                
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