#!/usr/bin/env python3

import os
import subprocess
import logging
import zipfile

repo_root = os.path.abspath(os.path.dirname(__file__))
os.chdir(repo_root)
os.makedirs(os.path.dirname("tests/"), exist_ok=True)
os.chdir(repo_root+"/tests/")

logging.basicConfig(filename='toolsTestSuite.log', filemode='w', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def log_subprocess_output(process):
    stdout, stderr = process.communicate()
    if stdout:
        logger.info("Subprocess output:\n%s", stdout.decode())
    if stderr:
        logger.error("Subprocess error output:\n%s", stderr.decode())

try:
    logger.info("Running tests...")
    print("Running tests...")
    # Download metadata
    logger.info("Downloading metadata from SRA...")
    print("Downloading metadata from SRA...")
    subprocess.run("pysradb metadata SRP254278 | grep -v -e 'Escherichia coli' -e 'trmK' -e 'miaA' -e 'ttcA' -e 'thiI' -e 'run_accession' | cut -f22 > accessions.tsv", shell=True, check=True)
    logger.info("Done.")
    print("Done.")
    # makdir output directory for fastq files
    logger.info("Downloading fastq files...")
    print("Downloading fastq files...")
    with open("accessions.tsv", "r") as f:
        accessions = f.read().splitlines()
        # Download fastq files using prefetch and fastq-dump
        for i in accessions:
            if os.path.exists(f"fastq_raw/{i}.fastq"):
                logger.info(f"{i}.fastq already exists, skipping download.")
            else:
                result = subprocess.run(f"prefetch {i} --output-file fastq_raw/{i}.sra", shell=True, check=True, capture_output=True, text=True)
                logger.info(f"Prefetch output for {i}:\n{result.stdout}\n{result.stderr}")
                result = subprocess.run(f"fastq-dump -O fastq_raw -X 100000 {i}", shell=True, check=True, capture_output=True, text=True)
                logger.info(f"fastq-dump output for {i}:\n{result.stdout}\n{result.stderr}")
            # Remove the .sra file to save space
            sra_file = f"fastq_raw/{i}.sra"
            if os.path.exists(sra_file):
                os.remove(sra_file)
    logger.info("Done.")
    print("Done.")
    # Download Vibrio cholerae tRNA sequences
    os.makedirs("vibrChol1-tRNAs", exist_ok=True)
    logger.info("Downloading Vibrio cholerae tRNA sequences...")
    print("Downloading Vibrio cholerae tRNA sequences...")
    if os.path.exists("vibrChol1-tRNAs/vibrChol1-tRNAs.fa"):
        logger.info("vibrChol1-tRNAs.fa already exists, skipping download.")
    else:
        # gtrnadb is currently using redirects that are difficult to handle with subprocess, so we are just including the tar.gz file directly in the repo for now
        # result = subprocess.run("curl -o vibrChol1-tRNAs.tar.gz -L https://gtrnadb.org/genomes/bacteria/Vibr_chol_O1_biovar_El_Tor_N16961/vibrChol1-tRNAs.tar.gz && tar xzvf vibrChol1-tRNAs.tar.gz -C vibrChol1-tRNAs/", shell=True, check=True, capture_output=True, text=True)
        result = subprocess.run("tar xzvf vibrChol1-tRNAs.tar.gz -C vibrChol1-tRNAs/", shell=True, check=True, capture_output=True, text=True)
        logger.info(f"Download output:\n{result.stdout}\n{result.stderr}")
    logger.info("Done.")
    print("Done.")
    # Download the Genome and GFF annotation for Vibrio cholerae
    logger.info("Downloading Vibrio cholerae genome and GFF annotation...")
    print("Downloading Vibrio cholerae genome and GFF annotation...")
    os.makedirs("vibrChol1", exist_ok=True)
    if os.path.exists("vibrChol1/genes.gtf"):
        logger.info("Vibrio cholerae genome already exists, skipping download.")
    else:
        result = subprocess.run('''curl -s -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000006745.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000006745.1.zip" -H "Accept: application/zip"''', shell=True, check=True, capture_output=True, text=True)
        logger.info(f"Download output:{result.stdout}{result.stderr}")
        with zipfile.ZipFile("GCF_000006745.1.zip", "r") as zip_ref:
            zip_ref.extractall("genomes")
        os.remove("GCF_000006745.1.zip")
        # Modify the fasta and gff files to have chrI and chrII instead of NC_002505.1 and NC_002506.1
        result = subprocess.run('''sed -i -e \"/NC_002505.1/c\\>chrI\" genomes/ncbi_dataset/data/GCF_000006745.1/GCF_000006745.1_ASM674v1_genomic.fna && sed -i -e \"/NC_002506.1/c\\>chrII\" genomes/ncbi_dataset/data/GCF_000006745.1/GCF_000006745.1_ASM674v1_genomic.fna''', shell=True, check=True, capture_output=True, text=True)
        logger.info(f"sed output:\n{result.stdout}\n{result.stderr}")
        # Create GTF file from GFF
        result = subprocess.run('''gffread -E genomes/ncbi_dataset/data/GCF_000006745.1/genomic.gff -T -o genomes/ncbi_dataset/data/GCF_000006745.1/genomic.gtf''', shell=True, check=True, capture_output=True, text=True)
        logger.info(f"gffread output:\n{result.stdout}\n{result.stderr}")
        # Modify the GTF file to have chrI and chrII instead of NC_002505.1 and NC_002506.1
        result = subprocess.run('''cat genomes/ncbi_dataset/data/GCF_000006745.1/genomic.gtf | sed 's/NC_002505.1/chrI/g' | sed 's/NC_002506.1/chrII/g' | grep -v '^#' > genomes/ncbi_dataset/data/GCF_000006745.1/genes.gtf''', shell=True, check=True, capture_output=True, text=True)
        logger.info(f"gene creation output:\n{result.stdout}\n{result.stderr}")
        # Move the files .fna and .gtf to vibrChol1 directory
        os.rename("genomes/ncbi_dataset/data/GCF_000006745.1/GCF_000006745.1_ASM674v1_genomic.fna", "vibrChol1/GCF_000006745.1_ASM674v1_genomic.fna")
        os.rename("genomes/ncbi_dataset/data/GCF_000006745.1/genes.gtf", "vibrChol1/genes.gtf")
        # Remove the genomes directory to save space
        subprocess.run("rm -rf genomes", shell=True, check=True)
    logger.info("Done.")
    print("Done.")
    # Trim the fastq files with fastp
    if not os.path.exists("fastq_trim/vibrChol1_trim_manifest_updated.txt"):
        logger.info("Trimming fastq files with fastp...")
        print("Trimming fastq files with fastp...")
        result = subprocess.run("python3 ../trnagraph.py preprocess trim -r vibrChol1 -i vibrChol1.manifest.txt -a1 ACTGTAGGCACCATCAATC", shell=True, check=True, capture_output=True, text=True)
        logger.info(f"Fastp trimming output:\n{result.stdout}\n{result.stderr}")
        logger.info("Done.")
        print("Done.")
    else:
        logger.info("Trimmed fastq files already exist, skipping trimming.")
        print("Trimmed fastq files already exist, skipping trimming.")
    # Create bowtie2 index
    if not os.path.exists("trnadb/vibrChol1_db-tRNAgenome.1.bt2l"):
        logger.info("Creating bowtie2 index...")
        print("Creating bowtie2 index...")
        result = subprocess.run("python3 ../trnagraph.py preprocess makedb -g vibrChol1/GCF_000006745.1_ASM674v1_genomic.fna -t vibrChol1-tRNAs/vibrChol1-tRNAs.out -r vibrChol1-tRNAs/vibrChol1-tRNAs.fa -m vibrChol1-tRNAs/vibrChol1-tRNAs_name_map.txt -s bact -o trnadb/vibrChol1_db", shell=True, check=True, capture_output=True, text=True)
        logger.info(f"Bowtie2 index creation output:\n{result.stdout}\n{result.stderr}")
        logger.info("Done.")
        print("Done.")
    else:
        logger.info("Bowtie2 index already exists, skipping creation.")
        print("Bowtie2 index already exists, skipping creation.")
    # Map reads to tRNA genes
    logger.info("Mapping reads to tRNA genes...")
    print("Mapping reads to tRNA genes...")
    result = subprocess.run("python3 ../trnagraph.py preprocess map -e vibrChol1-map -d trnadb/vibrChol1_db -s vibrChol1.metadata.txt --hub", shell=True, check=True, capture_output=True, text=True)
    logger.info(f"Mapping output:\n{result.stdout}\n{result.stderr}")
    logger.info("Done.")
    print("Done.")

    logger.info("All tests completed.")
    print("All tests completed.")

except subprocess.CalledProcessError as e:
    logger.error("Subprocess failed with exit code %d:\nStdout: %s\nStderr: %s", e.returncode, e.stdout, e.stderr)

except FileNotFoundError as e:
    logger.error("File not found: %s", e)