# Seq2VCF Pipeline

This repository contains a bioinformatics pipeline designed to process raw sequencing data (FASTQ) into high-quality Variant Call Format (VCF) files. The pipeline is fully automated, modular, and adaptable to high-throughput datasets, offering flexibility and reproducibility for various genotyping projects.

---

## Purpose and Functionality

The pipeline integrates essential steps for variant calling into a single streamlined workflow, automating tasks such as:
1. Quality Control (QC): Assesses sequencing quality.
2. Adapter Trimming: Removes adapters and filters low-quality reads.
3. Alignment: Maps reads to a reference genome.
4. Post-processing: Sorts, indexes, and cleans BAM files.
5. Variant Calling: Detects SNPs and indels, outputting filtered VCF files.
6. (Optional Bonus) Annotation: Adds functional insights to variants.

---

## Input and Output File Formats

### Input Files
- **Raw Sequencing Reads**:
  - Format: Compressed FASTQ (`*.fq.gz`)
  - Example: `data/raw_data.fq.gz`
- **Barcodes**:
  - Format: Tab-separated text file with columns: Barcode Sequence, Sample Name
  - Example: `data/barcodes.txt`
- **Reference Genome**:
  - Format: FASTA (`*.fa`)
  - Requirements: Indexed using `bwa index` and `samtools faidx`
  - Example: `refgenome/reference.fa`

### Output Structure
Pipeline results are stored in the following directories:
```bash
/results/
  ├── demultiplexed/          # Reads grouped by sample
  ├── qc_reports/             # Quality control reports
  ├── trimmed/                # Adapter-trimmed reads
  ├── alignment/              # Aligned and indexed BAM files
  ├── variants/               # Filtered VCF files
```
---

## Software Dependencies and Installation

### Dependencies

Ensure the following tools are installed:

- **Python (v3.5+)**
- **Java (OpenJDK 1.8.0_102+)**
- **Cutadapt (v2.1+)**
- **Sabre (v1.000+)**
- **BWA (v0.7.17+)**
- **SAMtools (v1.8+)**
- **BCFtools (v1.15+)**
- **FastQC (v0.11.2+)**
- **SnpEff (optional for annotation)**

Modules can be loaded via your system's environment manager.

### Installation

Clone the repository and ensure tools are accessible from your `$PATH`:
```bash
git clone https://github.com/yourgroup/Seq2VCF.git
cd Seq2VCF

---

## How to Run the Pipeline

### Prepare Input Files:
- Place raw sequencing reads in `data/`.
- Ensure the barcode file (`barcodes.txt`) is in the same directory.
- Verify the reference genome is indexed in `refgenome/`.

### Edit Parameters:
Update `parameters.conf` with paths and thresholds:
- **ADAPTER_SEQ**: Adapter sequence to trim.
- **MIN_LENGTH**: Minimum retained read length.
- **QUAL_THRESHOLD**, **DEPTH_THRESHOLD**, etc.: Variant filtering thresholds.

### Execute the Pipeline:
Run the pipeline script:
```bash
bash Seq2VCF_pipeline.sh

---
### Monitor Progress:
Use `squeue` (or equivalent) to track SLURM job status.

---

## Log Files and Troubleshooting

- Logs are generated at each step, stored in `logs/`.
- Check logs for error messages to troubleshoot issues:
  - Ensure input files are correctly formatted and indexed.
  - Verify dependencies are properly installed.

---

## Features and Flexibility

- **Reproducibility**: Relative paths ensure compatibility with any system.
- **Modularity**: Each step is parameterized and independent.
- **Flexibility**: The pipeline supports any FASTQ dataset conforming to the expected format.

---

## Bonus: Functional Annotation (Optional)

Variants can be annotated using SnpEff:
```bash
java -Xmx4g -jar snpEff.jar annotate -v GRCh37.75 variants/final_filtered.vcf.gz > variants/annotated.vcf

---

## Project Deliverables

- **Pipeline Script**: Automates all steps from raw reads to VCF generation.
- **Generated VCF**: Example output in `results/variants/`.
- **GitHub Repository**: Includes script, README.md, example files, and a log file.
- **Wiki**: Detailed documentation of each step and error-handling guidance.

---

## Authors

This pipeline was collaboratively developed by [Author 1], [Author 2], [Author 3], and [Author 4] as part of BVG-7003 Assignment 2.
