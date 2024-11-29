# BVG-7003-Devoir-2--Seq2VCF


## Seq2VCF Pipeline

This repository provides a bash-based pipeline for genotyping-by-sequencing (GBS) data analysis. It includes all steps from demultiplexing to variant calling, with filtering and quality control. The script is modular and designed for high-throughput batch processing.

## Dependencies

Before running the pipeline, ensure the following tools and versions are installed and available in your system:

- **Python**: v3.5+
- **Java**: OpenJDK 1.8.0_102+
- **Cutadapt**: v2.1+
- **Sabre**: v1.000+
- **BWA**: v0.7.17+
- **Samtools**: v1.8+
- **Bcftools**: v1.15+
- **FastQC**: v0.11.2+

Use the provided module commands in the script or load the required modules as per your system's module manager.

## Input Files

Organize your input files in the following structure:

```bash
/project/
  ├── data/
  │   ├── <raw_data.fq.gz>      # Raw sequencing reads
  │   ├── barcodes.txt          # Barcode file for demultiplexing
  ├── refgenome/
  │   ├── reference.fa          # Reference genome in FASTA format
```

### Input Details

- **Raw Data**:
  - File format: `*.fq.gz`
  - Location: `/data`
- **Barcodes**:
  - Tab-separated file where each line contains:
    - Barcode sequence
    - Sample name
- **Reference Genome**:
  - FASTA format
  - Ensure the reference is indexed using `bwa index` and `samtools faidx`.

## Output Files

The pipeline will generate the following directories and files under the `/results` directory:

```bash
/results/
  ├── demultiplexed/          # Demultiplexed reads (matched and unmatched)
  ├── qc_reports/             # Quality control reports
  ├── trimmed/                # Adapter-trimmed reads
  ├── alignment/              # Sorted BAM files with updated headers
  ├── variants/               # VCF files (raw, filtered, and sorted)
```

## How to Run

### Set up User Variables

Modify the following variables in the script to reflect your paths:

```bash
DATA_DIR="/path/to/data"
BARCODES_FILE="${DATA_DIR}/barcodes.txt"
REF_GENOME="/path/to/reference.fa"
RESULTS_DIR="/path/to/results"
```

### Submit the Job

Run the script using a scheduler like SLURM:

```bash
sbatch CanGBS.sh
```

### Output

The final VCF file will be available at:

```bash
${RESULTS_DIR}/variants/variants_sorted.vcf.gz
```
```