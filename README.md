# Seq2VCF Pipeline

This repository contains a bioinformatics pipeline for processing raw sequencing data into high-quality Variant Call Format (VCF) files. Designed for genotyping-by-sequencing (GBS) data, the pipeline is modular, efficient, and suitable for high-throughput batch processing.

## Pipeline Overview

The pipeline automates the following steps:
1. **Demultiplexing**: Splits sequencing reads based on barcodes.
2. **Quality Control**: Assesses sequencing quality metrics.
3. **Adapter Trimming**: Removes adapters and short reads.
4. **Alignment**: Maps reads to a reference genome.
5. **Variant Calling**: Identifies and filters genetic variants.

Each step is parameterized for reproducibility and includes robust error-handling mechanisms.

## System Requirements

### Software Dependencies
Ensure the following tools are installed and available:
- **Python (v3.5+)**
- **Java (OpenJDK 1.8.0_102+)**
- **Cutadapt (v2.1+)**
- **Sabre (v1.000+)**
- **BWA (v0.7.17+)**
- **SAMtools (v1.8+)**
- **BCFtools (v1.15+)**
- **FastQC (v0.11.2+)**

Modules can be loaded using your system's environment manager.

## Input Requirements

### Directory Structure
Organize your input files as follows:
```bash
/project/
  ├── data/
  │   ├── raw_data.fq.gz      # Raw sequencing reads
  │   ├── barcodes.txt        # Barcode file for demultiplexing
  ├── refgenome/
  │   ├── reference.fa        # Reference genome in FASTA format
  ├── results/                # Directory for outputs
```


## Input Files

### Raw Sequencing Data
- **File format**: `*.fq.gz`
- **Location**: `/data`

### Barcodes File
- **Format**: Tab-separated file with columns:
  - Barcode sequence
  - Sample name
- **Location**: `/data/barcodes.txt`

### Reference Genome
- **Format**: FASTA
- **Indexing**: Must be indexed using `bwa index` and `samtools faidx`.

## Parameters
Parameters are loaded from `parameters.conf`. Key variables include:
- **ADAPTER_SEQ**: Sequencing adapter to trim.
- **MIN_LENGTH**: Minimum length for retained reads post-trimming.
- **Variant calling thresholds**:
  - **QUAL_THRESHOLD**: Minimum variant quality score.
  - **DEPTH_THRESHOLD**: Minimum read depth.
  - **MAX_MISSING**: Maximum fraction of missing data per variant.
  - **MAF_THRESHOLD**: Minimum minor allele frequency.

## Output Structure
Pipeline results are stored under `results/`:
```bash
/results/
  ├── demultiplexed/          # Reads grouped by samples
  ├── qc_reports/             # Quality control reports
  ├── trimmed/                # Adapter-trimmed reads
  ├── alignment/              # Aligned BAM files
  ├── variants/               # Filtered VCF files
```

## Pipeline Steps

### Step 1: Demultiplexing
**Tool**: Sabre  
**Command**:
```bash
sabre se -f *.fq.gz -b barcodes.txt -u results/demultiplexed/unmatched.fq.gz > results/demultiplexed/demux.log 2>&1
```
**Output**: Reads split into sample-specific files. Logs are stored for debugging.

### Step 2: Quality Control
**Tool**: FastQC  
**Command**:
```bash
fastqc -t 10 -o results/qc_reports/ data/*.fq.gz
```
**Output**: Reports stored in `results/qc_reports/`. Key metrics include:
- Sequence quality scores.
- Adapter content.

### Step 3: Adapter Trimming
**Tool**: Cutadapt  
**Command**:
```bash
cutadapt -a $ADAPTER_SEQ -m $MIN_LENGTH -o results/trimmed/{sample}_trimmed.fq.gz {input_file}
```
**Output**: Trimmed reads stored in `results/trimmed/`.

### Step 4: Alignment
**Tool**: BWA and SAMtools  
**Commands**:
- **Align reads**:
  ```bash
  bwa mem -t 10 refgenome/reference.fa {input_file} > alignment/{sample}.sam
  ```
- **Convert and sort BAM**:
  ```bash
  samtools view -bS alignment/{sample}.sam | samtools sort > alignment/{sample}_sorted.bam
  samtools index alignment/{sample}_sorted.bam
  ```
**Output**: Indexed BAM files stored in `results/alignment/`.

### Step 5: Variant Calling
**Tool**: BCFtools  
**Commands**:
- **Generate raw variants**:
  ```bash
  samtools mpileup -g -f refgenome/reference.fa -b bam_list.txt | bcftools call -mv -Ov > variants/raw.vcf
  ```
- **Filter variants**:
  ```bash
  bcftools filter -e "QUAL<$QUAL_THRESHOLD || DP<$DEPTH_THRESHOLD" variants/raw.vcf > variants/filtered.vcf
  bcftools filter -e "F_MISSING>$MAX_MISSING || MAF<$MAF_THRESHOLD" variants/filtered.vcf > variants/final_filtered.vcf
  ```
- **Compress and index final VCF**:
  ```bash
  bgzip -c variants/final_filtered.vcf > variants/final_filtered.vcf.gz
  tabix -p vcf variants/final_filtered.vcf.gz
  ```
**Output**: Final filtered VCF stored in `results/variants/final_filtered.vcf.gz`.

## How to Run
1. **Set Parameters**: Update paths and thresholds in `parameters.conf`.
2. **Submit Job**: Run the pipeline via SLURM:
   ```bash
   sbatch Seq2VCF_pipeline.sh
   ```
3. **Monitor Progress**: Check job status using:
   ```bash
   squeue
   ```

## Troubleshooting
- Logs are generated for each step to identify and debug errors.
- Ensure all input files are correctly formatted and indexed.

