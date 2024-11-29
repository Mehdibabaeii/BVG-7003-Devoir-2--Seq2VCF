#!/bin/bash


# -----------------------------
# User-defined Variables
# -----------------------------
DATA_DIR="/path/to/data"                     # Path to the raw data directory
BARCODES_FILE="${DATA_DIR}/barcodes.txt"     # Path to the barcode file
REF_GENOME="/path/to/reference.fa"           # Path to the reference genome
RESULTS_DIR="/path/to/results"               # Path to the results directory

# Derived Paths
DEMUX_DIR="${RESULTS_DIR}/demultiplexed"
QC_DIR="${RESULTS_DIR}/qc_reports"
TRIM_DIR="${RESULTS_DIR}/trimmed"
ALIGN_DIR="${RESULTS_DIR}/alignment"
VCF_DIR="${RESULTS_DIR}/variants"

# Create required directories
mkdir -p ${DEMUX_DIR} ${QC_DIR} ${TRIM_DIR} ${ALIGN_DIR} ${VCF_DIR}

# -----------------------------
# Step 1: Demultiplexing
# -----------------------------
echo "Starting Demultiplexing..."
cd ${DATA_DIR}

sabre se -f *.fq.gz \
         -b ${BARCODES_FILE} \
         -u ${DEMUX_DIR}/unmatched.fq.gz \
         > ${DEMUX_DIR}/demux.log 2>&1

if [ $? -ne 0 ]; then
    echo "Error during Demultiplexing. Check logs."
    exit 1
fi

cd -
echo "Demultiplexing Complete."

# Verify demultiplexed files
if ls ${DATA_DIR}/*.fq 1> /dev/null 2>&1; then
    echo "Demultiplexed files found."
else
    echo "No demultiplexed files found. Exiting."
    exit 1
fi

# -----------------------------
# Step 2: Quality Control
# -----------------------------
echo "Starting Quality Control..."
fastqc -t 10 -o ${QC_DIR} ${DATA_DIR}/*.fq

if [ $? -ne 0 ]; then
    echo "Error during Quality Control."
    exit 1
fi
echo "Quality Control Complete."

# -----------------------------
# Step 3: Adapter Trimming
# -----------------------------
echo "Starting Adapter Trimming..."
for fq in ${DATA_DIR}/*.fq; do
    base=$(basename "$fq" .fq)
    cutadapt -a AGATCGGAAGAGCGGG -m 50 -o ${TRIM_DIR}/${base}_trimmed.fq.gz "$fq"
    if [ $? -ne 0 ]; then
        echo "Error during Adapter Trimming for $fq."
        exit 1
    fi
done
echo "Adapter Trimming Complete."

# -----------------------------
# Step 4: Alignment
# -----------------------------
echo "Starting Alignment..."
for fq in ${TRIM_DIR}/*_trimmed.fq.gz; do
    base=$(basename "$fq" _trimmed.fq.gz)
    bwa mem -t 10 ${REF_GENOME} "$fq" > ${ALIGN_DIR}/${base}.sam
    if [ $? -ne 0 ]; then
        echo "Error in Alignment for $fq."
        exit 1
    fi

    # Convert SAM to BAM, sort, and index
    samtools view -bS ${ALIGN_DIR}/${base}.sam > ${ALIGN_DIR}/${base}.bam
    samtools sort ${ALIGN_DIR}/${base}.bam -o ${ALIGN_DIR}/${base}_sorted.bam
    samtools index ${ALIGN_DIR}/${base}_sorted.bam

    # Clean up intermediate files
    rm ${ALIGN_DIR}/${base}.sam ${ALIGN_DIR}/${base}.bam
done
echo "Alignment Complete."

# -----------------------------
# Step 5: Variant Calling with Filters
# -----------------------------
echo "Starting Variant Calling..."
ls ${ALIGN_DIR}/*_sorted.bam > ${VCF_DIR}/bam_list.txt

samtools mpileup -g -f ${REF_GENOME} -b ${VCF_DIR}/bam_list.txt | \
bcftools call -mv -Ov -o ${VCF_DIR}/variants_raw.vcf

# Apply filters
bcftools filter -e 'QUAL<20 || DP<10' ${VCF_DIR}/variants_raw.vcf -o ${VCF_DIR}/variants_filtered.vcf
#bcftools filter -e 'F_MISSING > 0.2 || MAF < 0.01' ${VCF_DIR}/variants_filtered.vcf -o ${VCF_DIR}/variants_filtered_final.vcf
bcftools sort ${VCF_DIR}/variants_filtered.vcf -o ${VCF_DIR}/variants_sorted.vcf

bgzip -c ${VCF_DIR}/variants_sorted.vcf > ${VCF_DIR}/variants_sorted.vcf.gz
tabix -p vcf ${VCF_DIR}/variants_sorted.vcf.gz

if [ $? -ne 0 ]; then
    echo "Error in Variant Calling. Check logs."
    exit 1
fi

echo "Variant Calling Complete. Final VCF file is at: ${VCF_DIR}/variants_sorted.vcf.gz"