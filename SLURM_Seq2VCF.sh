#!/bin/bash

# -----------------------------
# Load Parameters
# -----------------------------
source parameters.conf

# Derived Paths
DEMUX_DIR="${RESULTS_DIR}/demultiplexed"  # Directory for demultiplexed files
QC_DIR="${RESULTS_DIR}/qc_reports"          # Directory for quality reports
TRIM_DIR="${RESULTS_DIR}/trimmed"           # Directory for trimmed files
ALIGN_DIR="${RESULTS_DIR}/alignment"         # Directory for aligned files
VCF_DIR="${RESULTS_DIR}/variants"            # Directory for VCF files

# Create required directories
mkdir -p ${DEMUX_DIR} ${QC_DIR} ${TRIM_DIR} ${ALIGN_DIR} ${VCF_DIR}  # Create necessary directories

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
fastqc -t 10 -o ${QC_DIR} ${DATA_DIR}/*.fq  # Run FastQC for quality control

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
    base=$(basename "$fq" .fq)  # Extract base name of the file
    cutadapt -a ${ADAPTER_SEQ} -m ${MIN_LENGTH} -o ${TRIM_DIR}/${base}_trimmed.fq.gz "$fq"  # Trim adapters
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
    base=$(basename "$fq" _trimmed.fq.gz)  # Extract base name of the file
    bwa mem -t 10 ${REF_GENOME} "$fq" > ${ALIGN_DIR}/${base}.sam  # Align with BWA
    if [ $? -ne 0 ]; then
        echo "Error in Alignment for $fq."
        exit 1
    fi

    # Convert SAM to BAM, sort, and index
    samtools view -bS ${ALIGN_DIR}/${base}.sam > ${ALIGN_DIR}/${base}.bam  # Convert SAM to BAM
    samtools sort ${ALIGN_DIR}/${base}.bam -o ${ALIGN_DIR}/${base}_sorted.bam  # Sort BAM files
    samtools index ${ALIGN_DIR}/${base}_sorted.bam  # Index the sorted BAM file

    # Clean up intermediate files
    rm ${ALIGN_DIR}/${base}.sam ${ALIGN_DIR}/${base}.bam  # Remove intermediate files
done
echo "Alignment Complete."

# -----------------------------
# Step 5: Variant Calling with Filters
# -----------------------------
echo "Starting Variant Calling..."
ls ${ALIGN_DIR}/*_sorted.bam > ${VCF_DIR}/bam_list.txt  # Create a list of BAM files

samtools mpileup -g -f ${REF_GENOME} -b ${VCF_DIR}/bam_list.txt | \
bcftools call -mv -Ov -o ${VCF_DIR}/variants_raw.vcf  # Call variants

# Apply filters
bcftools filter -e "QUAL<${QUAL_THRESHOLD} || DP<${DEPTH_THRESHOLD}" ${VCF_DIR}/variants_raw.vcf -o ${VCF_DIR}/variants_filtered.vcf  # Apply filters
bcftools filter -e "F_MISSING > ${MAX_MISSING} || MAF < ${MAF_THRESHOLD}" ${VCF_DIR}/variants_filtered.vcf -o ${VCF_DIR}/variants_filtered_final.vcf  # Apply additional filters
bcftools sort ${VCF_DIR}/variants_filtered_final.vcf -o ${VCF_DIR}/variants_sorted.vcf  # Sort variants

bgzip -c ${VCF_DIR}/variants_sorted.vcf > ${VCF_DIR}/variants_sorted.vcf.gz  # Compress the sorted VCF file
tabix -p vcf ${VCF_DIR}/variants_sorted.vcf.gz  # Index the compressed VCF file

if [ $? -ne 0 ]; then
    echo "Error in Variant Calling. Check logs."
    exit 1
fi

echo "Variant Calling Complete. Final VCF file is at: ${VCF_DIR}/variants_sorted.vcf.gz"

# -----------------------------
# Step 6: Functional Annotation with SnpEff (Optional)
# -----------------------------
echo "Starting Variant Annotation with SnpEff..."
java -Xmx4g -jar snpEff.jar annotate -v ${REF_GENOME} ${VCF_DIR}/variants_sorted.vcf.gz > ${VCF_DIR}/variants_annotated.vcf

if [ $? -ne 0 ]; then
    echo "Error during Variant Annotation with SnpEff."
    exit 1
fi

# Compress and index the annotated VCF
bgzip -c ${VCF_DIR}/variants_annotated.vcf > ${VCF_DIR}/variants_annotated.vcf.gz
tabix -p vcf ${VCF_DIR}/variants_annotated.vcf.gz

echo "Annotation Complete. Annotated VCF file is at: ${VCF_DIR}/variants_annotated.vcf.gz"
