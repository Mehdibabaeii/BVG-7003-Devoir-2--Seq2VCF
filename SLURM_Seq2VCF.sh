#!/bin/bash

# -----------------------------
# Load Modules
# -----------------------------
module load python/3.5
module load cutadapt/2.1
module load sabre/1.000
module load bwa/0.7.17
module load samtools/1.8
module load bcftools/1.15
module load fastqc/0.11.2
module load snpEff/5.2e
module load java/jdk/22.0.2


# -----------------------------
# Load Parameters
# -----------------------------
source parameters.conf

# Derived Paths
DEMUX_DIR="${RESULTS_DIR}/demultiplexed"
QC_DIR="${RESULTS_DIR}/qc_reports"
TRIM_DIR="${RESULTS_DIR}/trimmed"
ALIGN_DIR="${RESULTS_DIR}/alignment"
VCF_DIR="${RESULTS_DIR}/variants"
ANOT_DIR="${RESULTS_DIR}/annoted_variants"

# Create required directories
mkdir -p ${DEMUX_DIR} ${QC_DIR} ${TRIM_DIR} ${ALIGN_DIR} ${VCF_DIR} ${ANOT_DIR}

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
    cutadapt -a ${ADAPTER_SEQ} -m ${MIN_LENGTH} -o ${TRIM_DIR}/${base}_trimmed.fq.gz "$fq"
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
bcftools filter -e "QUAL<${QUAL_THRESHOLD} || DP<${DEPTH_THRESHOLD}" ${VCF_DIR}/variants_raw.vcf -o ${VCF_DIR}/variants_filtered.vcf
bcftools filter -e "F_MISSING > ${MAX_MISSING} || MAF < ${MAF_THRESHOLD}" ${VCF_DIR}/variants_filtered.vcf -o ${VCF_DIR}/variants_filtered_final.vcf
bcftools sort ${VCF_DIR}/variants_filtered_final.vcf -o ${VCF_DIR}/variants_sorted.vcf

bgzip -c ${VCF_DIR}/variants_sorted.vcf > ${VCF_DIR}/variants_sorted.vcf.gz
tabix -p vcf ${VCF_DIR}/variants_sorted.vcf.gz

if [ $? -ne 0 ]; then
    echo "Error in Variant Calling. Check logs."
    exit 1
fi

echo "Variant Calling Complete. Final VCF file is at: ${VCF_DIR}/variants_sorted.vcf.gz"

# -----------------------------
# Step 6: Functional Annotation with SnpEff (Optional)
# -----------------------------

if [ "${RUN_ANNOTATION}" = "true" ]; then
    echo "Annotation is enabled. Starting SnpEff annotation..."

    # Vérifiez si snpEff est disponible
    if [ ! -f ${SNP_EFF_PATH} ]; then
        echo "Error: snpEff.jar not found at ${SNP_EFF_PATH}!"
        exit 1
    fi

    # Vérifiez si le génome de référence est disponible dans le répertoire de données SnpEff
    if [ ! -d "${HOME}/.snpEff/data/${SNP_EFFECT_DB}" ]; then
        echo "Genome ${SNP_EFFECT_DB} not found. Attempting to download..."

        # Télécharger le génome en fonction du nom dans SNP_EFFECT_DB
        java -jar ${SNP_EFF_PATH} download ${SNP_EFFECT_DB}
        if [ $? -ne 0 ]; then
            echo "Error: Could not download genome ${SNP_EFFECT_DB}. Check your network or database name."
            exit 1
        fi
        echo "Genome ${SNP_EFFECT_DB} successfully downloaded."
    else
        echo "Genome ${SNP_EFFECT_DB} is already available. Proceeding with annotation."
    fi

    # Exécutez l'annotation
    java -jar ${SNP_EFF_PATH} ann ${SNP_EFFECT_DB} ${VCF_DIR}/variants_sorted.vcf.gz > ${ANOT_DIR}/variants_annotated.vcf
    bgzip -c ${ANOT_DIR}/variants_annotated.vcf > ${ANOT_DIR}/variants_annotated.vcf.gz
    tabix -p vcf ${ANOT_DIR}/variants_annotated.vcf.gz

    if [ $? -ne 0 ]; then
        echo "Error during SnpEff annotation. Check logs."
        exit 1
    fi
    echo "Annotation complete. Annotated variants stored in: ${ANOT_DIR}/variants_annotated.vcf.gz"
else
    echo "Annotation is disabled. Skipping SnpEff annotation."
fi