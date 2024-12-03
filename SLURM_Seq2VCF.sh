#!/bin/bash

# -----------------------------
# Load Parameters
# -----------------------------
source parameters.conf

# -----------------------------
# SLURM directives for resource allocation
# -----------------------------
#SBATCH --time=${SBATCH_TIME}           # Use the time allocated in parameters.conf
#SBATCH --mem=${SBATCH_MEM}             # Use the memory allocated in parameters.conf
#SBATCH --cpus-per-task=${THREADS}      # Use the number of threads allocated in parameters.conf

# -----------------------------
# Modules
# -----------------------------
module load Python/3.5
module load Java/1.8.0_102
module load Cutadapt/2.1
module load Sabre/1.000
module load BWA/v0.7.17
module load SAMtools/1.8
module load BCFtools/1.15
module load FastQC/0.11.2
module load SnpEff/5.2e

# -----------------------------
# Derived Paths
# -----------------------------
DEMUX_DIR="./results/demultiplexed"  # Directory for demultiplexed files
QC_DIR="./results/qc_reports"        # Directory for quality reports
TRIM_DIR="./results/trimmed"         # Directory for trimmed files
ALIGN_DIR="./results/alignment"      # Directory for aligned files
VCF_DIR="./results/variants"        # Directory for VCF files
LOG_DIR="./logs"            # Directory for logs

# -----------------------------
# Create required directories
# -----------------------------
mkdir -p results logs ${DEMUX_DIR} ${QC_DIR} ${TRIM_DIR} ${ALIGN_DIR} ${VCF_DIR} ${LOG_DIR}  # Create necessary directories

# -----------------------------
# Set SNP Effect Database based on REF_GENOME
# -----------------------------
SNP_EFFECT_DB=${SNP_EFFECT_DB}  # Use value deffine in parameters

# -----------------------------
# Step 1: Demultiplexing
# -----------------------------
LOG_FILE="${LOG_DIR}/demultiplexing_$(date +%Y%m%d%H%M%S).log"
echo "Starting Demultiplexing..." | tee -a ${LOG_FILE}
cd ${DATA_DIR}

# Check if barcodes file exists
if [ ! -f ${BARCODES_FILE} ]; then
    echo "Error: Barcodes file not found!" | tee -a ${LOG_FILE}
    exit 1
fi

sabre se -f *.fq.gz \
         -b ${BARCODES_FILE} \
         -u ${DEMUX_DIR}/unmatched.fq.gz \
         > ${DEMUX_DIR}/demux.log 2>&1

if [ $? -ne 0 ]; then
    echo "Error during Demultiplexing. Check logs." | tee -a ${LOG_FILE}
    exit 1
fi

cd -
echo "Demultiplexing Complete." | tee -a ${LOG_FILE}

# -----------------------------
# Step 2: Quality Control
# -----------------------------
LOG_FILE="${LOG_DIR}/qc_$(date +%Y%m%d%H%M%S).log"
echo "Starting Quality Control..." | tee -a ${LOG_FILE}
fastqc -t 10 -o ${QC_DIR} ${DATA_DIR}/*.fq  # Run FastQC for quality control

if [ $? -ne 0 ]; then
    echo "Error during Quality Control." | tee -a ${LOG_FILE}
    exit 1
fi
echo "Quality Control Complete." | tee -a ${LOG_FILE}

# -----------------------------
# Step 3: Adapter Trimming
# -----------------------------
LOG_FILE="${LOG_DIR}/adapter_trimming_$(date +%Y%m%d%H%M%S).log"
echo "Starting Adapter Trimming..." | tee -a ${LOG_FILE}

# Check if Adapter sequence is defined
if [ -z "${ADAPTER_SEQ}" ]; then
    echo "Error: Adapter sequence is not defined!" | tee -a ${LOG_FILE}
    exit 1
fi

# Parallelize trimming using GNU Parallel for efficiency
ls ${DATA_DIR}/*.fq | parallel -j ${THREADS} 'base=$(basename {} .fq); cutadapt -a ${ADAPTER_SEQ} -m ${MIN_LENGTH} -o ${TRIM_DIR}/{/.}_trimmed.fq.gz {}' 2>&1 | tee -a ${LOG_FILE}

if [ $? -ne 0 ]; then
    echo "Error during Adapter Trimming." | tee -a ${LOG_FILE}
    exit 1
fi
echo "Adapter Trimming Complete." | tee -a ${LOG_FILE}

# -----------------------------
# Step 4: Alignment
# -----------------------------
LOG_FILE="${LOG_DIR}/alignment_$(date +%Y%m%d%H%M%S).log"
echo "Starting Alignment..." | tee -a ${LOG_FILE}
for fq in ${TRIM_DIR}/*_trimmed.fq.gz; do
    base=$(basename "$fq" _trimmed.fq.gz)  # Extract base name of the file
    bwa mem -t ${THREADS} ${REF_GENOME} "$fq" > ${ALIGN_DIR}/${base}.sam  # Align with BWA

    if [ $? -ne 0 ]; then
        echo "Error in Alignment for $fq." | tee -a ${LOG_FILE}
        exit 1
    fi

    # Convert SAM to BAM, sort, and index
    samtools view -bS ${ALIGN_DIR}/${base}.sam > ${ALIGN_DIR}/${base}.bam  # Convert SAM to BAM
    samtools sort ${ALIGN_DIR}/${base}.bam -o ${ALIGN_DIR}/${base}_sorted.bam  # Sort BAM files
    samtools index ${ALIGN_DIR}/${base}_sorted.bam  # Index the sorted BAM file

    # Clean up intermediate files
    rm ${ALIGN_DIR}/${base}.sam ${ALIGN_DIR}/${base}.bam  # Remove intermediate files
done
echo "Alignment Complete." | tee -a ${LOG_FILE}

# -----------------------------
# Step 5: Variant Calling with Filters
# -----------------------------
LOG_FILE="${LOG_DIR}/variant_calling_$(date +%Y%m%d%H%M%S).log"
echo "Starting Variant Calling..." | tee -a ${LOG_FILE}
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
    echo "Error in Variant Calling. Check logs." | tee -a ${LOG_FILE}
    exit 1
fi

echo "Variant Calling Complete. Final VCF file is at: ${VCF_DIR}/variants_sorted.vcf.gz" | tee -a ${LOG_FILE}

# -----------------------------
# Step 6: Functional Annotation with SnpEff (Optional)
# -----------------------------
LOG_FILE="${LOG_DIR}/annotation_$(date +%Y%m%d%H%M%S).log"
echo "Starting Variant Annotation with SnpEff..." | tee -a ${LOG_FILE}

# Vérifie si SnpEff est installé
if ! command -v java &> /dev/null; then
    echo "Error: Java is not installed. Please install Java and try again." | tee -a ${LOG_FILE}
    exit 1
fi

# Annoter les variantes avec SnpEff en utilisant la base de données spécifiée
java -Xmx4g -jar snpEff.jar annotate -v ${SNP_EFFECT_DB} ${VCF_DIR}/variants_sorted.vcf.gz > ${VCF_DIR}/variants_annotated.vcf

if [ $? -ne 0 ]; then
    echo "Error during Variant Annotation with SnpEff." | tee -a ${LOG_FILE}
    exit 1
fi

# Compresser et indexer le VCF annoté
bgzip -c ${VCF_DIR}/variants_annotated.vcf > ${VCF_DIR}/variants_annotated.vcf.gz
tabix -p vcf ${VCF_DIR}/variants_annotated.vcf.gz

echo "Annotation Complete. Annotated VCF
