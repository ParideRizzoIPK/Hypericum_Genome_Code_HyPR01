#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --partition=cpu        
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16     # Queue-optimized: 16 cores find slots significantly faster
#SBATCH --mem=100G             # Mathematically safe buffer (35G Index + 60G Sort Cap)
#SBATCH --time=2-00:00:00        
#SBATCH --job-name=STAR_competitive_aligner
#SBATCH --output=logs/star_alignments%j.out
#SBATCH --error=logs/star_alignments%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  # Defer strict -u validation until after custom user profiles are fully loaded

# ==============================================================================
#                       ENVIRONMENT MODULE LOADING
# ==============================================================================
echo "[$(date)] Loading environment modules..."
source ~/.bashrc

# Safely enable strict variable evaluation now that user profile layers are established
set -u  

module purge 2>/dev/null || true
module load STAR/2.7.9a
module load samtools/1.23.1

export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# ==============================================================================
#                      PATH DEFINITIONS AND SANITY CHECKS
# ==============================================================================
OUTPUT_DIR="/path/to/output/directory"
INDEX_DIR="${OUTPUT_DIR}/STAR_combined_index"
BAM_RAW_DIR="${OUTPUT_DIR}/BAM_master_alignments"
BAM_FINAL_DIR="${OUTPUT_DIR}/BAM_haplotype_specific"

# FASTQ Source Directories
FASTQ_DIR_NEW="/path/to/input/directory"
FASTQ_DIR_OLD="/path/to/input/directory"

# Verify index directory integrity before launching resource-heavy alignment
if [ ! -f "${INDEX_DIR}/Genome" ]; then
    echo "ERROR: Static genome index folder missing or incomplete inside ${INDEX_DIR}." >&2
    exit 1
fi

mkdir -p "$BAM_RAW_DIR" "$BAM_FINAL_DIR"

# ==============================================================================
#            STEP 1: MULTI-SAMPLE COMPETITIVE ALIGNMENT LOOPS
# ==============================================================================
THREADS=${SLURM_CPUS_PER_TASK}
echo "[$(date)] Starting competitive read mapping loop using ${THREADS} threads..."

# --- ARRAY A: Process New 2x150 nt Sequencing Datasets ---
SAMPLES_NEW=("Leaves" "Roots")

for SAMPLE in "${SAMPLES_NEW[@]}"; do
    R1="${FASTQ_DIR_NEW}/${SAMPLE}_merged_R1.fastq.gz"
    R2="${FASTQ_DIR_NEW}/${SAMPLE}_merged_R2.fastq.gz"
    PREFIX="${BAM_RAW_DIR}/${SAMPLE}_"
    MASTER_BAM="${PREFIX}Aligned.sortedByCoord.out.bam"
    
    if [ ! -f "$MASTER_BAM" ]; then
        echo "  -> Mapping 2x150nt dataset: ${SAMPLE}..."
        STAR --runThreadN "$THREADS" \
             --genomeDir "$INDEX_DIR" \
             --readFilesIn "$R1" "$R2" \
             --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMstrandField intronMotif \
             --limitBAMsortRAM 60000000000 \
             --outFileNamePrefix "$PREFIX"
        
        samtools index "$MASTER_BAM"
    fi
done

# --- ARRAY B: Process Old 2x100 nt Sequencing Datasets ---
SAMPLES_OLD=("Stem" "Bud" "InPetal" "Leaf" "OutPetal" "Pistil")

for SAMPLE in "${SAMPLES_OLD[@]}"; do
    R1="${FASTQ_DIR_OLD}/${SAMPLE}_R1_merged.fastq.gz"
    R2="${FASTQ_DIR_OLD}/${SAMPLE}_R2_merged.fastq.gz"
    PREFIX="${BAM_RAW_DIR}/${SAMPLE}_2012_"
    MASTER_BAM="${PREFIX}Aligned.sortedByCoord.out.bam"
    
    if [ ! -f "$MASTER_BAM" ]; then
        echo "  -> Mapping 2x100nt dataset: ${SAMPLE}..."
        STAR --runThreadN "$THREADS" \
             --genomeDir "$INDEX_DIR" \
             --readFilesIn "$R1" "$R2" \
             --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMstrandField intronMotif \
             --limitBAMsortRAM 60000000000 \
             --outFileNamePrefix "$PREFIX"
        
        samtools index "$MASTER_BAM"
    fi
done

# ==============================================================================
#         STEP 2: COORDINATE-BASED DE-MULTIPLEXING (BAM DE-CONVOLUTION)
# ==============================================================================
echo "[$(date)] Commencing coordinate de-multiplexing via samtools headers matching..."

HAP1_BAM_FILES=()
HAP2_BAM_FILES=()
ALL_PROCESSED_BAMS=($(ls ${BAM_RAW_DIR}/*Aligned.sortedByCoord.out.bam))

for RAW_BAM in "${ALL_PROCESSED_BAMS[@]}"; do
    BASE_NAME=$(basename "$RAW_BAM" "Aligned.sortedByCoord.out.bam")
    
    HAP1_SPLIT_BAM="${BAM_FINAL_DIR}/${BASE_NAME}Hap1.bam"
    HAP2_SPLIT_BAM="${BAM_FINAL_DIR}/${BASE_NAME}Hap2.bam"
    
    echo "  -> Splitting target mapping array: ${BASE_NAME}"
    
    HAP1_CHRS=$(samtools idxstats "$RAW_BAM" | cut -f1 | grep '_Hap1$' || true)
    HAP2_CHRS=$(samtools idxstats "$RAW_BAM" | cut -f1 | grep '_Hap2$' || true)
    
    if [ ! -f "$HAP1_SPLIT_BAM" ] && [ -n "$HAP1_CHRS" ]; then
        samtools view -b -@ "$THREADS" "$RAW_BAM" ${HAP1_CHRS} > "$HAP1_SPLIT_BAM"
        samtools index "$HAP1_SPLIT_BAM"
    fi
    
    if [ ! -f "$HAP2_SPLIT_BAM" ] && [ -n "$HAP2_CHRS" ]; then
        samtools view -b -@ "$THREADS" "$RAW_BAM" ${HAP2_CHRS} > "$HAP2_SPLIT_BAM"
        samtools index "$HAP2_SPLIT_BAM"
    fi
    
    HAP1_BAM_FILES+=("$HAP1_SPLIT_BAM")
    HAP2_BAM_FILES+=("$HAP2_SPLIT_BAM")
done

# ==============================================================================
#                    DOWNSTREAM BRAKER3 ARGS RECORD SUMMARY
# ==============================================================================
echo "========================================================================"
echo " PHASE 2 SUCCESS: UNCONTAMINATED COMPETITIVE TRACKS GENERATED"
echo "========================================================================"
echo "Hap1 Comma-Separated Input String for Phase 3 Script:"
(IFS=,; echo "${HAP1_BAM_FILES[*]}")
echo ""
echo "Hap2 Comma-Separated Input String for Phase 3 Script:"
(IFS=,; echo "${HAP2_BAM_FILES[*]}")
echo "========================================================================"
echo "[$(date)] Pipeline workflow finished execution successfully."