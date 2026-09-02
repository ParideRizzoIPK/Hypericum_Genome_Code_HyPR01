#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --job-name=HyPR01_BAM_purge
#SBATCH --partition=cpu                       # Standard CPU partition
#SBATCH --cpus-per-task=8                     # Multithreaded BAM encoding compression
#SBATCH --mem=40G                             # Safe memory headroom for on-the-fly streaming
#SBATCH --time=08:00:00                       # Generous execution cushion
#SBATCH --output=logs/bam_purge_%j.out
#SBATCH --error=logs/bam_purge_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  # Defer -u flag until cluster profile layers load cleanly

# ==============================================================================
#                       ENVIRONMENT MODULE LOADING
# ==============================================================================
echo "[$(date)] Loading environment profiles..."
source ~/.bashrc

set -u  # Activate strict variable enforcement after profile validation

module purge 2>/dev/null || true
module load samtools/1.23.1

# Grant the containerized samtools module full system-level visibility to shared NFS storage
export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# ==============================================================================
#                      PATH DEFINITIONS AND SANITY CHECKS
# ==============================================================================
BAM_RAW_DIR="/path/to/your/directory"
FIXED_OUT_DIR="/path/to/output/directory"

FASTA_HAP1="/path/to/input/directory/hap1.masked_unchr.fa"
FASTA_HAP2="/path/to/input/directory/hap2.masked_unchr.fa"

mkdir -p "$FIXED_OUT_DIR"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Verify reference genomes and index arrays exist before beginning stream conversions
echo "[$(date)] Indexing and measuring true reference assembly sizes..."
[ -f "${FASTA_HAP1}.fai" ] || samtools faidx "$FASTA_HAP1"
[ -f "${FASTA_HAP2}.fai" ] || samtools faidx "$FASTA_HAP2"

TMP_DIR="${FIXED_OUT_DIR}/tmp_hdr_purge_${SLURM_JOB_ID}"
mkdir -p "$TMP_DIR"

echo "========================================================================"
echo "          STARTING TRUE PURGE REMEDIATION REHEADER ENGINE              "
echo "========================================================================"
echo "Input Directory  : $BAM_RAW_DIR"
echo "Output Directory : $FIXED_OUT_DIR"
echo "========================================================================"

cd "$BAM_RAW_DIR"

# ==============================================================================
#               PURGE ENGINE: HAPLOTYPE 1 PROCESSING LAYER
# ==============================================================================
for bam in *Hap1.bam; do
    [ -f "$bam" ] || continue
    echo "[$(date)] Compiling clean dictionary and streaming alignment body for: $bam"
    
    HDR_FILE="${TMP_DIR}/${bam}_pure_h1.sam"
    
    # 1. Build a new header string using ONLY target haplotype coordinates
    echo "@HD	VN:1.6	SO:coordinate" > "$HDR_FILE"
    awk '{print "@SQ\tSN:"$1"\tLN:"$2}' "${FASTA_HAP1}.fai" >> "$HDR_FILE"
    
    # Append alternative header metadata blocks (Read groups, Program tags) while omitting old sequence records
    samtools view -H "$bam" | grep -E -v '^@HD|^@SQ' >> "$HDR_FILE"
    
    # 2. Stream, strip suffixes from text alignment fields, and encode a clean BAM binary
    (
        cat "$HDR_FILE"
        samtools view "$bam" | awk -v sep="\t" 'BEGIN {FS=OFS=sep} {gsub(/_Hap1/, "", $3); gsub(/_Hap1/, "", $7); print}'
    ) | samtools view -b -@ "$THREADS" - > "${FIXED_OUT_DIR}/$bam"
    
    # Generate clean coordinate index arrays
    samtools index -@ "$THREADS" "${FIXED_OUT_DIR}/$bam"
    echo "  -> Finished writing: ${FIXED_OUT_DIR}/$bam"
done

# ==============================================================================
#               PURGE ENGINE: HAPLOTYPE 2 PROCESSING LAYER
# ==============================================================================
for bam in *Hap2.bam; do
    [ -f "$bam" ] || continue
    echo "[$(date)] Compiling clean dictionary and streaming alignment body for: $bam"
    
    HDR_FILE="${TMP_DIR}/${bam}_pure_h2.sam"
    
    # 1. Build a new header string using ONLY target haplotype coordinates
    echo "@HD	VN:1.6	SO:coordinate" > "$HDR_FILE"
    awk '{print "@SQ\tSN:"$1"\tLN:"$2}' "${FASTA_HAP2}.fai" >> "$HDR_FILE"
    samtools view -H "$bam" | grep -E -v '^@HD|^@SQ' >> "$HDR_FILE"
    
    # 2. Stream, strip suffixes from text alignment fields, and encode a clean BAM binary
    (
        cat "$HDR_FILE"
        samtools view "$bam" | awk -v sep="\t" 'BEGIN {FS=OFS=sep} {gsub(/_Hap2/, "", $3); gsub(/_Hap2/, "", $7); print}'
    ) | samtools view -b -@ "$THREADS" - > "${FIXED_OUT_DIR}/$bam"
    
    samtools index -@ "$THREADS" "${FIXED_OUT_DIR}/$bam"
    echo "  -> Finished writing: ${FIXED_OUT_DIR}/$bam"
done

# ==============================================================================
#                             WORKSPACE CLEANUP
# ==============================================================================
rm -rf "$TMP_DIR"
echo "========================================================================"
echo " 🎉 SUCCESS: All target BAM coordinates have been purified and realigned."
echo " Total BAM files ready for BRAKER3: $(ls -1 ${FIXED_OUT_DIR}/*.bam | wc -l)"
echo " Finished at: $(date)"
echo "========================================================================"