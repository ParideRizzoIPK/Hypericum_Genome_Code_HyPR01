#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --job-name=HyPR01_BAM_audit
#SBATCH --partition=cpu                       # CPU compute queue
#SBATCH --cpus-per-task=8                     # Multithreaded calculations
#SBATCH --mem=40G                             # Safe memory headroom for array lookups
#SBATCH --time=02:00:00                       # Queue-optimized walltime for instant backfill slots
#SBATCH --output=logs/bam_audit_%j.out
#SBATCH --error=logs/bam_audit_%j.err

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

# Grant containerized samtools visibility to your storage arrays
export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# ==============================================================================
#                      PATH DEFINITIONS AND SANITY CHECKS
# ==============================================================================
BAM_DIR="/path/to/your/directory"
FASTA_HAP1="/path/to/input/directory/hap1.masked_unchr.fa"
FASTA_HAP2="/path/to/input/directory/hap2.masked_unchr.fa"

REPORT_FILE="${BAM_DIR}/BAM_evidence_audit_report.txt"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Initialize a clean audit report file directly inside the target BAM directory
cat << EOF > "$REPORT_FILE"
================================================================================
            HYPERICUM GENOME RNA-SEQ EVIDENCE QC & INTEGRITY REPORT
================================================================================
Generated on : $(date)
Job ID       : ${SLURM_JOB_ID:-N/A}
Target Dir   : $BAM_DIR
================================================================================

EOF

# ==============================================================================
#          STEP 1: GENERATE TRUTHS FROM HARD-MASKED REFERENCE FASTAS
# ==============================================================================
echo "[$(date)] Measuring true reference assembly sizes..."

# Ensure FASTA files have active tracking indexes (.fai)
[ -f "${FASTA_HAP1}.fai" ] || samtools faidx "$FASTA_HAP1"
[ -f "${FASTA_HAP2}.fai" ] || samtools faidx "$FASTA_HAP2"

# Parse total nucleotide sizes from index arrays
HAP1_FASTA_SIZE=$(awk '{sum+=$2} END {print sum}' "${FASTA_HAP1}.fai")
HAP2_FASTA_SIZE=$(awk '{sum+=$2} END {print sum}' "${FASTA_HAP2}.fai")

cat << EOF >> "$REPORT_FILE"
### 1. REFERENCE GENOME BASELINES (TRUE ASSEMBLY SIZE)
* Haplotype 1 True Nucleotide Count : $HAP1_FASTA_SIZE bp
* Haplotype 2 True Nucleotide Count : $HAP2_FASTA_SIZE bp

================================================================================
### 2. INDIVIDUAL PURGED FILE AUDITS & BREADTH OF COVERAGE METRICS
================================================================================
EOF

# ==============================================================================
#          STEP 2: LOOP AND AUDIT PURGED FILE METRICS
# ==============================================================================
echo "[$(date)] Starting file header dictionary and coverage metrics scans..."

cd "$BAM_DIR"

# Loop over all purged BAM files in alphabetical order
for bam in *.bam; do
    [ -f "$bam" ] || continue
    echo "Auditing: $bam"
    
    # Identify which haplotype reference target this file maps to
    if [[ "$bam" == *Hap1* ]]; then
        TRUE_FASTA_SIZE="$HAP1_FASTA_SIZE"
        HAP_TAG="Hap1"
        OPPOSITE_TAG="_Hap2"
    else
        TRUE_FASTA_SIZE="$HAP2_FASTA_SIZE"
        HAP_TAG="Hap2"
        OPPOSITE_TAG="_Hap1"
    fi
    
    # A. PURGED HEADER DICTIONARY AUDIT
    # Check for any lingering references from BOTH the cross-haplotype tag or the alternative haplotype space
    LINGERING_TAGS=$(samtools view -H "$bam" | grep "@SQ" | grep -E "(_Hap1|_Hap2)" || true)
    CROSS_HAP_ENTRY=$(samtools view -H "$bam" | grep "@SQ" | grep "$OPPOSITE_TAG" || true)
    
    if [ -n "$LINGERING_TAGS" ] || [ -n "$CROSS_HAP_ENTRY" ]; then
        HEADER_STATUS="❌ FAILED (Extraneous sequences or unpurged cross-haplotype tags found!)"
    else
        HEADER_STATUS="✔ PASSED (BAM dictionary completely purified and confined to target assembly space)"
    fi
    
    # Capture unique target chromosome names listed in BAM to output to report
    CHR_LIST=$(samtools view -H "$bam" | grep "@SQ" | cut -f2 | sed 's/SN://g' | tr '\n' ' ' | sed 's/ $//')
    
    # B. GENOMIC BREADTH OF COVERAGE CALCULATION
    # Stream depth summary tables, filtering out header formatting
    COVERAGE_DATA=$(samtools coverage "$bam" | grep -v "^#" || true)
    
    # Sum total reference lengths declared in BAM, and total covered bases (at least 1 read)
    BAM_DECLARED_LEN=$(echo "$COVERAGE_DATA" | awk '{sum+=$3} END {print sum}')
    BAM_COVERED_BASES=$(echo "$COVERAGE_DATA" | awk '{sum+=$5} END {print sum}')
    
    # Set safe fallbacks if files are empty
    BAM_DECLARED_LEN=${BAM_DECLARED_LEN:-0}
    BAM_COVERED_BASES=${BAM_COVERED_BASES:-0}
    
    # Calculate exact genomic breadth coverage percentage against the TRUE assembly length
    PERCENT_COVERAGE=$(awk -v cov="$BAM_COVERED_BASES" -v tot="$TRUE_FASTA_SIZE" 'BEGIN {printf "%.4f", (cov/tot)*100}')
    
    # Write entries out to the centralized text report file
    cat << EOF >> "$REPORT_FILE"
File Name           : $bam
Target Haplotype    : $HAP_TAG
Header Audit Status : $HEADER_STATUS
Mapped Chromosomes  : $CHR_LIST
Metrics:
  -> Reference FASTA Size       : $TRUE_FASTA_SIZE bp
  -> Size Declared inside BAM   : $BAM_DECLARED_LEN bp
  -> Distinct Bases with Reads  : $BAM_COVERED_BASES bp
  -> True Genomic Breadth (%)   : $PERCENT_COVERAGE %
--------------------------------------------------------------------------------
EOF
done

echo "========================================================================" >> "$REPORT_FILE"
echo "                      END OF AUDIT REPORT                               " >> "$REPORT_FILE"
echo "========================================================================" >> "$REPORT_FILE"

echo "[$(date)] Audit pipeline complete. Summary report compiled inside your folder."