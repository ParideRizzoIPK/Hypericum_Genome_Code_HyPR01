#!/bin/bash
#SBATCH --job-name=HyPR01_DualHap_InterProScan
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=logs/InterProScan_DualHap_%j.out
#SBATCH --error=logs/InterProScan_DualHap_%j.err

# ────────────────────────────────────────────────────────
# SCRIPT METADATA & VERSION CONTROL
# Pipeline Version: 1.3.2 (Pure Search Mode)
# Last Modified: 2026-07-13
# Description: Dual-haplotype InterProScan compute engine.
#              Peptide extraction via gffread with standalone
#              HMM domain database execution profiles.
# ────────────────────────────────────────────────────────

set -e

# ────────────────────────────────────────────────────────
# 1. MODULE LOADING & ENVIRONMENT
# ────────────────────────────────────────────────────────
. /etc/profile.d/modules.sh
module load gffread/0.12.6
module load interproscan/5.57-90.0

# --- ASSEMBLY & DATA PATHS ---
GENOME_HAP1="/path/to/input/directory/hap1.masked_unchr.fa"
GENOME_HAP2="/path/to/input/directory/hap2.masked_unchr.fa"

# GFF3 Coordinate Source Directory (EggNOG directory tree)
GFF3_SRC_DIR="/path/to/input/directory"

# Dedicated InterProScan Directory Tree
BASE_INTERPRO_DIR="/path/to/your/directory"

# ────────────────────────────────────────────────────────
# 2. RUN LOOP PROCESSING STRATEGY
# ────────────────────────────────────────────────────────
declare -A HAP_CONFIGS
HAP_CONFIGS[Hap1]="$GENOME_HAP1"
HAP_CONFIGS[Hap2]="$GENOME_HAP2"

for HAP in Hap1 Hap2; do
    GENOME_FASTA=${HAP_CONFIGS[$HAP]}
    
    # Define and create unique InterProScan subdirectories
    TARGET_DIR="${BASE_INTERPRO_DIR}/HyPR01_${HAP}"
    mkdir -p "$TARGET_DIR"
    
    # Source GFF3 from the EGGNOG directory path where they live
    INPUT_GFF="${GFF3_SRC_DIR}/HyPR01_${HAP}/HyPR01_${HAP}_final_annotations.gff3"
    
    EXTRACTED_FAA="${TARGET_DIR}/HyPR01_${HAP}_interpro_proteins_raw.faa"
    CLEANED_FAA="${TARGET_DIR}/HyPR01_${HAP}_interpro_proteins_cleaned.faa"
    OUT_PREFIX="${TARGET_DIR}/HyPR01_${HAP}_interproscan"
    
    TEMP_FASTA="${TARGET_DIR}/temp_genome_sorted_headers.fa"

    echo "========================================================================"
    echo "[$(date)] STARTING INTERPROSCAN WORKFLOW FOR TRACK: $HAP"
    echo "========================================================================"
    echo "InterProScan Target Dir: $TARGET_DIR"
    echo "Genome Reference:        $GENOME_FASTA"
    echo "Coordinate Source GFF3:  $INPUT_GFF"
    
    if [ ! -f "$INPUT_GFF" ]; then
        echo "[ERROR] Missing input GFF3 coordinates: $INPUT_GFF. Skipping."
        continue
    fi

    # Phase A: Re-extracting and sanitizing proteins to eliminate bad translation tokens
    echo "[-] Step 1: Normalizing fasta headers and running gffread translation..."
    sed 's/^>UnChr/>ChrUn/' "$GENOME_FASTA" > "$TEMP_FASTA"
    gffread "$INPUT_GFF" -g "$TEMP_FASTA" -V -y "$EXTRACTED_FAA"
    rm -f "$TEMP_FASTA"

    # Secondary deep clean pass: Strip any trailing or internal '.' or '*' characters
    sed '/^[^>]/s/[\.\*]//g' "$EXTRACTED_FAA" > "$CLEANED_FAA"
    rm -f "$EXTRACTED_FAA"

    # Quantify sequence totals
    PROTEIN_COUNT=$(grep -c "^>" "$CLEANED_FAA" || true)
    echo "[+] Extraction complete. Total proteins available for InterProScan: $PROTEIN_COUNT"

    # Phase B: Launching InterProScan Core Engine (ProDom removed)
    echo "[-] Step 2: Launching InterProScan Core Engine..."
    interproscan.sh \
        -i "$CLEANED_FAA" \
        -b "$OUT_PREFIX" \
        -f TSV,GFF3,XML \
        --goterms \
        --pathways \
        -appl Pfam,PANTHER,SMART,PRINTS,SUPERFAMILY,TIGRFAM,Gene3D,CDD \
        --cpu "$SLURM_CPUS_PER_TASK"

    # Clean up local peptide fasta to save space
    rm -f "$CLEANED_FAA"
    echo "[+] InterProScan run successfully finished for $HAP."
    echo ""
done

echo "[*] Dual Haplotype InterProScan Run completely finished at $(date)"