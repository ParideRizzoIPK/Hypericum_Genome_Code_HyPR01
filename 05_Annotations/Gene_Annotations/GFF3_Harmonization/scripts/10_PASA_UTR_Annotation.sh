#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32                    # Maximize multithreading for GMAP alignment loops
#SBATCH --mem=150G                            # High-capacity buffer for loading deep transcript assemblies
#SBATCH --time=3-00:00:00                     # 3 days walltime limit allocation
#SBATCH --job-name=pasa_utr
#SBATCH --output=logs/pasa_utr_%j.out
#SBATCH --error=logs/pasa_utr_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  # Defer strict -u variable enforcement until profile environment variables load cleanly

# ==============================================================================
#                       ENVIRONMENT SETUP & MODULES
# ==============================================================================
echo "[$(date)] Loading user profile layers and cluster modules..."
source ~/.bashrc

set -u  # Activate strict variable enforcement safely after environment profile is loaded

module purge 2>/dev/null || true
module load gmap/2019-09-12                    # Load verified standalone cluster mapping utility
module load samtools/1.23.1

# Grant Apptainer runtime visibility to the cluster storage arrays
export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# Define immutable core pathing arrays
CONTAINER_IMG="/opt/Bio/PASApipeline/2.5.2/bin/pasapipeline.sif"
BASE_OUT="/path/to/output/directory"

HAPLOTYPES=("Hap1" "Hap2")

# Verify container presence before initiating job track
if [ ! -f "$CONTAINER_IMG" ]; then
    echo "ERROR: Target Apptainer image footprint missing at: $CONTAINER_IMG" >&2
    exit 1
fi

# ==============================================================================
#                         PIPELINE LOOP TRACK
# ==============================================================================
for HAP in "${HAPLOTYPES[@]}"; do
    echo "======================================================================"
    echo "[$(date)] INITIALIZING PASA UTR ANNOTATION ENGINE FOR: ${HAP}"
    echo "======================================================================"
    
    # Establish local workspace
    WDIR="${BASE_OUT}/${HAP}"
    mkdir -p "$WDIR"
    cd "$WDIR"

    # Define haplotype-specific asset paths
    if [ "$HAP" == "Hap1" ]; then
        GENOME="/path/to/input/directory/hap1.masked_unchr.fa"
        TRANSCRIPTS="/path/to/input/directory/transcripts_merged.fasta"
        HARMONIZED_GFF="/path/to/input/directory/harmonized_consensus_Hap1.gff3"
    else
        GENOME="/path/to/input/directory/hap2.masked_unchr.fa"
        TRANSCRIPTS="/path/to/input/directory/transcripts_merged.fasta"
        HARMONIZED_GFF="/path/to/input/directory/harmonized_consensus_Hap2.gff3"
    fi

    # Verify input integrity
    for asset in "$GENOME" "$TRANSCRIPTS" "$HARMONIZED_GFF"; do
        if [ ! -f "$asset" ]; then
            echo "ERROR: Required input asset missing for ${HAP}: $asset" >&2
            exit 1
        fi
    done

    # --------------------------------------------------------------------------
    # STEP 1: FIX - THOROUGH PRE-FLIGHT PURGING (CLEARING HIDDEN CHECKPOINTS)
    # --------------------------------------------------------------------------
    echo "[$(date)] Performing automated workspace cleansing of old failed trackers..."
    # Explicitly including the double-underscore checkpoint folder to prevent false resume skips
    rm -rf pasa_db.sqlite* __pasa_pasa_db.sqlite_SQLite_chkpts *.status pasa_project.* tmp.* GMAP_dir *.gene_structures_post_*_updates.*.gff3

    # --------------------------------------------------------------------------
    # STEP 2: DYNAMIC ABSOLUTE PASA CONFIGURATION INJECTION
    # --------------------------------------------------------------------------
    echo "[$(date)] Injecting custom localized SQLite absolute path parameters..."
    cat << EOF > pasa.config
# PASA Pipeline Configuration for Local SQLite Database
DATABASE=${WDIR}/pasa_db.sqlite
EOF

    # --------------------------------------------------------------------------
    # STEP 3: PHASE A - TRANSCRIPT ALIGNMENT ASSEMBLY
    # --------------------------------------------------------------------------
    echo "[$(date)] Launching Phase A: Multi-threaded GMAP alignment assembly..."
    
    apptainer exec "$CONTAINER_IMG" /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
        -c pasa.config \
        -C -R \
        -g "$GENOME" \
        -t "$TRANSCRIPTS" \
        --ALIGNERS gmap \
        --CPU "${SLURM_CPUS_PER_TASK}"

    echo "[$(date)] Phase A completed successfully. Transcript segments grouped into splice graphs."

    # --------------------------------------------------------------------------
    # STEP 4: PHASE B - REFERENCE COMPARISON & UTR EXPANSION
    # --------------------------------------------------------------------------
    echo "[$(date)] Launching Phase B: Injecting UTR vectors into harmonized coding sequences..."
    
    apptainer exec "$CONTAINER_IMG" /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
        -c pasa.config \
        -A \
        -g "$GENOME" \
        -t "$TRANSCRIPTS" \
        -L --annots "$HARMONIZED_GFF"

    # --------------------------------------------------------------------------
    # STEP 5: COMPILING FINAL STANDARDIZED OUTPUT
    # --------------------------------------------------------------------------
    echo "[$(date)] Extracting and formatting final GFF3 curation tracks..."
    
    FINAL_PASA_GFF=$(ls -t ${WDIR}/pasa_db.sqlite.gene_structures_post_[pP][aS][sS][aA]_updates.*.gff3 2>/dev/null | head -n 1 || true)
    
    if [ -n "$FINAL_PASA_GFF" ] && [ -f "$FINAL_PASA_GFF" ]; then
        cp "$FINAL_PASA_GFF" "${WDIR}/final_harmonized_with_utr_${HAP}.gff3"
        echo "✔ SUCCESS: Finalized UTR-enriched track compiled at: ${WDIR}/final_harmonized_with_utr_${HAP}.gff3"
    else
        echo "WARNING: Could not automatically map internal database output file names. Manual file verification required in: $WDIR"
    fi

done

echo "======================================================================"
echo "[$(date)] PASA COMPREHENSIVE PIPELINE EXECUTION FINISHED CLEANLY"
echo "======================================================================"