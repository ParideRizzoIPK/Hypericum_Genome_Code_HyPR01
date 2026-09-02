#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G            
#SBATCH --time=48:00:00       # Extended to safely cover both EDTA execution and the secondary LAI mappings
#SBATCH --job-name=EDTA_LAI
#SBATCH --output=logs/edta_%j.out
#SBATCH --error=logs/edta_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  

echo "#########################################################"
echo "Pipeline: Step 1 (Script R1) - EDTA TE Discovery & LAI"
echo "Date: $(date)"
echo "#########################################################"

# --- 1. ENVIRONMENT & MODULE INITIALIZATION ---
source ~/.bashrc

module load EDTA/2.2.2

# Safely engage strict validation for unset variables post-load
set -u            

export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# --- 2. PATHS & RESOURCE ALLOCATION ---
MASTER_OUT_DIR="/path/to/output/directory"
mkdir -p "${MASTER_OUT_DIR}"

HAP1_IN_FA="/path/to/input/directory/Hap1_curated_chr_un_gapclosed.fasta"
HAP2_IN_FA="/path/to/input/directory/Hap2_curated_chr_un_gapclosed.fasta"

TOTAL_THREADS=${SLURM_CPUS_PER_TASK}
HALF_THREADS=$(( TOTAL_THREADS / 2 ))

echo "Allocated total CPUs: ${TOTAL_THREADS} cores."
echo "Splitting resource footprint: ${HALF_THREADS} threads per haplotype run."

# ==============================================================================
# PIPELINE EXECUTION LOOP: DUAL BACKGROUND RUNS
# ==============================================================================

# --- RUNTRACK 1: HAPLOTYPE 1 ---
echo ">>>> INITIALIZING HAPLOTYPE 1 EDTA & LAI BACKGROUND INSTANCE <<<<"
HAP1_WSPACE="${MASTER_OUT_DIR}/Hap1_EDTA_wspace"
mkdir -p "${HAP1_WSPACE}"

(
    cd "${HAP1_WSPACE}"
    # Dynamically capture and preserve the exact assembly filename
    HAP1_BASE=$(basename "${HAP1_IN_FA}")
    ln -sf "${HAP1_IN_FA}" ./"${HAP1_BASE}"
    
    echo "[$(date)] Starting EDTA for ${HAP1_BASE}..." > edta_run_hap1.log
    EDTA.pl \
        --genome "${HAP1_BASE}" \
        --species others \
        --step all \
        --anno 1 \
        --overwrite 1 \
        --threads "${HALF_THREADS}" >> edta_run_hap1.log 2>&1
        
    echo "[$(date)] EDTA complete. Mapping LTRs for LAI calculation..." >> edta_run_hap1.log
    RepeatMasker -pa "${HALF_THREADS}" -q -no_is -norna -nolow -div 40 \
        -lib "${HAP1_BASE}.mod.EDTA.raw/LTR/${HAP1_BASE}.mod.LTRlib.fa" \
        "${HAP1_BASE}" >> edta_run_hap1.log 2>&1
        
    echo "[$(date)] Calculating LAI..." >> edta_run_hap1.log
    /opt/Bio/EDTA/2.2.2/share/LTR_retriever/LAI \
        -genome "${HAP1_BASE}" \
        -intact "${HAP1_BASE}.mod.EDTA.raw/LTR/${HAP1_BASE}.mod.pass.list" \
        -all "${HAP1_BASE}.out" >> edta_run_hap1.log 2>&1
        
    echo "[$(date)] Haplotype 1 block completely finished." >> edta_run_hap1.log
) &
HAP1_PID=$!
echo "Haplotype 1 process successfully detached under PID: ${HAP1_PID}"


# --- RUNTRACK 2: HAPLOTYPE 2 ---
echo ">>>> INITIALIZING HAPLOTYPE 2 EDTA & LAI BACKGROUND INSTANCE <<<<"
HAP2_WSPACE="${MASTER_OUT_DIR}/Hap2_EDTA_wspace"
mkdir -p "${HAP2_WSPACE}"

(
    cd "${HAP2_WSPACE}"
    HAP2_BASE=$(basename "${HAP2_IN_FA}")
    ln -sf "${HAP2_IN_FA}" ./"${HAP2_BASE}"
    
    echo "[$(date)] Starting EDTA for ${HAP2_BASE}..." > edta_run_hap2.log
    EDTA.pl \
        --genome "${HAP2_BASE}" \
        --species others \
        --step all \
        --anno 1 \
        --overwrite 1 \
        --threads "${HALF_THREADS}" >> edta_run_hap2.log 2>&1
        
    echo "[$(date)] EDTA complete. Mapping LTRs for LAI calculation..." >> edta_run_hap2.log
    RepeatMasker -pa "${HALF_THREADS}" -q -no_is -norna -nolow -div 40 \
        -lib "${HAP2_BASE}.mod.EDTA.raw/LTR/${HAP2_BASE}.mod.LTRlib.fa" \
        "${HAP2_BASE}" >> edta_run_hap2.log 2>&1
        
    echo "[$(date)] Calculating LAI..." >> edta_run_hap2.log
    /opt/Bio/EDTA/2.2.2/share/LTR_retriever/LAI \
        -genome "${HAP2_BASE}" \
        -intact "${HAP2_BASE}.mod.EDTA.raw/LTR/${HAP2_BASE}.mod.pass.list" \
        -all "${HAP2_BASE}.out" >> edta_run_hap2.log 2>&1
        
    echo "[$(date)] Haplotype 2 block completely finished." >> edta_run_hap2.log
) &
HAP2_PID=$!
echo "Haplotype 2 process successfully detached under PID: ${HAP2_PID}"


# ==============================================================================
# WAITING BARRIER & INTEGRITY CHECKS
# ==============================================================================
echo ">>>> Synchronization barrier engaged: Awaiting both background runs... <<<<"

wait "${HAP1_PID}"
echo "  -> Process Haplotype 1 (${HAP1_PID}) terminated."

wait "${HAP2_PID}"
echo "  -> Process Haplotype 2 (${HAP2_PID}) terminated."

echo ">>>> Running validation audits on generated outputs... <<<<"

HAP1_BASE=$(basename "${HAP1_IN_FA}")
if [ ! -f "${HAP1_WSPACE}/${HAP1_BASE}.mod.EDTA.TElib.fa" ]; then
    echo "CRITICAL WARNING: Haplotype 1 failed to produce the required TE Library. Check log at: ${HAP1_WSPACE}/edta_run_hap1.log"
else
    echo "  -> Haplotype 1 completed seamlessly."
fi

HAP2_BASE=$(basename "${HAP2_IN_FA}")
if [ ! -f "${HAP2_WSPACE}/${HAP2_BASE}.mod.EDTA.TElib.fa" ]; then
    echo "CRITICAL WARNING: Haplotype 2 failed to produce the required TE Library. Check log at: ${HAP2_WSPACE}/edta_run_hap2.log"
else
    echo "  -> Haplotype 2 completed seamlessly."
fi

echo "#########################################################"
echo "EDTA PIPELINE COMPLETE."
echo "Isolated output directories ready for Script R2 (Softmasking):"
echo "  - Hap1 Workspace: ${HAP1_WSPACE}"
echo "  - Hap2 Workspace: ${HAP2_WSPACE}"
echo "#########################################################"