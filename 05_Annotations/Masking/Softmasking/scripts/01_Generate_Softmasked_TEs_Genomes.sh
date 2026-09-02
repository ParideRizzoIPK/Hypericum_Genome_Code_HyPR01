#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G             # Moderate memory footprint since libraries are pre-compiled
#SBATCH --time=04:00:00       # 4-hour window is highly conservative for a targeted library scan
#SBATCH --job-name=R2_softmask
#SBATCH --output=logs/softmask_%j.out
#SBATCH --error=logs/softmask_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  

echo "#########################################################"
echo "Pipeline: Step 1 (Script R2) - Targeted Soft-Masking Generation"
echo "Date: $(date)"
echo "#########################################################"

# --- 1. ENVIRONMENT & MODULE INITIALIZATION ---
source ~/.bashrc

# Load the core repeat masking software suite
module load EDTA/2.2.2

set -u            

export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# --- 2. CONFIGURATION AND DIRECTORIES ---
MASTER_OUT_DIR="/path/to/output/directory"
mkdir -p "${MASTER_OUT_DIR}"

# Source locations from Step R1 Workspaces
R1_DIR="/path/to/your/directory"
HAP1_R1_Wspace="${R1_DIR}/Hap1_EDTA_wspace"
HAP2_R1_Wspace="${R1_DIR}/Hap2_EDTA_wspace"

# Define exact original base file handles used in Step R1
HAP1_BASE="Hap1_curated_chr_un_gapclosed.fasta"
HAP2_BASE="Hap2_curated_chr_un_gapclosed.fasta"

THREADS=${SLURM_CPUS_PER_TASK}

# ==============================================================================
# MODULE R2.1: TARGETED SOFT-MASKING FOR HAPLOTYPE 1
# ==============================================================================
echo ">>>> PROCESSING HAPLOTYPE 1 SOFT-MASKING <<<<"
HAP1_LOCAL_WSPACE="${MASTER_OUT_DIR}/Hap1_Softmask_wspace"
mkdir -p "${HAP1_LOCAL_WSPACE}"
cd "${HAP1_LOCAL_WSPACE}"

# Symlink the unmasked, header-normalized genome sequence and custom library from R1
ln -sf "${HAP1_R1_Wspace}/${HAP1_BASE}.mod" ./"${HAP1_BASE}.mod"
ln -sf "${HAP1_R1_Wspace}/${HAP1_BASE}.mod.EDTA.TElib.fa" ./"${HAP1_BASE}.mod.EDTA.TElib.fa"

echo "  -> Running RepeatMasker with -xsmall against Haplotype 1 custom TE library..."
RepeatMasker \
    -pa "${THREADS}" \
    -xsmall \
    -dir . \
    -lib "${HAP1_BASE}.mod.EDTA.TElib.fa" \
    "${HAP1_BASE}.mod"

# Standardize output naming convention for production tracking
mv "${HAP1_BASE}.mod.masked" "${MASTER_OUT_DIR}/hap1.softmasked.fa"


# ==============================================================================
# MODULE R2.2: TARGETED SOFT-MASKING FOR HAPLOTYPE 2
# ==============================================================================
echo ">>>> PROCESSING HAPLOTYPE 2 SOFT-MASKING <<<<"
HAP2_LOCAL_WSPACE="${MASTER_OUT_DIR}/Hap2_Softmask_wspace"
mkdir -p "${HAP2_LOCAL_WSPACE}"
cd "${HAP2_LOCAL_WSPACE}"

# Symlink the unmasked, header-normalized genome sequence and custom library from R1
ln -sf "${HAP2_R1_Wspace}/${HAP2_BASE}.mod" ./"${HAP2_BASE}.mod"
ln -sf "${HAP2_R1_Wspace}/${HAP2_BASE}.mod.EDTA.TElib.fa" ./"${HAP2_BASE}.mod.EDTA.TElib.fa"

echo "  -> Running RepeatMasker with -xsmall against Haplotype 2 custom TE library..."
RepeatMasker \
    -pa "${THREADS}" \
    -xsmall \
    -dir . \
    -lib "${HAP2_BASE}.mod.EDTA.TElib.fa" \
    "${HAP2_BASE}.mod"

# Standardize output naming convention for production tracking
mv "${HAP2_BASE}.mod.masked" "${MASTER_OUT_DIR}/hap2.softmasked.fa"


# ==============================================================================
# MODULE R2.3: COMPLIANCE & INTEGRITY AUDIT
# ==============================================================================
echo ">>>> RUNNING PRODUCTION COMPLIANCE AUDITS <<<<"
cd "${MASTER_OUT_DIR}"

echo "  -> Verifying soft-mask case conversion ratios..."
HAP1_COUNT=$(grep -v "^>" hap1.softmasked.fa | grep -o "[a-z]" | head -n 100 | wc -l || true)
HAP2_COUNT=$(grep -v "^>" hap2.softmasked.fa | grep -o "[a-z]" | head -n 100 | wc -l || true)

if [ "${HAP1_COUNT}" -gt 0 ] && [ "${HAP2_COUNT}" -gt 0 ]; then
    echo "SUCCESS: Valid soft-masked lowercase sequence structures verified."
    # Clean large intermediate track directories if validation clears
    rm -rf "${HAP1_LOCAL_WSPACE}" "${HAP2_LOCAL_WSPACE}"
else
    echo "ERROR: Softmasked verification failed. Lowercase nucleotide array is empty."
    exit 1
fi

echo "#########################################################"
echo "SCRIPT R2 COMPLETE. PRODUCTION SEQUENCES GENERATED:"
echo "  - Hap1 Production: ${MASTER_OUT_DIR}/hap1.softmasked.fa"
echo "  - Hap2 Production: ${MASTER_OUT_DIR}/hap2.softmasked.fa"
echo "#########################################################"