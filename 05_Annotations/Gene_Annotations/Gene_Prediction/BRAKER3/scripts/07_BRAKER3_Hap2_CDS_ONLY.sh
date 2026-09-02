#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32                    # Maximize multithreading across alignment bins
#SBATCH --mem=150G                            # Safe buffer allocation for deep hint maps
#SBATCH --time=2-00:00:00                     # Streamlined 2-day ceiling cushion (will finish quickly)
#SBATCH --job-name=braker3_cds_h2
#SBATCH --output=logs/braker3_cds_h2_%j.out
#SBATCH --error=logs/braker3_cds_h2_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  # Defer -u flag until cluster profile layers load cleanly

# ==============================================================================
#                       ENVIRONMENT MODULE LOADING
# ==============================================================================
echo "[$(date)] Loading environment profiles..."
source ~/.bashrc

set -u  # Activate strict variable enforcement after profile validation

module purge 2>/dev/null || true
module load STAR/2.7.9a
module load BRAKER/3.0.8
module load samtools/1.23.1

# Grant the underlying core tools full visibility to the shared storage arrays
export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# ==============================================================================
#                      PATH DEFINITIONS AND SANITY CHECKS
# ==============================================================================
BASE_OUT="/path/to/output/directory"
BRAKER_OUT_DIR="${BASE_OUT}/Hap2/BRAKER3_results"
AUG_DIR="${BASE_OUT}/Hap2/augustus_config_dir"

mkdir -p "$BRAKER_OUT_DIR" "$AUG_DIR"

# Core Structural Assets (Targeting your softmasked Haplotype 2 assemblies)
GENOME="/path/to/input/directory/hap2.softmasked.fa"
PROT_DB="${BASE_OUT}/Hypericum_Comprehensive_Prot_Evidence.fasta"
CLEANED_BAM_DIR="/path/to/your/directory"

SPECIES_NAME="hypericum_perforatum_hap2_v2"
THREADS="${SLURM_CPUS_PER_TASK}"

# Compile comma-separated alignment string using the validated, purified Haplotype 2 header tracks
BAM_LIST="${CLEANED_BAM_DIR}/Bud_2012_Hap2.bam,${CLEANED_BAM_DIR}/InPetal_2012_Hap2.bam,${CLEANED_BAM_DIR}/Leaf_2012_Hap2.bam,${CLEANED_BAM_DIR}/Leaves_Hap2.bam,${CLEANED_BAM_DIR}/OutPetal_2012_Hap2.bam,${CLEANED_BAM_DIR}/Pistil_2012_Hap2.bam,${CLEANED_BAM_DIR}/Roots_Hap2.bam,${CLEANED_BAM_DIR}/Stem_2012_Hap2.bam"

# Verify all required physical assets exist on `/path/to/your/directory` before executing pipeline
for asset in "$GENOME" "$PROT_DB"; do
    if [ ! -f "$asset" ]; then
        echo "ERROR: Missing primary verification asset: $asset" >&2
        exit 1
    fi
done

# ==============================================================================
#             PRE-FLIGHT RUNTIME CONFIGURATION STAGING
# ==============================================================================
# Setup writable Augustus configuration layout to prevent system write locks inside container
if [ ! -d "${AUG_DIR}/config" ]; then
    echo "Staging local Augustus config structures..."
    cp -r "$(dirname "$(which augustus)")/../config" "$AUG_DIR/"
fi
export AUGUSTUS_CONFIG_PATH="${AUG_DIR}/config"

# Validate local GeneMark execution license configuration
if [ ! -f ~/.gm_key ]; then
    echo "Staging GeneMark execution authorization tracking licenses..."
    zcat /path/to/your/directory/GeneMark_ETP_Software_Data/gm_key.gz > ~/.gm_key
fi

# Re-map tool paths to bypass cluster configuration layer abstractions
unset GENEMARK_PATH
export PATH="/path/to/your/directory:$PATH"
GM_BIN_DIR="/path/to/your/directory"

# ==============================================================================
#                          BRAKER3 COMPILATION ENGINE
# ==============================================================================
echo "[$(date)] Starting clean Haplotype 2 BRAKER3 CDS-only annotation run..."

braker.pl \
    --genome="$GENOME" \
    --bam="$BAM_LIST" \
    --prot_seq="$PROT_DB" \
    --species="$SPECIES_NAME" \
    --softmasking \
    --gff3 \
    --workingdir="$BRAKER_OUT_DIR" \
    --threads="$THREADS" \
    --AUGUSTUS_CONFIG_PATH="$AUGUSTUS_CONFIG_PATH" \
    --GENEMARK_PATH="$GM_BIN_DIR"

echo "[$(date)] Pipeline successfully completed execution."