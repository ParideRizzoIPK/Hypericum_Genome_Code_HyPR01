#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --job-name=HyPR01_ANNEVO_tracks
#SBATCH --partition=gpu                       
#SBATCH --account=your_slurm_account          # set to your cluster's GPU allocation/account name
#SBATCH --auks=yes                            # Enforce AUKS authentication ticket passing
#SBATCH --gres=gpu:1                          
#SBATCH --exclude=slurm-gpu-01,slurm-gpu-05   # Target cluster node exclusion filters
#SBATCH --cpus-per-task=12                    
#SBATCH --mem=64G                             
#SBATCH --time=48:00:00                       
#SBATCH --output=logs/annevo_tracks_%j.out
#SBATCH --error=logs/annevo_tracks_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  # Defer -u until after cluster environment module initializes paths

# ==============================================================================
#                       ENVIRONMENT MODULE LOADING
# ==============================================================================
echo "[$(date)] Loading cluster environment modules..."
source ~/.bashrc

module purge 2>/dev/null || true
module load ANNEVO/2.3.1                       # Load containerized environment module
module load gffread/0.12.6                    # Load protein extraction framework

set -u  # Activate strict unbound variable validation

# Grant the Apptainer container full visibility to your storage arrays
export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# True container path verified via your local filesystem investigation
TRUE_CONTAINER_PATH="/opt/Bio/ANNEVO/2.3.1/bin/container"
EXEC_CMD="apptainer run --nv $TRUE_CONTAINER_PATH"

# ==============================================================================
#                 TUNABLE WORKFLOW CONFIGURATION BLOCK
# ==============================================================================
BASE_OUT="/path/to/output/directory"

# Explicitly targeting your 15-class boundary-aware plant model file
MODEL_NAME="ANNEVO_Magnoliopsida.pt"
LINEAGE="Magnoliopsida" 

MODEL_PATH="${BASE_OUT}/${MODEL_NAME}"

echo "[$(date)] Running architecture layout: ${MODEL_NAME} with ${LINEAGE} parameter space"

# ==============================================================================
#                      PATH DEFINITIONS AND SANITY CHECKS
# ==============================================================================
# Haplotype-Specific Separate Subfolders
OUT_DIR_H1="${BASE_OUT}/Hap1"
OUT_DIR_H2="${BASE_OUT}/Hap2"

HAP1_GENOME="/path/to/input/directory/hap1.masked_unchr.fa"
HAP2_GENOME="/path/to/input/directory/hap2.masked_unchr.fa"

# Verify Essential Shared Files and uploaded weights exist before computing
for asset in "$MODEL_PATH" "$HAP1_GENOME" "$HAP2_GENOME"; do
    if [ ! -f "$asset" ]; then
        echo "ERROR: Required pipeline asset missing on cluster: $asset" >&2
        exit 1
    fi
done

# Initialize output subdirectories
mkdir -p "$OUT_DIR_H1" "$OUT_DIR_H2"
THREADS=${SLURM_CPUS_PER_TASK:-12}             

# Move directly into base workspace
cd "$BASE_OUT"

# ==============================================================================
#                 RUN LOOPS: SEQUENTIAL HAPLOTYPE PREDICTIONS
# ==============================================================================

# --- LOOP 1: Haplotype 1 (Hap1) Execution ---
echo "==================================================="
echo "[$(date)] STARTING STEP 1: ANNEVO Haplotype 1"
echo "==================================================="
H1_GFF="${OUT_DIR_H1}/annevo_raw_hap1.gff3"
H1_FAA="${OUT_DIR_H1}/annevo_raw_hap1_proteins.faa"

# --nv forwards host NVIDIA GPU drivers to compute models inside the container
$EXEC_CMD \
    --genome "$HAP1_GENOME" \
    --model_path "$MODEL_PATH" \
    --lineage "$LINEAGE" \
    --output "$H1_GFF" \
    --threads "$THREADS" \
    --batch_size 16                            

echo "[$(date)] Step 1 Complete. Extracting Hap1 translations..."
gffread "$H1_GFF" -g "$HAP1_GENOME" -y "$H1_FAA" 


# --- LOOP 2: Haplotype 2 (Hap2) Execution ---
echo "==================================================="
echo "[$(date)] STARTING STEP 2: ANNEVO Haplotype 2"
echo "==================================================="
H2_GFF="${OUT_DIR_H2}/annevo_raw_hap2.gff3"
H2_FAA="${OUT_DIR_H2}/annevo_raw_hap2_proteins.faa"

$EXEC_CMD \
    --genome "$HAP2_GENOME" \
    --model_path "$MODEL_PATH" \
    --lineage "$LINEAGE" \
    --output "$H2_GFF" \
    --threads "$THREADS" \
    --batch_size 16                            

echo "[$(date)] Step 2 Complete. Extracting Hap2 translations..."
gffread "$H2_GFF" -g "$HAP2_GENOME" -y "$H2_FAA" 


# ==============================================================================
#                         INTEGRATED RUN STATISTICS
# ==============================================================================
echo "==================================================="
echo "              ANNEVO SUMMARY METRICS               "
echo "==================================================="
echo "Haplotype 1 (Hap1) Output Location: ${OUT_DIR_H1}"
echo "  -> Raw Genes Predicted   : $(grep -c $'\tgene\t' "$H1_GFF")" 
echo "  -> Proteins Extracted    : $(grep -c "^>" "$H1_FAA")"     
echo ""
echo "Haplotype 2 (Hap2) Output Location: ${OUT_DIR_H2}"
echo "  -> Raw Genes Predicted   : $(grep -c $'\tgene\t' "$H2_GFF")" 
echo "  -> Proteins Extracted    : $(grep -c "^>" "$H2_FAA")"     
echo "Finished at: $(date)"                                        
echo "==================================================="