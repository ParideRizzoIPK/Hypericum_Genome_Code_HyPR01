#!/bin/bash
#
#SBATCH --job-name=HyPR01_QUAST_Scaffolds
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=5:00:00
#SBATCH --output=logs/quast_scaffolds_%j.out
#SBATCH --error=logs/quast_scaffolds_%j.err

# Strict error handling: exit on command failure (-e), unset variables (-u), pipeline failure (-o pipefail)
set -euo pipefail

# ==============================================================================
# 1. ENVIRONMENT & PATH SETUP
# ==============================================================================
export SINGULARITY_BIND="/path/to/your/directory"

# Load required module
module load quast/5.2.0

# Base & Main Output Directories
BASE_DIR="/path/to/your/directory"
OUT_DIR="${BASE_DIR}/HyPR01_Genome_QUAST"

# Directory containing pre-processed FASTAs from the previous run
PREP_DIR="${OUT_DIR}/prepared_fastas"

# NEW Output Subdirectory to keep previous 'multi_assembly_comparison' dataset intact
QUAST_OUT="${OUT_DIR}/multi_assembly_comparison_scaffolds"

# Raw Genome Input Paths (with ChrUn)
RAW_HAP1_GENOME="${BASE_DIR}/hap1.masked_ChrUn.fa"
RAW_HAP2_GENOME="${BASE_DIR}/hap2.masked_ChrUn.fa"

# Prepared Input Paths (reused directly from previous run)
HAP1_NO_CHRUN="${PREP_DIR}/hap1_noChrUn.fa"
HAP2_NO_CHRUN="${PREP_DIR}/hap2_noChrUn.fa"
SUBA_FA="${PREP_DIR}/subgenome_A.fa"
SUBB_FA="${PREP_DIR}/subgenome_B.fa"
SUBC_FA="${PREP_DIR}/subgenome_C.fa"
SUBD_FA="${PREP_DIR}/subgenome_D.fa"

# Threads
THREADS=16

# Create the new output subdirectory
mkdir -p "${QUAST_OUT}"

# ==============================================================================
# 2. INPUT VALIDATION
# ==============================================================================
echo "[$(date)] ==== STEP 1: VALIDATING PREVIOUSLY GENERATED FASTAs ===="

INPUT_FILES=(
  "${RAW_HAP1_GENOME}"
  "${RAW_HAP2_GENOME}"
  "${HAP1_NO_CHRUN}"
  "${HAP2_NO_CHRUN}"
  "${SUBA_FA}"
  "${SUBB_FA}"
  "${SUBC_FA}"
  "${SUBD_FA}"
)

# Verify all 8 FASTA files exist before beginning
for file in "${INPUT_FILES[@]}"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required input FASTA file does not exist: $file" >&2
        exit 1
    fi
    echo "  [OK] Verified: $file"
done

echo "All 8 input FASTA files successfully verified."

# ==============================================================================
# 3. RUN MULTI-ASSEMBLY QUAST EVALUATION WITH --split-scaffolds
# ==============================================================================
echo -e "\n[$(date)] ==== STEP 2: EXECUTING QUAST EVALUATION (--split-scaffolds) ===="

# Labels for QUAST report columns
LABELS="Hap1_with_ChrUn,Hap2_with_ChrUn,Hap1_noChrUn,Hap2_noChrUn,Subgenome_A,Subgenome_B,Subgenome_C,Subgenome_D"

echo "----------------------------------------------------------------------"
echo "[$(date)] Running side-by-side QUAST comparison"
echo "Output Directory: ${QUAST_OUT}"
echo "Threads:           ${THREADS}"
echo "Assemblies:        ${LABELS}"
echo "----------------------------------------------------------------------"

# Corrected parameter: --split-scaffolds splits scaffolds at runs of Ns >= 10
quast.py \
  -o "${QUAST_OUT}" \
  -t "${THREADS}" \
  -l "${LABELS}" \
  --large \
  --split-scaffolds \
  "${RAW_HAP1_GENOME}" \
  "${RAW_HAP2_GENOME}" \
  "${HAP1_NO_CHRUN}" \
  "${HAP2_NO_CHRUN}" \
  "${SUBA_FA}" \
  "${SUBB_FA}" \
  "${SUBC_FA}" \
  "${SUBD_FA}"

echo -e "\n[$(date)] ==== QUAST EVALUATION COMPLETED SUCCESSFULLY ===="
echo "New results saved independently without affecting prior runs:"
echo "  - HTML Report: ${QUAST_OUT}/report.html"
echo "  - Text Summary: ${QUAST_OUT}/report.txt"
echo "  - PDF Report:  ${QUAST_OUT}/report.pdf"