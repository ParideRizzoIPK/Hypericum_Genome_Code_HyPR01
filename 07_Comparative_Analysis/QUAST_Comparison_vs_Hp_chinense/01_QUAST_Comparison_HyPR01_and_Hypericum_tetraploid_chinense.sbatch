#!/bin/bash
#
#SBATCH --job-name=HyPR01_Genome_QUAST
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=5:00:00
#SBATCH --output=logs/quast_%j.out
#SBATCH --error=logs/quast_%j.err

# Strict error handling: exit on command failure (-e), unset variables (-u), pipeline failure (-o pipefail)
set -euo pipefail

# ==============================================================================
# 1. ENVIRONMENT & PATH SETUP
# ==============================================================================
export SINGULARITY_BIND="/path/to/your/directory"

# Load required cluster modules
module load quast/5.2.0
module load samtools/1.23.1

# Base & Output Directories
BASE_DIR="/path/to/your/directory"
OUT_DIR="${BASE_DIR}/HyPR01_Genome_QUAST"

# Working Subdirectories
PREP_DIR="${OUT_DIR}/prepared_fastas"
SPLIT_DIR="${BASE_DIR}/Hypericum_perforatum_subsp_chinense/HyPR01_synteny_v2/split"

# Raw Genome Input Paths (with ChrUn)
RAW_HAP1_GENOME="${BASE_DIR}/hap1.masked_ChrUn.fa"
RAW_HAP2_GENOME="${BASE_DIR}/hap2.masked_ChrUn.fa"

# Threads
THREADS=16

# Initialize directories
mkdir -p "${OUT_DIR}"
mkdir -p "${PREP_DIR}"

# ==============================================================================
# 2. INPUT VALIDATION
# ==============================================================================
echo "[$(date)] ==== STEP 1: VALIDATING INPUT FILES ===="

# Check raw haplotype files
for file in "${RAW_HAP1_GENOME}" "${RAW_HAP2_GENOME}"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required raw haplotype FASTA does not exist: $file" >&2
        exit 1
    fi
done

# Check 32 split chromosome files for subgenomes A, B, C, D
for S in A B C D; do
    for CHR in {1..8}; do
        FILE="${SPLIT_DIR}/ref_Chr${CHR}${S}.fa"
        if [ ! -f "$FILE" ]; then
            echo "ERROR: Missing split chromosome file required for subgenome ${S}: $FILE" >&2
            exit 1
        fi
    done
done

echo "All required input files successfully verified."

# ==============================================================================
# 3. PRE-PROCESSING: FILTERING (noChrUn) & CONCATENATION (Subgenomes)
# ==============================================================================
echo -e "\n[$(date)] ==== STEP 2: PRE-PROCESSING FASTAs ===="

# A. Extract Nuclear Chromosomes (Chr01 - Chr08) for Hap1 and Hap2
HAP1_NO_CHRUN="${PREP_DIR}/hap1_noChrUn.fa"
HAP2_NO_CHRUN="${PREP_DIR}/hap2_noChrUn.fa"

extract_nuclear_chromosomes () {
    local input_fa="$1"
    local output_fa="$2"
    
    echo "Indexing ${input_fa}..."
    samtools faidx "${input_fa}"
    
    echo "Extracting Chr01..Chr08 into ${output_fa}..."
    samtools faidx "${input_fa}" Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 > "${output_fa}"
}

echo "[$(date)] Extracting nuclear chromosomes for Hap1 (noChrUn)..."
extract_nuclear_chromosomes "${RAW_HAP1_GENOME}" "${HAP1_NO_CHRUN}"

echo "[$(date)] Extracting nuclear chromosomes for Hap2 (noChrUn)..."
extract_nuclear_chromosomes "${RAW_HAP2_GENOME}" "${HAP2_NO_CHRUN}"

# B. Concatenate Split Subgenomes A, B, C, D
SUBA_FA="${PREP_DIR}/subgenome_A.fa"
SUBB_FA="${PREP_DIR}/subgenome_B.fa"
SUBC_FA="${PREP_DIR}/subgenome_C.fa"
SUBD_FA="${PREP_DIR}/subgenome_D.fa"

for S in A B C D; do
    TARGET_FA="${PREP_DIR}/subgenome_${S}.fa"
    echo "[$(date)] Concatenating chromosomes 1-8 into subgenome_${S}.fa..."
    > "${TARGET_FA}"
    for CHR in {1..8}; do
        cat "${SPLIT_DIR}/ref_Chr${CHR}${S}.fa" >> "${TARGET_FA}"
    done
done

# ==============================================================================
# 4. RUN COMPREHENSIVE MULTI-ASSEMBLY QUAST EVALUATION
# ==============================================================================
echo -e "\n[$(date)] ==== STEP 3: EXECUTING COMPREHENSIVE QUAST ANALYSIS ===="

QUAST_OUT="${OUT_DIR}/multi_assembly_comparison"

# Labels for QUAST report columns
LABELS="Hap1_with_ChrUn,Hap2_with_ChrUn,Hap1_noChrUn,Hap2_noChrUn,Subgenome_A,Subgenome_B,Subgenome_C,Subgenome_D"

echo "----------------------------------------------------------------------"
echo "[$(date)] Running side-by-side QUAST comparison"
echo "Output Directory: ${QUAST_OUT}"
echo "Threads:           ${THREADS}"
echo "Assemblies:        ${LABELS}"
echo "----------------------------------------------------------------------"

quast.py \
  -o "${QUAST_OUT}" \
  -t "${THREADS}" \
  -l "${LABELS}" \
  --large \
  "${RAW_HAP1_GENOME}" \
  "${RAW_HAP2_GENOME}" \
  "${HAP1_NO_CHRUN}" \
  "${HAP2_NO_CHRUN}" \
  "${SUBA_FA}" \
  "${SUBB_FA}" \
  "${SUBC_FA}" \
  "${SUBD_FA}"

echo -e "\n[$(date)] ==== QUAST EVALUATION COMPLETED SUCCESSFULLY ===="
echo "Comparative summary reports are available at:"
echo "  - HTML Report: ${QUAST_OUT}/report.html"
echo "  - Text Summary: ${QUAST_OUT}/report.txt"
echo "  - PDF Report:  ${QUAST_OUT}/report.pdf"