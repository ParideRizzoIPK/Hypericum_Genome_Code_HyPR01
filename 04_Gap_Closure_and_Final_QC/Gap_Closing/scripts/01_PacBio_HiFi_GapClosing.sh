#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G            # High memory allocation matching your 33.4 GB long-read dataset criteria
#SBATCH --time=24:00:00
#SBATCH --job-name=tgs_gapcloser
#SBATCH --output=logs/gapclose_%j.out
#SBATCH --error=logs/gapclose_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  

echo "#########################################################"
echo "Pipeline: Step 5 - Whole-Genome HiFi Gap Closing (Workspace Isolated)"
echo "Date: $(date)"
echo "#########################################################"

# --- 1. ENVIRONMENT & MODULE INITIALIZATION ---
source ~/.bashrc
set -u            

# Load required cluster modules
module load samtools/1.23.1
module load TGS-GapCloser/1.2.1

# Container Bind Fixes 
export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# --- 2. AUTOMATED EXECUTABLE NAME DETECTION ---
if command -v TGS-GapCloser.sh &> /dev/null; then
    GAPCLOSER_CMD="TGS-GapCloser.sh"
elif command -v tgsgapcloser &> /dev/null; then
    GAPCLOSER_CMD="tgsgapcloser"
else
    echo "CRITICAL ERROR: Neither TGS-GapCloser.sh nor tgsgapcloser was found after loading the module."
    exit 1
fi
echo "Using verified executable alias: ${GAPCLOSER_CMD}"

# --- 3. CONFIGURATION AND DIRECTORIES ---
OUTPUT_DIR="/path/to/output/directory"
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

# Raw Input Data Vectors
HIC_FASTQ_READS="/path/to/input/directory/hypr01_hifi_reads.fastq"
HAP1_INPUT_FA="/path/to/input/directory/Hap1_curated_chr_un.fasta"
HAP2_INPUT_FA="/path/to/input/directory/Hap2_curated_chr_un_RENAMED.fasta"

THREADS=${SLURM_CPUS_PER_TASK}
FASTA_READS="${OUTPUT_DIR}/hypr01_hifi_reads.fasta"

# ==============================================================================
# STEP 5.1: MEMORY-SAFE FASTQ TO FASTA CONVERSION STREAM
# ==============================================================================
if [ ! -f "${FASTA_READS}" ]; then
    echo ">>>> Streaming 33.4 GB FASTQ data to valid FASTA footprint... <<<<"
    awk 'NR%4==1{sub(/^@/,">");print;next} NR%4==2{print}' "${HIC_FASTQ_READS}" > "${FASTA_READS}"
    echo "  -> FASTA long-read conversion complete."
else
    echo ">>>> Target long-read FASTA found. Skipping conversion step. <<<<"
fi

# ==============================================================================
# STEP 5.2: EXECUTE GAP-CLOSING FOR HAPLOTYPE 1 (ISOLATED WORKSPACE)
# ==============================================================================
echo ">>>> PROCESSING HAPLOTYPE 1 <<<<"
mkdir -p "${OUTPUT_DIR}/Hap1_Workspace"
cd "${OUTPUT_DIR}/Hap1_Workspace"

echo "  -> Isolating Haplotype 1 chromosomal sequences from UnChr..."
awk '/^>UnChr/{p=0;next} /^>/{p=1} p' "${HAP1_INPUT_FA}" > hap1_chroms_only.fasta
awk '/^>UnChr/{p=1;print;next} /^>/{p=0} p' "${HAP1_INPUT_FA}" > hap1_unchr_only.fasta

echo "  -> Running ${GAPCLOSER_CMD} on Haplotype 1 Chromosomes..."
${GAPCLOSER_CMD} \
    --scaff hap1_chroms_only.fasta \
    --reads "${FASTA_READS}" \
    --output hap1_gapclosed_run \
    --tgstype pb \
    --ne \
    --thread "${THREADS}"

echo "  -> Recombining Haplotype 1 gap-closed chromosomes with intact UnChr..."
cat hap1_gapclosed_run.scaff_seqs hap1_unchr_only.fasta > "${OUTPUT_DIR}/Hap1_curated_chr_un_gapclosed.fasta"


# ==============================================================================
# STEP 5.3: EXECUTE GAP-CLOSING FOR HAPLOTYPE 2 (ISOLATED WORKSPACE)
# ==============================================================================
echo ">>>> PROCESSING HAPLOTYPE 2 <<<<"
mkdir -p "${OUTPUT_DIR}/Hap2_Workspace"
cd "${OUTPUT_DIR}/Hap2_Workspace"

echo "  -> Isolating Haplotype 2 chromosomal sequences from UnChr..."
awk '/^>UnChr/{p=0;next} /^>/{p=1} p' "${HAP2_INPUT_FA}" > hap2_chroms_only.fasta
awk '/^>UnChr/{p=1;print;next} /^>/{p=0} p' "${HAP2_INPUT_FA}" > hap2_unchr_only.fasta

echo "  -> Running ${GAPCLOSER_CMD} on Haplotype 2 Chromosomes..."
# Checkpoints will generate locally here, completely blind to Hap1's previous metrics
${GAPCLOSER_CMD} \
    --scaff hap2_chroms_only.fasta \
    --reads "${FASTA_READS}" \
    --output hap2_gapclosed_run \
    --tgstype pb \
    --ne \
    --thread "${THREADS}"

echo "  -> Recombining Haplotype 2 gap-closed chromosomes with intact UnChr..."
cat hap2_gapclosed_run.scaff_seqs hap2_unchr_only.fasta > "${OUTPUT_DIR}/Hap2_curated_chr_un_gapclosed.fasta"


# ==============================================================================
# STEP 5.4: POST-PROCESSING WORKSPACE CLEANUP
# ==============================================================================
echo ">>>> Running workspace sanitation and cleaning temporary sandboxes... <<<<"
cd "${OUTPUT_DIR}"
rm -rf "${OUTPUT_DIR}/Hap1_Workspace" "${OUTPUT_DIR}/Hap2_Workspace"

echo "#########################################################"
echo "GAP CLOSING COMPLETED SUCCESSFULLY."
echo "Final Curation Targets Generated:"
echo "  - Hap1 Final: ${OUTPUT_DIR}/Hap1_curated_chr_un_gapclosed.fasta"
echo "  - Hap2 Final: ${OUTPUT_DIR}/Hap2_curated_chr_un_gapclosed.fasta"
echo "#########################################################"