#!/bin/bash
#SBATCH --job-name=QUAST_Eval
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/HyPR01_quast_hap2_%j.out
#SBATCH --error=logs/HyPR01_quast_hap2_%j.err
#SBATCH --partition=cpu

# --- Load necessary modules ---
module load gfatools
module load quast
module load python/3.12.0

# --- Define variables ---
GFA_INPUT_FILE="/path/to/your/directory/HyPR01_genome_hic_assembly.hic.hap2.p_ctg.gfa"
ASSEMBLY_FASTA_FILE="/path/to/your/directory/HyPR01_genome_hic_assembly.hic.hap2.p_ctg.fa"
ASSEMBLY_FILE="${ASSEMBLY_FASTA_FILE}"
OUTPUT_DIR="QUAST_results_HyPR01_hap2"
THREADS=8

# --- Step 1: Convert GFA to FASTA ---
echo "------------------------------------"
echo "Starting GFA to FASTA conversion..."
echo "Input GFA: ${GFA_INPUT_FILE}"
echo "Output FASTA: ${ASSEMBLY_FASTA_FILE}"
echo "------------------------------------"

if [ -f "${GFA_INPUT_FILE}" ]; then
    gfatools gfa2fa "${GFA_INPUT_FILE}" > "${ASSEMBLY_FASTA_FILE}"
    echo "GFA to FASTA conversion completed."
else
    echo "ERROR: Input GFA file not found: ${GFA_INPUT_FILE}"
    exit 1
fi

# Check if FASTA file was created and is not empty
if [ ! -s "${ASSEMBLY_FASTA_FILE}" ]; then
    echo "ERROR: FASTA file was not created or is empty after conversion: ${ASSEMBLY_FASTA_FILE}"
    exit 1
fi

# --- Run QUAST ---
echo "Starting QUAST evaluation..."
quast.py -o ${OUTPUT_DIR} -t ${THREADS} ${ASSEMBLY_FILE}

echo "QUAST evaluation finished. Results in ${OUTPUT_DIR}"
