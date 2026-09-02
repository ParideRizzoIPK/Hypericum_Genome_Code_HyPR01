#!/bin/bash
#SBATCH --job-name=GFA2FA_BUSCO_Eval
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16      # Adjust based on BUSCO documentation and available resources
#SBATCH --mem=256G                 # Adjust based on genome size and BUSCO needs (conversion also uses some memory)
#SBATCH --time=2-00:30:00         # 2 days and 30 minutes (added a bit for conversion)
#SBATCH --output=gfa2fa_busco_%j.out
#SBATCH --error=gfa2fa_busco_%j.err
#SBATCH --partition=cpu # Specify your Slurm partition

# Exit script on any error
set -e

# --- Load necessary modules ---
module load gfatools  # Or the specific module name for gfatools on your cluster
module load BUSCO/5.2.2     # Or the specific module name for BUSCO
module load augustus/3.3.1  # BUSCO often uses Augustus for gene prediction (ensure it's compatible)
#module load python    # BUSCO Python dependencies

# --- Define variables ---
GFA_INPUT_FILE="/path/to/your/directory"
ASSEMBLY_FASTA_FILE="/path/to/input/directory/HyPR01_genome_hic_assembly.hic.bothhap.p_ctg.fa" # Output of conversion, input for BUSCO

BUSCO_OUTPUT_NAME="BUSCO_HyPR01_bothhap"
THREADS=16


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

# --- Step 2: Run BUSCO ---
echo "------------------------------------"
echo "Starting BUSCO evaluation with auto-lineage detection..."
echo "Input FASTA: ${ASSEMBLY_FASTA_FILE}"
echo "Output Name: ${BUSCO_OUTPUT_NAME}"
echo "Threads: ${THREADS}"
echo "------------------------------------"

busco -i "${ASSEMBLY_FASTA_FILE}" \
      -o "${BUSCO_OUTPUT_NAME}" \
      -m genome \
      -c "${THREADS}" \
      --auto-lineage

echo "BUSCO evaluation finished."
echo "Results are in the directory: ${BUSCO_OUTPUT_NAME}/"
echo "------------------------------------"
echo "Job completed."