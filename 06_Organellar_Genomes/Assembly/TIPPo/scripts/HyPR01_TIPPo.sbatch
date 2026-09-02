#!/bin/bash

# --- SBATCH Directives for SLURM Cluster ---
#SBATCH --job-name=TIPPo_HyPR01_organelle_final
#SBATCH --auks=yes
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=600G
#SBATCH --nodes=1
#SBATCH --output=logs/tippo_organelle_final_%j.out
#SBATCH --error=logs/tippo_organelle_final_%j.err
#SBATCH --time=13-00:00:00

# --- Script Best Practices ---
set -e -u -o pipefail

# --- Environment Setup ---
echo "INFO: Setting up module environment..."
. /etc/profile.d/modules.sh
module purge
module switch TIPP TIPP/1.3.0
# IMPORTANT: Bind paths for your home directory and the raw reads directory
export SINGULARITY_BIND="/path/to/your/directory,/path/to/your/other/directory"
echo "INFO: Environment setup complete."

# --- Define File and Directory Paths ---
READS_FILE_ORIGINAL="/path/to/input/directory/hypr01_hifi_reads.fastq"
BASE_OUTPUT_DIR="/path/to/output/directory"
ORGANELLE_OUTPUT_DIR="${BASE_OUTPUT_DIR}/organelle_assembly_final"
LOG_DIR="${BASE_OUTPUT_DIR}/slurm_logs"

READS_FILE_DOWNSAMPLED="${ORGANELLE_OUTPUT_DIR}/HyPR01_downsampled_1M_reads.fastq"

# --- Pre-run Setup ---
echo "INFO: Creating output directories..."
mkdir -p "${ORGANELLE_OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"

# --- STEP 1: Downsample Reads to 1 Million Reads (Memory-Efficiently) ---
echo "==============================================="
echo "INFO: Starting memory-efficient downsampling of reads to 1 million..."

# Use 'cat' for uncompressed .fastq files.
# Sample down to 1,000,000 reads
cat "${READS_FILE_ORIGINAL}" | tipp seqtk sample -s100 - 1000000 > "${READS_FILE_DOWNSAMPLED}"

echo "INFO: Downsampling complete. New reads file is ${READS_FILE_DOWNSAMPLED}"
echo "==============================================="

# --- STEP 2: Run TIPPo Assembly on the Downsampled Reads ---
echo "INFO: Starting DE NOVO mitochondrial assembly on downsampled reads..."
cd "${ORGANELLE_OUTPUT_DIR}"

tipp TIPPo.v2.4.pl \
  -f "${READS_FILE_DOWNSAMPLED}" \
  -g organelle \
  -p hifi \
  -t "${SLURM_CPUS_PER_TASK}"

echo "INFO: Mitochondrial assembly command finished."
echo "==============================================="
echo "Script finished. Results are in ${ORGANELLE_OUTPUT_DIR}"