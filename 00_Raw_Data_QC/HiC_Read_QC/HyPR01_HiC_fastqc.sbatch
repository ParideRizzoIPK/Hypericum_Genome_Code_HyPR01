#!/bin/bash
#
#SBATCH --job-name=fastqc_HyPR01_HiC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/HyPR01_HiC_fastqc_%x_%j.out
#SBATCH --error=logs/HyPR01_HiC_fastqc_%x_%j.err

## -- ENVIRONMENT SETUP --
export SINGULARITY_BIND="/path/to/your/directory"
module load fastqc/0.12.1

THREADS=2

# Directories
DATA_DIR="/path/to/your/directory"
OUT_DIR="${DATA_DIR}/HyPR01_HiC_fastqc_output"

mkdir -p "${OUT_DIR}"

# Input FastQ files
R1="${DATA_DIR}/2037511_Hypericum_PR1_HiC166_S68_L002_R1_001.fastq.gz"
R2="${DATA_DIR}/2037511_Hypericum_PR1_HiC166_S68_L002_R2_001.fastq.gz"

echo "[$(date)] ==== STARTING FASTQC ANALYSIS ===="
echo "[$(date)] Running FastQC on HiC paired-end reads..."

fastqc -t "${THREADS}" -o "${OUT_DIR}" "${R1}" "${R2}"

echo "[$(date)] ==== ALL TASKS COMPLETED ===="
