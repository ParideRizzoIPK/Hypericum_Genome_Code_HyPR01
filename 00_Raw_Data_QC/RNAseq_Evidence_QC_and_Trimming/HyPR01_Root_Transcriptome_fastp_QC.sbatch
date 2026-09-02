#!/bin/bash
#SBATCH --job-name=fastp_roots
#SBATCH --output=logs/fastp_roots_%j.out
#SBATCH --error=logs/fastp_roots_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00

# Exit immediately if a command exits with a non-zero status
set -euo pipefail

echo "Starting fastp trimming pipeline at $(date)"

# ==========================================
# 1. DEFINE PATHS & VARIABLES
# ==========================================
ROOT_DIR="/path/to/your/directory"
OUTPUT_DIR="${ROOT_DIR}/trimmed_reads"

# The SLURM_CPUS_PER_TASK variable automatically matches your --cpus-per-task setting (16)
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# ==========================================
# 2. SETUP DIRECTORIES & MODULES
# ==========================================
echo "Creating output directory at $OUTPUT_DIR..."
mkdir -p "$OUTPUT_DIR"

module purge
# Loading the newest available version
module load fastp/1.0.1

# ==========================================
# 3. RUN FASTP ON ALL SAMPLES
# ==========================================
echo "Searching for fastq files and starting trimming..."

# Find all R1 files in the subdirectories
find "$ROOT_DIR" -type f -name "*_R1_001.fastq.gz" | while read -r R1_FILE; do
    
    # Automatically determine the matching R2 file name by replacing "R1" with "R2"
    R2_FILE="${R1_FILE/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    
    # Extract just the base name of the file (e.g., "2962406_R1_S87_L002") to use for output naming
    BASENAME=$(basename "$R1_FILE" "_R1_001.fastq.gz")
    
    # Define the output files
    OUT_R1="${OUTPUT_DIR}/${BASENAME}_trimmed_R1.fastq.gz"
    OUT_R2="${OUTPUT_DIR}/${BASENAME}_trimmed_R2.fastq.gz"
    HTML_REPORT="${OUTPUT_DIR}/${BASENAME}_fastp.html"
    JSON_REPORT="${OUTPUT_DIR}/${BASENAME}_fastp.json"
    
    # Check if the output already exists so it doesn't repeat work if restarted
    if [ ! -f "$OUT_R1" ]; then
        echo "------------------------------------------------"
        echo "Processing: $BASENAME"
        
        fastp \
            -i "$R1_FILE" \
            -I "$R2_FILE" \
            -o "$OUT_R1" \
            -O "$OUT_R2" \
            --html "$HTML_REPORT" \
            --json "$JSON_REPORT" \
            --thread "$THREADS"
    else
        echo "Skipping $BASENAME: Trimmed files already exist."
    fi
done

echo "------------------------------------------------"
echo "All samples successfully trimmed at $(date)"