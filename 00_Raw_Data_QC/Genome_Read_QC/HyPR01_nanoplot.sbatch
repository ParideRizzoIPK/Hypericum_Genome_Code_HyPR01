#!/bin/bash
#SBATCH --job-name=HyPR01_NanoPlot_QC
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/nanoplot_run_%j.out
#SBATCH --error=logs/nanoplot_run_%j.err

# ────────────────────────────────────────────────────────
# SCRIPT METADATA & VERSION CONTROL
# Tool Version: NanoPlot (Installed via micromamba)
# Last Modified: 2026-07-15
# Description: Generates quality control statistics, plots,
#              and HTML summaries from PacBio CCS consensus.
# ────────────────────────────────────────────────────────

set -e

# --- MASTER WORKING DIRECTORIES & FILES ---
INPUT_FASTQ="/path/to/input/directory/hypr01_hifi_reads.fastq"
OUTPUT_DIR="/path/to/output/directory"

# 1. Create output directory dynamically before writing logs
mkdir -p "$OUTPUT_DIR"

echo "[$(date)] Pipeline initiated."
echo "[INFO] Input File: $INPUT_FASTQ"
echo "[INFO] Output Directory: $OUTPUT_DIR"

# 2. Load Systems Environment & Activate Micromamba
echo "[$(date)] Initializing micromamba shell hook..."

# Source your personal profile configuration
if [ -f "$HOME/.bashrc" ]; then
    source "$HOME/.bashrc"
fi

# Hardcode the exact root prefix on the filer mount
export MAMBA_ROOT_PREFIX="/path/to/your/directory/micromamba"

# Target the micromamba binary at its true path
MAMBA_BIN="/path/to/your/directory/micromamba/bin/micromamba"

# Evaluate shell hook using the explicit path
eval "$($MAMBA_BIN shell hook --shell=bash)"

echo "[$(date)] Activating 'nanoplot' environment..."
micromamba activate nanoplot

# Verify installation in log
echo "[INFO] NanoPlot Version:"
NanoPlot --version

# 3. Execute NanoPlot Analysis
echo "[$(date)] Launching NanoPlot processing loop..."
NanoPlot \
  --threads 8 \
  --fastq "$INPUT_FASTQ" \
  --outdir "$OUTPUT_DIR" \
  --prefix HyPR01_PacBio_QC_

echo "[$(date)] NanoPlot execution completed successfully."