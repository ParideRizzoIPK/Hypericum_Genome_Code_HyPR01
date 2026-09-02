#!/bin/bash

#SBATCH --job-name=fastp_qc_leaves_by_id
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --output=logs/job_fastp_%j.out
#SBATCH --error=logs/job_fastp_%j.err

# Load environment
. /etc/profile.d/modules.sh
module load fastp/1.0.1

# Define directories
IN_DIR="/path/to/input/directory"
OUT_DIR="/path/to/output/directory"
mkdir -p "${OUT_DIR}"

# ---------------------------------------------------------
# Sample 2778380 (S99)
# ---------------------------------------------------------
cat \
  "${IN_DIR}/2778380/1/2778380_221_S99_L001_R1_001.fastq.gz" \
  "${IN_DIR}/2778380/2/2778380_221_S99_L002_R1_001.fastq.gz" \
  > "${OUT_DIR}/2778380_merged_R1.fastq.gz"

cat \
  "${IN_DIR}/2778380/1/2778380_221_S99_L001_R2_001.fastq.gz" \
  "${IN_DIR}/2778380/2/2778380_221_S99_L002_R2_001.fastq.gz" \
  > "${OUT_DIR}/2778380_merged_R2.fastq.gz"

fastp \
  --in1 "${OUT_DIR}/2778380_merged_R1.fastq.gz" \
  --in2 "${OUT_DIR}/2778380_merged_R2.fastq.gz" \
  --disable_adapter_trimming \
  --disable_quality_filtering \
  --disable_length_filtering \
  --thread 8 \
  --html "${OUT_DIR}/2778380_raw_fastp_report.html" \
  --json "${OUT_DIR}/2778380_raw_fastp_report.json"

# ---------------------------------------------------------
# Sample 2778381 (S106)
# ---------------------------------------------------------
cat \
  "${IN_DIR}/2778381/1/2778381_222_S106_L001_R1_001.fastq.gz" \
  "${IN_DIR}/2778381/2/2778381_222_S106_L002_R1_001.fastq.gz" \
  > "${OUT_DIR}/2778381_merged_R1.fastq.gz"

cat \
  "${IN_DIR}/2778381/1/2778381_222_S106_L001_R2_001.fastq.gz" \
  "${IN_DIR}/2778381/2/2778381_222_S106_L002_R2_001.fastq.gz" \
  > "${OUT_DIR}/2778381_merged_R2.fastq.gz"

fastp \
  --in1 "${OUT_DIR}/2778381_merged_R1.fastq.gz" \
  --in2 "${OUT_DIR}/2778381_merged_R2.fastq.gz" \
  --disable_adapter_trimming \
  --disable_quality_filtering \
  --disable_length_filtering \
  --thread 8 \
  --html "${OUT_DIR}/2778381_raw_fastp_report.html" \
  --json "${OUT_DIR}/2778381_raw_fastp_report.json"

# ---------------------------------------------------------
# Sample 2778382 (S113)
# ---------------------------------------------------------
cat \
  "${IN_DIR}/2778382/1/2778382_223_S113_L001_R1_001.fastq.gz" \
  "${IN_DIR}/2778382/2/2778382_223_S113_L002_R1_001.fastq.gz" \
  > "${OUT_DIR}/2778382_merged_R1.fastq.gz"

cat \
  "${IN_DIR}/2778382/1/2778382_223_S113_L001_R2_001.fastq.gz" \
  "${IN_DIR}/2778382/2/2778382_223_S113_L002_R2_001.fastq.gz" \
  > "${OUT_DIR}/2778382_merged_R2.fastq.gz"

fastp \
  --in1 "${OUT_DIR}/2778382_merged_R1.fastq.gz" \
  --in2 "${OUT_DIR}/2778382_merged_R2.fastq.gz" \
  --disable_adapter_trimming \
  --disable_quality_filtering \
  --disable_length_filtering \
  --thread 8 \
  --html "${OUT_DIR}/2778382_raw_fastp_report.html" \
  --json "${OUT_DIR}/2778382_raw_fastp_report.json"

# ---------------------------------------------------------
# Sample 2778383 (S64)
# ---------------------------------------------------------
cat \
  "${IN_DIR}/2778383/1/2778383_221_E_S64_L001_R1_001.fastq.gz" \
  "${IN_DIR}/2778383/2/2778383_221_E_S64_L002_R1_001.fastq.gz" \
  > "${OUT_DIR}/2778383_merged_R1.fastq.gz"

cat \
  "${IN_DIR}/2778383/1/2778383_221_E_S64_L001_R2_001.fastq.gz" \
  "${IN_DIR}/2778383/2/2778383_221_E_S64_L002_R2_001.fastq.gz" \
  > "${OUT_DIR}/2778383_merged_R2.fastq.gz"

fastp \
  --in1 "${OUT_DIR}/2778383_merged_R1.fastq.gz" \
  --in2 "${OUT_DIR}/2778383_merged_R2.fastq.gz" \
  --disable_adapter_trimming \
  --disable_quality_filtering \
  --disable_length_filtering \
  --thread 8 \
  --html "${OUT_DIR}/2778383_raw_fastp_report.html" \
  --json "${OUT_DIR}/2778383_raw_fastp_report.json"