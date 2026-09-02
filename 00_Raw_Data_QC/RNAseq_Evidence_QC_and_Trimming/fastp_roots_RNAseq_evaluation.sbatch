#!/bin/bash

#SBATCH --job-name=fastp_qc_roots_by_id
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
# Sample 2962406 (R1_S87: L002, L003, L004)
# ---------------------------------------------------------
cat \
  "${IN_DIR}/2962406/2/2962406_R1_S87_L002_R1_001.fastq.gz" \
  "${IN_DIR}/2962406/3/2962406_R1_S87_L003_R1_001.fastq.gz" \
  "${IN_DIR}/2962406/4/2962406_R1_S87_L004_R1_001.fastq.gz" \
  > "${OUT_DIR}/2962406_merged_R1.fastq.gz"

cat \
  "${IN_DIR}/2962406/2/2962406_R1_S87_L002_R2_001.fastq.gz" \
  "${IN_DIR}/2962406/3/2962406_R1_S87_L003_R2_001.fastq.gz" \
  "${IN_DIR}/2962406/4/2962406_R1_S87_L004_R2_001.fastq.gz" \
  > "${OUT_DIR}/2962406_merged_R2.fastq.gz"

fastp \
  --in1 "${OUT_DIR}/2962406_merged_R1.fastq.gz" \
  --in2 "${OUT_DIR}/2962406_merged_R2.fastq.gz" \
  --disable_adapter_trimming \
  --disable_quality_filtering \
  --disable_length_filtering \
  --thread 8 \
  --html "${OUT_DIR}/2962406_raw_fastp_report.html" \
  --json "${OUT_DIR}/2962406_raw_fastp_report.json"

# ---------------------------------------------------------
# Sample 2962407 (R5_S99: L002, L003, L004)
# ---------------------------------------------------------
cat \
  "${IN_DIR}/2962407/2/2962407_R5_S99_L002_R1_001.fastq.gz" \
  "${IN_DIR}/2962407/3/2962407_R5_S99_L003_R1_001.fastq.gz" \
  "${IN_DIR}/2962407/4/2962407_R5_S99_L004_R1_001.fastq.gz" \
  > "${OUT_DIR}/2962407_merged_R1.fastq.gz"

cat \
  "${IN_DIR}/2962407/2/2962407_R5_S99_L002_R2_001.fastq.gz" \
  "${IN_DIR}/2962407/3/2962407_R5_S99_L003_R2_001.fastq.gz" \
  "${IN_DIR}/2962407/4/2962407_R5_S99_L004_R2_001.fastq.gz" \
  > "${OUT_DIR}/2962407_merged_R2.fastq.gz"

fastp \
  --in1 "${OUT_DIR}/2962407_merged_R1.fastq.gz" \
  --in2 "${OUT_DIR}/2962407_merged_R2.fastq.gz" \
  --disable_adapter_trimming \
  --disable_quality_filtering \
  --disable_length_filtering \
  --thread 8 \
  --html "${OUT_DIR}/2962407_raw_fastp_report.html" \
  --json "${OUT_DIR}/2962407_raw_fastp_report.json"

# ---------------------------------------------------------
# Sample 2962425 (R39_S125: L002, L003, L004)
# ---------------------------------------------------------
cat \
  "${IN_DIR}/2962425/2/2962425_R39_S125_L002_R1_001.fastq.gz" \
  "${IN_DIR}/2962425/3/2962425_R39_S125_L003_R1_001.fastq.gz" \
  "${IN_DIR}/2962425/4/2962425_R39_S125_L004_R1_001.fastq.gz" \
  > "${OUT_DIR}/2962425_merged_R1.fastq.gz"

cat \
  "${IN_DIR}/2962425/2/2962425_R39_S125_L002_R2_001.fastq.gz" \
  "${IN_DIR}/2962425/3/2962425_R39_S125_L003_R2_001.fastq.gz" \
  "${IN_DIR}/2962425/4/2962425_R39_S125_L004_R2_001.fastq.gz" \
  > "${OUT_DIR}/2962425_merged_R2.fastq.gz"

fastp \
  --in1 "${OUT_DIR}/2962425_merged_R1.fastq.gz" \
  --in2 "${OUT_DIR}/2962425_merged_R2.fastq.gz" \
  --disable_adapter_trimming \
  --disable_quality_filtering \
  --disable_length_filtering \
  --thread 8 \
  --html "${OUT_DIR}/2962425_raw_fastp_report.html" \
  --json "${OUT_DIR}/2962425_raw_fastp_report.json"

# ---------------------------------------------------------
# Sample 2962430 (R58_S146: L002, L003, L004)
# ---------------------------------------------------------
cat \
  "${IN_DIR}/2962430/2/2962430_R58_S146_L002_R1_001.fastq.gz" \
  "${IN_DIR}/2962430/3/2962430_R58_S146_L003_R1_001.fastq.gz" \
  "${IN_DIR}/2962430/4/2962430_R58_S146_L004_R1_001.fastq.gz" \
  > "${OUT_DIR}/2962430_merged_R1.fastq.gz"

cat \
  "${IN_DIR}/2962430/2/2962430_R58_S146_L002_R2_001.fastq.gz" \
  "${IN_DIR}/2962430/3/2962430_R58_S146_L003_R2_001.fastq.gz" \
  "${IN_DIR}/2962430/4/2962430_R58_S146_L004_R2_001.fastq.gz" \
  > "${OUT_DIR}/2962430_merged_R2.fastq.gz"

fastp \
  --in1 "${OUT_DIR}/2962430_merged_R1.fastq.gz" \
  --in2 "${OUT_DIR}/2962430_merged_R2.fastq.gz" \
  --disable_adapter_trimming \
  --disable_quality_filtering \
  --disable_length_filtering \
  --thread 8 \
  --html "${OUT_DIR}/2962430_raw_fastp_report.html" \
  --json "${OUT_DIR}/2962430_raw_fastp_report.json"

# ---------------------------------------------------------
# Sample 2962448 (Q71_S157: L002, L003, L004)
# ---------------------------------------------------------
cat \
  "${IN_DIR}/2962448/2/2962448_Q71_S157_L002_R1_001.fastq.gz" \
  "${IN_DIR}/2962448/3/2962448_Q71_S157_L003_R1_001.fastq.gz" \
  "${IN_DIR}/2962448/4/2962448_Q71_S157_L004_R1_001.fastq.gz" \
  > "${OUT_DIR}/2962448_merged_R1.fastq.gz"

cat \
  "${IN_DIR}/2962448/2/2962448_Q71_S157_L002_R2_001.fastq.gz" \
  "${IN_DIR}/2962448/3/2962448_Q71_S157_L003_R2_001.fastq.gz" \
  "${IN_DIR}/2962448/4/2962448_Q71_S157_L004_R2_001.fastq.gz" \
  > "${OUT_DIR}/2962448_merged_R2.fastq.gz"

fastp \
  --in1 "${OUT_DIR}/2962448_merged_R1.fastq.gz" \
  --in2 "${OUT_DIR}/2962448_merged_R2.fastq.gz" \
  --disable_adapter_trimming \
  --disable_quality_filtering \
  --disable_length_filtering \
  --thread 8 \
  --html "${OUT_DIR}/2962448_raw_fastp_report.html" \
  --json "${OUT_DIR}/2962448_raw_fastp_report.json"

# ---------------------------------------------------------
# Sample 2962460 (Q83_S140: L002, L003, L004)
# ---------------------------------------------------------
cat \
  "${IN_DIR}/2962460/2/2962460_Q83_S140_L002_R1_001.fastq.gz" \
  "${IN_DIR}/2962460/3/2962460_Q83_S140_L003_R1_001.fastq.gz" \
  "${IN_DIR}/2962460/4/2962460_Q83_S140_L004_R1_001.fastq.gz" \
  > "${OUT_DIR}/2962460_merged_R1.fastq.gz"

cat \
  "${IN_DIR}/2962460/2/2962460_Q83_S140_L002_R2_001.fastq.gz" \
  "${IN_DIR}/2962460/3/2962460_Q83_S140_L003_R2_001.fastq.gz" \
  "${IN_DIR}/2962460/4/2962460_Q83_S140_L004_R2_001.fastq.gz" \
  > "${OUT_DIR}/2962460_merged_R2.fastq.gz"

fastp \
  --in1 "${OUT_DIR}/2962460_merged_R1.fastq.gz" \
  --in2 "${OUT_DIR}/2962460_merged_R2.fastq.gz" \
  --disable_adapter_trimming \
  --disable_quality_filtering \
  --disable_length_filtering \
  --thread 8 \
  --html "${OUT_DIR}/2962460_raw_fastp_report.html" \
  --json "${OUT_DIR}/2962460_raw_fastp_report.json"