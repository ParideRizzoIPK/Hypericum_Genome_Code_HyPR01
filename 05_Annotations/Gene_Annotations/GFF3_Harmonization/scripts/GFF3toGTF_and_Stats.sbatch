#!/bin/bash
#SBATCH --job-name=GFF3toGTF_and_Stats
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2     # awk, gffread, and gt stat are generally efficient
#SBATCH --mem=10G
#SBATCH --time=0-03:00:00     # 3 hours should be sufficient
#SBATCH --output=gff3_to_gtf_stats_%j.out
#SBATCH --error=gff3_to_gtf_stats_%j.err
#SBATCH --partition=cpu # IMPORTANT: Specify your Slurm partition

# Exit immediately if a command exits with a non-zero status
set -e

# --- Load necessary modules ---
echo "Loading modules..."
module load gffread       # Or module for Cufflinks/StringTie which includes gffread
module load genometools   # Or module for GenomeTools (gt)
echo "Modules loaded."

# --- Define Input and Output File Paths (Using Absolute Paths) ---

# 1. Original GFF3 file from AUGUSTUS
ORIGINAL_AUGUSTUS_GFF3_FILE="/path/to/input/directory/augustus_hap2.gff3"

# 2. Directory where AUGUSTUS outputs are located and where intermediate/new files will be written
PROCESSING_DIR="/path/to/your/directory" # Using this as a base

# 3. Name for the intermediate, *filtered* GFF3 file (only AUGUSTUS predictions and ## directives)
FILTERED_GFF3_FILE="${PROCESSING_DIR}/augustus_hap2.AUGUSTUS_only.gff3"

# 4. Name for the output GTF file
OUTPUT_GTF_FILE="${PROCESSING_DIR}/augustus_hap2.AUGUSTUS_only.gtf"

# 5. GenomeTools (gt stat) output directory and files using the GTF
GT_STATS_OUTPUT_DIR="/path/to/output/directory" # New directory for these stats
mkdir -p "${GT_STATS_OUTPUT_DIR}"

STATS_SUMMARY_FILE_FROM_GTF="${GT_STATS_OUTPUT_DIR}/augustus_hap2.from_gtf.gt_summary_statistics.txt"
# Add other gt stat output files if needed, similar to previous script

# --- Verify Original Input GFF3 File Exists ---
echo "Checking for original AUGUSTUS GFF3 input file..."
if [ ! -f "${ORIGINAL_AUGUSTUS_GFF3_FILE}" ]; then
    echo "ERROR: Original AUGUSTUS GFF3 file not found at ${ORIGINAL_AUGUSTUS_GFF3_FILE}"
    exit 1
fi
echo "Original GFF3 input file found: ${ORIGINAL_AUGUSTUS_GFF3_FILE}"

# --- Step 1: Filter the original GFF3 file to keep only AUGUSTUS entries and ## GFF3 directives ---
echo "------------------------------------"
echo "Step 1: Filtering original GFF3 file..."
echo "Input GFF3: ${ORIGINAL_AUGUSTUS_GFF3_FILE}"
echo "Output Filtered GFF3: ${FILTERED_GFF3_FILE}"

# Revised awk command:
# - Prints lines starting with "##" (GFF3 directives).
# - OR prints lines that DO NOT start with "#" (i.e., are data/feature lines) AND where the 2nd column (source) is "AUGUSTUS".
awk 'BEGIN{OFS="\t"} $0 ~ /^##/ || ($1 !~ /^#/ && $2 == "AUGUSTUS")' "${ORIGINAL_AUGUSTUS_GFF3_FILE}" > "${FILTERED_GFF3_FILE}"

echo "GFF3 filtering complete. Filtered file created at ${FILTERED_GFF3_FILE}"
if [ ! -s "${FILTERED_GFF3_FILE}" ]; then
    echo "ERROR: Filtered GFF3 file (${FILTERED_GFF3_FILE}) is empty or was not created. Check original file and awk command."
    exit 1
fi
echo "------------------------------------"


# --- Step 2: Convert the Filtered GFF3 to GTF using gffread ---
echo ""
echo "------------------------------------"
echo "Step 2: Converting filtered GFF3 to GTF format..."
echo "Input Filtered GFF3: ${FILTERED_GFF3_FILE}"
echo "Output GTF: ${OUTPUT_GTF_FILE}"

# The -T option tells gffread to output in GTF format.
# gffread will attempt to create valid GTF attributes, including gene_id and transcript_id.
gffread "${FILTERED_GFF3_FILE}" -T -o "${OUTPUT_GTF_FILE}" \
    || { echo "ERROR: gffread conversion to GTF failed."; exit 1; }

echo "Conversion to GTF complete. Output GTF file: ${OUTPUT_GTF_FILE}"
if [ ! -s "${OUTPUT_GTF_FILE}" ]; then
    echo "ERROR: Output GTF file (${OUTPUT_GTF_FILE}) is empty or was not created."
    exit 1
fi
echo "------------------------------------"


# --- Step 3: Generate Annotation Statistics from the GTF file using GenomeTools (gt stat) ---
echo ""
echo "------------------------------------"
echo "Step 3: Running GenomeTools (gt stat) on the generated GTF file..."
echo "Input GTF: ${OUTPUT_GTF_FILE}"
echo "Output statistics file: ${STATS_SUMMARY_FILE_FROM_GTF}"
echo "------------------------------------"

# gt stat can also read GTF files. The output might be focused on features typically present in GTF
# (e.g., exons, CDS, start_codon, stop_codon associated with transcript_id and gene_id).
gt stat -force -o "${STATS_SUMMARY_FILE_FROM_GTF}" "${OUTPUT_GTF_FILE}" \
    || { echo "ERROR: gt stat failed on the GTF file. Please check error log and the GTF file: ${OUTPUT_GTF_FILE}"; exit 1; }

# You can add other gt stat options like -genelengthdistri etc. if desired,
# but they will be based on the features present in the GTF file.
# For example, GTF files don't typically have "gene" features, but rather exons/CDS grouped by transcript_id and gene_id.
# gt stat will report based on the feature types it finds in the GTF.

echo "GenomeTools (gt stat) processing on GTF completed."
echo "Statistics from GTF are in: ${STATS_SUMMARY_FILE_FROM_GTF}"
echo "------------------------------------"
echo "Script finished."