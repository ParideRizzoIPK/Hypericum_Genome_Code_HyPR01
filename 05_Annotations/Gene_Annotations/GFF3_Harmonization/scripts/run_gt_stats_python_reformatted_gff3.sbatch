#!/bin/bash
#SBATCH --job-name=GT_Stats_PythonReformatted
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2     # gt stat is generally efficient
#SBATCH --mem=8G              # Adjust if your GFF file is very large
#SBATCH --time=0-02:00:00     # 2 hours should be sufficient
#SBATCH --output=gt_stats_python_reformatted_hap2_%j.out
#SBATCH --error=gt_stats_python_reformatted_hap2_%j.err
#SBATCH --partition=cpu # IMPORTANT: Specify your Slurm partition

# Exit immediately if a command exits with a non-zero status
set -e

# --- Load necessary modules ---
echo "Loading GenomeTools module..."
module load genometools   # Or your cluster's specific module name for GenomeTools (gt)
echo "Modules loaded."

# --- Define File Paths ---

# 1. Path to the REFORMATTED GFF3 file (output from your Python script)
#    IMPORTANT: Update this path if your Python script saved the file with a different name or location.
REFORMATTED_GFF3_FILE="/path/to/input/directory/new2.gff3"

# 2. Output directory for GenomeTools statistics
GT_STATS_OUTPUT_DIR="/path/to/output/directory" # Or a new desired directory
mkdir -p "${GT_STATS_OUTPUT_DIR}"

# 3. Output files for GenomeTools statistics
STATS_SUMMARY_FILE="${GT_STATS_OUTPUT_DIR}/augustus_hap2.python_reformatted.gt_summary_statistics.txt"
GENE_LENGTH_DIST_FILE="${GT_STATS_OUTPUT_DIR}/augustus_hap2.python_reformatted.gene_length_distribution.txt"
EXON_LENGTH_DIST_FILE="${GT_STATS_OUTPUT_DIR}/augustus_hap2.python_reformatted.exon_length_distribution.txt"
INTRON_LENGTH_DIST_FILE="${GT_STATS_OUTPUT_DIR}/augustus_hap2.python_reformatted.intron_length_distribution.txt"

# --- Verify Input Reformatted GFF3 File Exists ---
echo "Checking for the reformatted GFF3 input file: ${REFORMATTED_GFF3_FILE}"
if [ ! -f "${REFORMATTED_GFF3_FILE}" ]; then
    echo "ERROR: Reformatted GFF3 file not found at ${REFORMATTED_GFF3_FILE}"
    echo "Please ensure the Python script ran successfully and created this file, and that the path is correct."
    exit 1
fi
if [ ! -s "${REFORMATTED_GFF3_FILE}" ]; then
    echo "ERROR: Reformatted GFF3 file (${REFORMATTED_GFF3_FILE}) is empty."
    exit 1
fi
# Quick check to see if it starts with the GFF3 version directive
GFF_HEADER_CHECK=$(head -n 1 "${REFORMATTED_GFF3_FILE}")
if [[ "${GFF_HEADER_CHECK}" != "##gff-version 3" ]]; then
    echo "WARNING: The reformatted GFF3 file (${REFORMATTED_GFF3_FILE}) does not start with '##gff-version 3'. Please verify its content."
    # You might choose to exit here if this is a strict requirement for your gt stat version
    # exit 1
fi
echo "Reformatted GFF3 file found and seems to have a valid header."

# --- Run GenomeTools (gt stat) ---
echo "------------------------------------"
echo "Running GenomeTools (gt stat) on the Python-reformatted GFF3 file..."
echo "Input GFF3: ${REFORMATTED_GFF3_FILE}"
echo "Output directory for statistics: ${GT_STATS_OUTPUT_DIR}"

# Basic summary statistics (counts of feature types, lengths, etc.)
# The -source augustus filter can be used if desired, but since the file should only contain AUGUSTUS features,
# it might not be strictly necessary. It won't hurt, though.
gt stat -force -o "${STATS_SUMMARY_FILE}" "${REFORMATTED_GFF3_FILE}" \
    || { echo "ERROR: gt stat (summary) failed. Please check the Slurm error log (${SLURM_JOB_ID}.err) and the reformatted GFF3 file: ${REFORMATTED_GFF3_FILE}"; exit 1; }

# Generate gene length distribution data
gt stat -force -genelengthdistri -o "${GENE_LENGTH_DIST_FILE}" "${REFORMATTED_GFF3_FILE}" \
    || { echo "ERROR: gt stat (genelengthdistri) failed."; exit 1; }

# Generate exon length distribution data
gt stat -force -exonlengthdistri -o "${EXON_LENGTH_DIST_FILE}" "${REFORMATTED_GFF3_FILE}" \
    || { echo "ERROR: gt stat (exonlengthdistri) failed."; exit 1; }

# Generate intron length distribution data
gt stat -force -intronlengthdistri -o "${INTRON_LENGTH_DIST_FILE}" "${REFORMATTED_GFF3_FILE}" \
    || { echo "ERROR: gt stat (intronlengthdistri) failed."; exit 1; }

echo "GenomeTools (gt stat) processing completed successfully."
echo "General summary statistics are in: ${STATS_SUMMARY_FILE}"
echo "Gene length distribution data: ${GENE_LENGTH_DIST_FILE}"
echo "Exon length distribution data: ${EXON_LENGTH_DIST_FILE}"
echo "Intron length distribution data: ${INTRON_LENGTH_DIST_FILE}"
echo "------------------------------------"
echo "Script finished."