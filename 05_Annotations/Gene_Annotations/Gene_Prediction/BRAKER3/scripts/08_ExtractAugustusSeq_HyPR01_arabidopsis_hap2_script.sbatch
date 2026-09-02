#!/bin/bash
#SBATCH --job-name=ExtractAugustusSeqs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2     # awk and gffread are generally fast and efficient
#SBATCH --mem=8G              # Sufficient for filtering and gffread
#SBATCH --time=0-02:00:00     # 2 hours, typically much faster
#SBATCH --output=extract_augustus_HyPR01_hap2_arabidopsis_sequences_%j.out
#SBATCH --error=extract_augustus_HyPR01_hap2_arabidopsis_sequences_%j.err
#SBATCH --partition=cpu # IMPORTANT: Specify your Slurm partition

# Exit immediately if a command exits with a non-zero status
set -e

# --- Load necessary modules ---
# gffread is often bundled with Cufflinks or StringTie, or available as a standalone module.
echo "Loading modules..."
module load gffread       # Or: module load cufflinks; or: module load stringtie;
echo "Modules loaded."

# --- Define Input and Output File Paths (Using Absolute Paths from your info) ---

# 1. Original GFF3 file from AUGUSTUS (this is your primary input)
ORIGINAL_AUGUSTUS_GFF3_FILE="/path/to/input/directory/augustus_hap2.gff3"

# 2. Directory where AUGUSTUS outputs are located and where the filtered GFF3 and new sequence files will be saved
AUGUSTUS_OUTPUT_DIR="/path/to/output/directory"

# 3. Name for the intermediate, *filtered* GFF3 file.
#    This file will be created in the AUGUSTUS_OUTPUT_DIR and will contain only AUGUSTUS gene predictions and the ##gff-version 3 header.
FILTERED_GFF3_FOR_PROCESSING="${AUGUSTUS_OUTPUT_DIR}/augustus_hap2.AUGUSTUS_only_for_seq_extraction.gff3"

# 4. Path to the soft-masked genome FASTA file (required by gffread as the reference sequence)
MASKED_GENOME_FASTA="/path/to/input/directory"

# 5. Desired output FASTA files for proteins, CDS, and transcripts (will be saved in AUGUSTUS_OUTPUT_DIR)
PROTEIN_FAA_OUTPUT="${AUGUSTUS_OUTPUT_DIR}/augustus_hap2.protein.faa"
CDS_FNA_OUTPUT="${AUGUSTUS_OUTPUT_DIR}/augustus_hap2.cds.fna"
TRANSCRIPT_FNA_OUTPUT="${AUGUSTUS_OUTPUT_DIR}/augustus_hap2.transcripts.fna"

# --- Verify Original Input Files Exist ---
echo "Checking for necessary input files..."
if [ ! -f "${ORIGINAL_AUGUSTUS_GFF3_FILE}" ]; then
    echo "ERROR: Original AUGUSTUS GFF3 file not found at ${ORIGINAL_AUGUSTUS_GFF3_FILE}"
    exit 1
fi
if [ ! -f "${MASKED_GENOME_FASTA}" ]; then
    echo "ERROR: Masked genome FASTA file not found at ${MASKED_GENOME_FASTA}"
    exit 1
fi
echo "Input files verified."

# --- Step 1: Filter the original GFF3 file to keep only AUGUSTUS entries and ## GFF3 directives ---
echo "------------------------------------"
echo "Step 1: Filtering original GFF3 file specifically for gffread..."
echo "Input GFF3: ${ORIGINAL_AUGUSTUS_GFF3_FILE}"
echo "Output Filtered GFF3 for sequence extraction: ${FILTERED_GFF3_FOR_PROCESSING}"

# This awk command prints:
#   - Lines starting with "##" (GFF3 directives like ##gff-version 3)
#   - OR lines where the second column (source field) is exactly "AUGUSTUS".
# This effectively skips other comment lines (single #) and the 'softmask' lines from your original file.
awk 'BEGIN{OFS="\t"} $1 ~ /^##/ || $2 == "AUGUSTUS"' "${ORIGINAL_AUGUSTUS_GFF3_FILE}" > "${FILTERED_GFF3_FOR_PROCESSING}"

echo "GFF3 filtering complete. Filtered file for sequence extraction created at ${FILTERED_GFF3_FOR_PROCESSING}"

# Basic check to see if the filtered file has content
if [ ! -s "${FILTERED_GFF3_FOR_PROCESSING}" ]; then
    AUGUSTUS_LINES_IN_ORIGINAL=$(grep -c $'\tAUGUSTUS\t' "${ORIGINAL_AUGUSTUS_GFF3_FILE}")
    if [ "${AUGUSTUS_LINES_IN_ORIGINAL}" -eq "0" ]; then
        echo "ERROR: No lines with 'AUGUSTUS' as source found in the original GFF3: ${ORIGINAL_AUGUSTUS_GFF3_FILE}"
        echo "The filtered GFF3 file is empty. Cannot proceed with sequence extraction."
    else
        echo "ERROR: Filtered GFF3 file (${FILTERED_GFF3_FOR_PROCESSING}) is empty or not created, but AUGUSTUS lines were present in the original."
        echo "Please check the awk command and the original GFF3 content near AUGUSTUS entries."
    fi
    exit 1
fi
echo "------------------------------------"


# --- Step 2: Extract Protein, CDS, and Transcript sequences using gffread ---
# This step uses the *FILTERED* GFF3 file created in Step 1.
echo ""
echo "------------------------------------"
echo "Step 2: Extracting sequences using gffread..."
echo "Using Filtered GFF3: ${FILTERED_GFF3_FOR_PROCESSING}"
echo "Using Genome FASTA: ${MASKED_GENOME_FASTA}"
echo "------------------------------------"

gffread "${FILTERED_GFF3_FOR_PROCESSING}" \
         -g "${MASKED_GENOME_FASTA}" \
         -x "${CDS_FNA_OUTPUT}" \
         -y "${PROTEIN_FAA_OUTPUT}" \
         -w "${TRANSCRIPT_FNA_OUTPUT}" \
         || { echo "ERROR: gffread failed. Please check the filtered GFF3 file (${FILTERED_GFF3_FOR_PROCESSING}), the genome FASTA file (${MASKED_GENOME_FASTA}), and ensure gffread is loaded correctly."; exit 1; }

echo "Sequence extraction completed successfully."
echo "  CDS sequences (DNA) are in: ${CDS_FNA_OUTPUT}"
echo "  Protein sequences (Amino Acid) are in: ${PROTEIN_FAA_OUTPUT}"
echo "  Transcript sequences (cDNA/mRNA) are in: ${TRANSCRIPT_FNA_OUTPUT}"
echo "------------------------------------"
echo "Script finished."