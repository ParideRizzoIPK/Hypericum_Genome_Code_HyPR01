#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00       # Post-processing text streams is rapid; 4 hours is highly safe
#SBATCH --job-name=define_chromosomes
#SBATCH --output=logs/defchr_%j.out
#SBATCH --error=logs/defchr_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  

echo "#########################################################"
echo "Pipeline: Step 3 - Define Chromosomes & Construct UnChr"
echo "Date: $(date)"
echo "#########################################################"

# --- 1. ENVIRONMENT INITIALIZATION ---
source ~/.bashrc
set -u            

# --- 2. DIRECTORIES & INPUT FILE PATHS ---
OUTPUT_DIR="/path/to/output/directory"
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

# Explicit inputs matching our normalized footprint adjustments
HAP1_IN="/path/to/input/directory/HyPR-01_Hap1_ragtag_scaffold.fasta"
HAP2_IN="/path/to/input/directory/Hap2_YAHS_min100kb.fa"

# ==============================================================================
# STEP 3.1: EMBEDDED PYTHON PROCESSING ENGINE
# ==============================================================================
# This script handles multi-line fasta records, ranks by length, translates 
# identifiers, inserts the 10kb N spacer padding, and records 1-based coordinates.
cat << 'EOF' > process_haplotype.py
import sys

input_fasta = sys.argv[1]
output_prefix = sys.argv[2]

print(f"Processing: {input_fasta} with prefix {output_prefix}...")

# 1. Parse FASTA records into memory
records = []
current_header = None
current_seq = []

with open(input_fasta, 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if current_header:
                records.append((current_header, "".join(current_seq)))
            current_header = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    if current_header:
        records.append((current_header, "".join(current_seq)))

# 2. Sort all scaffolds by length descending
records.sort(key=lambda x: len(x[1]), reverse=True)

# 3. Initialize Output Handles
map_file = f"{output_prefix}_chr_map.tsv"
fasta_file = f"{output_prefix}_curated_chr_un.fasta"
index_file = f"{output_prefix}_UnChr_index.tsv"

with open(map_file, 'w') as mf, open(fasta_file, 'w') as ff, open(index_file, 'w') as idx_f:
    # Write Index Header using standard structural naming conventions
    idx_f.write("original_scaffold\tUnChr_start\tUnChr_end\n")
    
    # Process top 8 records as clean chromosomes
    for i in range(min(8, len(records))):
        orig_header = records[i][0]
        clean_id = orig_header.split()[0]
        chr_name = f"Chr{i+1:02d}"
        seq = records[i][1]
        
        # Record tracking map and stream fasta blocks
        mf.write(f"{clean_id}\t{chr_name}\n")
        ff.write(f">{chr_name}\n")
        for j in range(0, len(seq), 80):
            ff.write(seq[j:j+80] + "\n")
            
    # Process all remaining records as concatenated unplaced fragments
    unplaced_records = records[8:]
    if unplaced_records:
        unchr_seq_parts = []
        current_pos = 1
        gap_string = "N" * 10000  # Formally hardcoded 10kb safe gene-prediction buffer
        
        for k, (orig_header, seq) in enumerate(unplaced_records):
            clean_id = orig_header.split()[0]
            mf.write(f"{clean_id}\tUNPLACED\n")
            
            # Calculate standard 1-based inclusive coordinates
            start_pos = current_pos
            end_pos = current_pos + len(seq) - 1
            idx_f.write(f"{clean_id}\t{start_pos}\t{end_pos}\n")
            
            unchr_seq_parts.append(seq)
            if k < len(unplaced_records) - 1:
                unchr_seq_parts.append(gap_string)
                current_pos = end_pos + 1 + 10000
            else:
                current_pos = end_pos + 1
                
        # Write unified UnChr molecule to FASTA
        unchr_total_seq = "".join(unchr_seq_parts)
        ff.write(">UnChr\n")
        for j in range(0, len(unchr_total_seq), 80):
            ff.write(unchr_total_seq[j:j+80] + "\n")

print(f"Successfully generated outputs for {output_prefix}.")
EOF

# ==============================================================================
# STEP 3.2: RUN THE TRANSFORMATION FOR BOTH HAPLOTYPES
# ==============================================================================
echo ">>>> Starting Post-Processing for Haplotype 1... <<<<"
python3 process_haplotype.py "${HAP1_IN}" "Hap1"

echo ">>>> Starting Post-Processing for Haplotype 2... <<<<"
python3 process_haplotype.py "${HAP2_IN}" "Hap2"

# Clean up the execution script from the output directory
rm -f process_haplotype.py

echo "#########################################################"
echo "PIPELINE STEP COMPLETE."
echo "Outputs successfully generated in: ${OUTPUT_DIR}"
echo "  - Hap1 Fasta: Hap1_curated_chr_un.fasta"
echo "  - Hap2 Fasta: Hap2_curated_chr_un.fasta"
echo "#########################################################"