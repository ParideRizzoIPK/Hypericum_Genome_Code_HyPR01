#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00       # This step is fast; 2 hours is more than enough
#SBATCH --job-name=ragtag_filter_100kb
#SBATCH --output=logs/filter_%j.out
#SBATCH --error=logs/filter_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  

echo "#########################################################"
echo "Pipeline: Step 2.1 - Filter RagTag Queries (>= 100 kb)"
echo "Date: $(date)"
echo "#########################################################"

# --- 1. ENVIRONMENT INITIALIZATION ---
source ~/.bashrc
set -u            

# --- 2. DIRECTORIES & FILES ---
OUTPUT_DIR="/path/to/output/directory"
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

HAP1_IN="/path/to/input/directory/HyPR-01_v10_telo.hap1_scaffolds_final.fa"
HAP2_IN="/path/to/input/directory/HyPR-01_v10_telo.hap2_scaffolds_final.fa"

HAP1_OUT="${OUTPUT_DIR}/Hap1_YAHS_min100kb.fa"
HAP2_OUT="${OUTPUT_DIR}/Hap2_YAHS_min100kb.fa"

# ==============================================================================
# STEP 2.1.1: CREATE THE INLINE PYTHON FILTERING ENGINE
# ==============================================================================
# This script reads FASTA records on the fly, counts their length, 
# and writes them out only if they meet the threshold criteria.
cat << 'EOF' > stream_filter.py
import sys

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]
min_size = int(sys.argv[3])

kept_count = 0
dropped_count = 0

with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
    current_header = ""
    current_seq = []
    
    def process_record(header, seq_list):
        global kept_count, dropped_count
        seq_str = "".join(seq_list)
        if len(seq_str) >= min_size:
            outfile.write(f"{header}\n")
            # Wrap FASTA sequence rows at a standard 80 characters
            for i in range(0, len(seq_str), 80):
                outfile.write(seq_str[i:i+80] + "\n")
            kept_count += 1
        else:
            dropped_count += 1

    for line in infile:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header:
                process_record(current_header, current_seq)
            current_header = line
            current_seq = []
        else:
            current_seq.append(line)
            
    # Process the final trailing record in the file
    if current_header:
        process_record(current_header, current_seq)

print(f"  Total Scaffolds Retained (>= {min_size:,} bp): {kept_count}")
print(f"  Total Scaffolds Filtered Out (< {min_size:,} bp): {dropped_count}")
EOF

# ==============================================================================
# STEP 2.1.2: EXECUTE FILTERING SEQUENTIALLY
# ==============================================================================
echo ">>>> Filtering Haplotype 1... <<<<"
python3 stream_filter.py "${HAP1_IN}" "${HAP1_OUT}" 100000

echo ">>>> Filtering Haplotype 2... <<<<"
python3 stream_filter.py "${HAP2_IN}" "${HAP2_OUT}" 100000

# Clean up the helper python script
rm -f stream_filter.py

echo "#########################################################"
echo "FILTERING COMPLETE."
echo "Filtered Outputs Saved in: ${OUTPUT_DIR}"
echo "#########################################################"