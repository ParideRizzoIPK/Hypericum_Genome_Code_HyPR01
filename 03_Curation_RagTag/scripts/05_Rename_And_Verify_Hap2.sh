#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --job-name=rename_verify_hap2
#SBATCH --output=logs/rename_%j.out
#SBATCH --error=logs/rename_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  

echo "#########################################################"
echo "Pipeline: Step 4 Continuation - Rename & Verify Haplotype 2"
echo "Date: $(date)"
echo "#########################################################"

# --- 1. ENVIRONMENT & MODULE INITIALIZATION ---
source ~/.bashrc
set -u            

# Load required cluster modules
module load minimap2/2.24
module load samtools/1.23.1
module load matplotlib/3.7.1

# Container Bind Fixes for Apptainer filesystem virtualization
export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# --- 2. DIRECTORIES & INPUT FILE PATHS ---
OUTPUT_DIR="/path/to/output/directory"
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

HAP1_SOURCE="/path/to/input/directory/Hap1_curated_chr_un.fasta"
HAP2_SOURCE="/path/to/input/directory/Hap2_curated_chr_un.fasta"

THREADS=${SLURM_CPUS_PER_TASK}

# ==============================================================================
# STEP 4.3: CENTRALIZATION & HAPLOTYPE 2 RE-ORDERING/RENAMING
# ==============================================================================
echo ">>>> Centralizing Files: Copying Haplotype 1 locally... <<<<"
cp "${HAP1_SOURCE}" ./Hap1_curated_chr_un.fasta

echo ">>>> Running Python Sequence Renaming and Sorting Engine... <<<<"
cat << 'EOF' > rename_and_sort_hap2.py
import sys

input_fa = "Hap2_curated_chr_un.fasta"
output_fa = "Hap2_curated_chr_un_RENAMED.fasta"

# Define the collinear mapping dictionary decoded from synteny report
mapping_matrix = {
    'Chr01': 'Chr01',
    'Chr05': 'Chr02',
    'Chr02': 'Chr03',
    'Chr03': 'Chr04',
    'Chr04': 'Chr05',
    'Chr06': 'Chr06',
    'Chr07': 'Chr07',
    'Chr08': 'Chr08',
    'UnChr': 'UnChr'
}

# Desired chronological structural order for final fasta file output
target_order = [f"Chr{i:02d}" for i in range(1, 9)] + ["UnChr"]

# Parse input fasta records into memory
records = {}
current_header = None
current_seq = []

with open(input_fa, 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if current_header:
                records[current_header] = "".join(current_seq)
            current_header = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    if current_header:
        records[current_header] = "".join(current_seq)

# Build the new sequence map dictionary using updated names
renamed_records = {}
for old_name, seq in records.items():
    clean_old_name = old_name.split()[0]
    if clean_old_name in mapping_matrix:
        new_name = mapping_matrix[clean_old_name]
        renamed_records[new_name] = seq
    else:
        # Fallback if unexpected minor fragments exist
        renamed_records[clean_old_name] = seq

# Write records out into final sequence destination matching target_order
with open(output_fa, 'w') as out:
    for name in target_order:
        if name in renamed_records:
            out.write(f">{name}\n")
            seq_str = renamed_records[name]
            # Maintain structural 80 character wrap boundaries
            for i in range(0, len(seq_str), 80):
                out.write(seq_str[i:i+80] + "\n")

print("Haplotype 2 sequence header renaming and sorting completed.")
EOF

# Link source Hap2 temporarily for script pipeline visibility
ln -sf "${HAP2_SOURCE}" ./Hap2_curated_chr_un.fasta
python3 rename_and_sort_hap2.py
rm -f rename_and_sort_hap2.py ./Hap2_curated_chr_un.fasta

# ==============================================================================
# STEP 4.4: RE-INDEX AND NEW PAF ALIGNMENT GENERATION
# ==============================================================================
echo ">>>> Re-indexing clean localized fasta sets... <<<<"
samtools faidx Hap1_curated_chr_un.fasta
samtools faidx Hap2_curated_chr_un_RENAMED.fasta

echo ">>>> Running Verification Minimap2 Alignment (Hap2 Renamed = Target)... <<<<"
minimap2 -x asm5 -t "${THREADS}" Hap2_curated_chr_un_RENAMED.fasta Hap1_curated_chr_un.fasta > Hap1_vs_Hap2_renamed_hap2.paf

# ==============================================================================
# STEP 4.5: RERUN PYTHON PLOTTING ENGINE (RENAMED LABEL OUTPUTS)
# ==============================================================================
echo ">>>> Compiling Verified Metrics & Plotting Chromosome Diagnostics... <<<<"

cat << 'EOF' > generate_renamed_outputs.py
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def get_merged_coverage(intervals):
    if not intervals:
        return 0
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]
    for current in intervals[1:]:
        prev_start, prev_end = merged[-1]
        curr_start, curr_end = current
        if curr_start <= prev_end:
            merged[-1] = (prev_start, max(prev_end, curr_end))
        else:
            merged.append(current)
    return sum(end - start + 1 for start, end in merged)

def parse_fai(fai_path):
    lengths = {}
    with open(fai_path, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.split('\t')
                lengths[parts[0]] = int(parts[1])
    return lengths

hap1_lens = parse_fai("Hap1_curated_chr_un.fasta.fai")
hap2_lens = parse_fai("Hap2_curated_chr_un_RENAMED.fasta.fai")

chr_order = [f"Chr{i:02d}" for i in range(1, 9)]
all_seqs_order = chr_order + ["UnChr"]

alignment_data = []
pair_aligned_bases = {}
hap1_intervals = {}
hap2_intervals = {}

with open("Hap1_vs_Hap2_renamed_hap2.paf", 'r') as f:
    for line in f:
        if not line.strip():
            continue
        p = line.split('\t')
        qname, qstart, qend, strand, tname, tstart, tend, block_len = p[0], int(p[2]), int(p[3]), p[4], p[5], int(p[7]), int(p[8]), int(p[10])
        
        if block_len >= 50000:
            alignment_data.append((qname, qstart, qend, strand, tname, tstart, tend))
            pair = (qname, tname)
            pair_aligned_bases[pair] = pair_aligned_bases.get(pair, 0) + block_len
            if pair not in hap1_intervals:
                hap1_intervals[pair] = []
                hap2_intervals[pair] = []
            hap1_intervals[pair].append((qstart, qend))
            hap2_intervals[pair].append((tstart, tend))

report_lines = []
for h1 in all_seqs_order:
    best_h2 = "None"
    max_bases = -1
    for h2 in all_seqs_order:
        bases = pair_aligned_bases.get((h1, h2), 0)
        if bases > max_bases:
            max_bases = bases
            best_h2 = h2
            
    if max_bases > 0:
        best_pair = (h1, best_h2)
        h1_cov_bp = get_merged_coverage(hap1_intervals[best_pair])
        h2_cov_bp = get_merged_coverage(hap2_intervals[best_pair])
        h1_pct = (h1_cov_bp / hap1_lens[h1]) * 100
        h2_pct = (h2_cov_bp / hap2_lens[best_h2]) * 100
    else:
        h1_pct, h2_pct = 0.0, 0.0
        
    report_lines.append(f"{h1:<15}{best_h2:<15}{h1_pct:<20.2f}{h2_pct:.2f}\n")

with open("curated_synteny_report_renamed_hap2.txt", 'w') as out:
    out.write("========================================================================\n")
    out.write("Verified Synteny Report: Renamed Haplotype 2 Check (Blocks >= 50kb)\n")
    out.write("========================================================================\n")
    out.write(f"{'Hap1 Chrom':<15}{'Best Hap2 Match':<15}{'Hap1 Coverage (%)':<20}{'Hap2 Coverage (%)'}\n")
    out.write("------------------------------------------------------------------------\n")
    out.writelines(report_lines)
    out.write("========================================================================\n")

def build_dotplot(seq_list, filename, title_label):
    def get_offsets(seq_list, lens):
        offsets = {}
        curr = 0
        for name in seq_list:
            offsets[name] = curr
            curr += lens.get(name, 0)
        return offsets, curr
        
    x_offsets, total_x = get_offsets(seq_list, hap2_lens)
    y_offsets, total_y = get_offsets(seq_list, hap1_lens)
    
    fig, ax = plt.subplots(figsize=(10, 10))
    for qname, qstart, qend, strand, tname, tstart, tend in alignment_data:
        if qname in seq_list and tname in seq_list:
            global_x1 = x_offsets[tname] + tstart
            global_x2 = x_offsets[tname] + tend
            global_y1 = y_offsets[qname] + qstart
            global_y2 = y_offsets[qname] + qend
            
            if strand == '-':
                ax.plot([global_x1, global_x2], [global_y2, global_y1], color='black', linewidth=1.2)
            else:
                ax.plot([global_x1, global_x2], [global_y1, global_y2], color='black', linewidth=1.2)
                
    for name in seq_list:
        if x_offsets[name] > 0:
            ax.axvline(x_offsets[name], color='gray', linestyle='--', linewidth=0.5)
        if y_offsets[name] > 0:
            ax.axhline(y_offsets[name], color='gray', linestyle='--', linewidth=0.5)
            
    x_ticks = [x_offsets[n] + hap2_lens.get(n, 0)/2 for n in seq_list]
    y_ticks = [y_offsets[n] + hap1_lens.get(n, 0)/2 for n in seq_list]
    
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(seq_list, rotation=45)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(seq_list)
    
    ax.set_xlim(0, total_x)
    ax.set_ylim(0, total_y)
    ax.set_xlabel("Haplotype 2 Assemblies (Target Coordinate Space)")
    ax.set_ylabel("Haplotype 1 Assemblies (Query Coordinate Space)")
    ax.set_title(title_label, fontsize=14, fontweight='bold', pad=15)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

build_dotplot(chr_order, "Hap1_vs_Hap2_Chromosomes_Only_renamed_hap2.png", "Synteny Map Verification: Renamed Clean Chromosomes")
build_dotplot(all_seqs_order, "Hap1_vs_Hap2_Whole_Genome_renamed_hap2.png", "Synteny Map Verification: Renamed Whole Genome")
print("Verification processing successful.")
EOF

python3 generate_renamed_outputs.py
rm -f generate_renamed_outputs.py

echo "#########################################################"
echo "PIPELINE COMPLETED."
echo "Centralized Outputs written to: ${OUTPUT_DIR}"
echo "  - Unified Hap1 Copy:  Hap1_curated_chr_un.fasta"
echo "  - Sorted Renamed Hap2: Hap2_curated_chr_un_RENAMED.fasta"
echo "  - Verification Summary: curated_synteny_report_renamed_hap2.txt"
echo "  - Check Figures:       *_renamed_hap2.png"
echo "#########################################################"