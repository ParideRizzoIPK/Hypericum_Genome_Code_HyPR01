#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=12:00:00       
#SBATCH --job-name=synteny_assessment
#SBATCH --output=logs/synteny_%j.out
#SBATCH --error=logs/synteny_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  

echo "#########################################################"
echo "Pipeline: Step 4 - Whole-Genome Synteny and Dotplot (Path Fixed)"
echo "Date: $(date)"
echo "#########################################################"

# --- 1. ENVIRONMENT & MODULE INITIALIZATION ---
source ~/.bashrc
set -u            

# Load required cluster modules
module load minimap2/2.24
module load samtools/1.23.1
module load matplotlib/3.7.1

# Container Bind Fixes to ensure Apptainer modules can access filesystem paths
export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# --- 2. DIRECTORIES & INPUT FILE PATHS ---
OUTPUT_DIR="/path/to/output/directory"
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

HAP1_FA="/path/to/input/directory/Hap1_curated_chr_un.fasta"
HAP2_FA="/path/to/input/directory/Hap2_curated_chr_un.fasta"

THREADS=${SLURM_CPUS_PER_TASK}

# ==============================================================================
# STEP 4.1: WHOLE-GENOME ALIGNMENT & INDEXING
# ==============================================================================
echo ">>>> Generating Assembly Indexes via Samtools... <<<<"
samtools faidx "${HAP1_FA}"
samtools faidx "${HAP2_FA}"

# FIX: Explicitly copy the generated indexes locally so Python can read them
cp "${HAP1_FA}.fai" ./Hap1_curated_chr_un.fasta.fai
cp "${HAP2_FA}.fai" ./Hap2_curated_chr_un.fasta.fai

echo ">>>> Running Minimap2 Alignment (Hap2 = Target/X, Hap1 = Query/Y)... <<<<"
minimap2 -x asm5 -t "${THREADS}" "${HAP2_FA}" "${HAP1_FA}" > Hap1_vs_Hap2.paf

# ==============================================================================
# STEP 4.2: EMBEDDED PYTHON METRICS COMPILER & PLOTTING ENGINE
# ==============================================================================
echo ">>>> Synthesizing Metrics and Rendering Dotplots via Matplotlib... <<<<"

cat << 'EOF' > generate_synteny_outputs.py
import os
import matplotlib
matplotlib.use('Agg')  # Force headless rendering on cluster environment
import matplotlib.pyplot as plt

# 1. Helper function for Coordinate-Collapsing Array Technique (Interval Merging)
def get_merged_coverage(intervals):
    if not intervals:
        return 0
    # Sort by start position
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

# 2. Parse Fasta Index (.fai) files to obtain true sequence lengths
def parse_fai(fai_path):
    lengths = {}
    with open(fai_path, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.split('\t')
                lengths[parts[0]] = int(parts[1])
    return lengths

hap1_lens = parse_fai("Hap1_curated_chr_un.fasta.fai")
hap2_lens = parse_fai("Hap2_curated_chr_un.fasta.fai")

# Define target tracking arrays
chr_order = [f"Chr{i:02d}" for i in range(1, 9)]
all_seqs_order = chr_order + ["UnChr"]

# Initialize structures to capture data pairs
alignment_data = []
pair_aligned_bases = {}  # (h1, h2) -> sum of raw alignment block lengths
hap1_intervals = {}      # (h1, h2) -> list of (start, end)
hap2_intervals = {}      # (h1, h2) -> list of (start, end)

# 3. Parse PAF and filter elements >= 50 kb
with open("Hap1_vs_Hap2.paf", 'r') as f:
    for line in f:
        if not line.strip():
            continue
        p = line.split('\t')
        
        qname = p[0]
        qstart, qend = int(p[2]), int(p[3])
        strand = p[4]
        tname = p[5]
        tstart, tend = int(p[7]), int(p[8])
        block_len = int(p[10])
        
        # Filter alignment noise below 50,000 bp
        if block_len >= 50000:
            alignment_data.append((qname, qstart, qend, strand, tname, tstart, tend))
            
            pair = (qname, tname)
            pair_aligned_bases[pair] = pair_aligned_bases.get(pair, 0) + block_len
            
            if pair not in hap1_intervals:
                hap1_intervals[pair] = []
                hap2_intervals[pair] = []
            
            hap1_intervals[pair].append((qstart, qend))
            hap2_intervals[pair].append((tstart, tend))

# 4. Determine best match chromosome pairs and calculate metrics
report_lines = []
for h1 in all_seqs_order:
    best_h2 = "None"
    max_bases = -1
    
    # Identify best matching target molecule based on total block length
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

# Write out the structural curation metrics table
with open("curated_synteny_report.txt", 'w') as out:
    out.write("========================================================================\n")
    out.write("Whole-Genome Synteny & Structural Difference Report (Blocks >= 50kb)\n")
    out.write("========================================================================\n")
    out.write(f"{'Hap1 Chrom':<15}{'Best Hap2 Match':<15}{'Hap1 Coverage (%)':<20}{'Hap2 Coverage (%)'}\n")
    out.write("------------------------------------------------------------------------\n")
    out.writelines(report_lines)
    out.write("========================================================================\n")

# 5. Visual Rendering Module: Plotting Function
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
    
    # Plot filtered line segments onto canvas space
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
                
    # Draw chromosome boundary grid dividers
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

# Render both structural visual modes
build_dotplot(chr_order, "Hap1_vs_Hap2_Chromosomes_Only.png", "Synteny Map: Clean Chromosomes (Chr01 - Chr08)")
build_dotplot(all_seqs_order, "Hap1_vs_Hap2_Whole_Genome.png", "Synteny Map: Whole Genome (Chromosomes + UnChr)")
print("Outputs processed successfully.")
EOF

python3 generate_synteny_outputs.py

# Clean up local indices and script copies
rm -f generate_synteny_outputs.py Hap1_curated_chr_un.fasta.fai Hap2_curated_chr_un.fasta.fai

echo "#########################################################"
echo "PIPELINE COMPLETED SUCCESSFULLY."
echo "Generated Outputs in: ${OUTPUT_DIR}"
echo "#########################################################"