#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G            
#SBATCH --time=36:00:00       # Adjusted to 36 hours to bypass cluster maintenance locks
#SBATCH --job-name=final_qc_pipeline
#SBATCH --output=logs/qc_%j.out
#SBATCH --error=logs/qc_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  

echo "#########################################################"
echo "Pipeline: Step 6 - Final Whole-Genome Quality Control (BUSCO Fixed)"
echo "Date: $(date)"
echo "#########################################################"

# --- 1. ENVIRONMENT & MODULE INITIALIZATION ---
source ~/.bashrc
set -u            

# Load verified native cluster modules
module load quast/5.2.0
module load BUSCO/5.8.2
module load merqury/1.3
module load samtools/1.23.1
module load minimap2/2.24
module load matplotlib/3.7.1

# Container Bind Safety variables
export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# --- 2. CONFIGURATION AND DIRECTORIES ---
OUTPUT_DIR="/path/to/output/directory"
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

# FIX: Establish a centralized, shared download storage bucket for BUSCO lineages
BUSCO_DOWNLOAD_DIR="${OUTPUT_DIR}/busco_downloads"
mkdir -p "${BUSCO_DOWNLOAD_DIR}"

# Input Curated & Gapclosed Assemblies from Step 5
HAP1_GAPCLOSED="/path/to/input/directory/Hap1_curated_chr_un_gapclosed.fasta"
HAP2_GAPCLOSED="/path/to/input/directory/Hap2_curated_chr_un_gapclosed.fasta"
HIFI_READS="/path/to/input/directory/hypr01_hifi_reads.fastq"

THREADS=${SLURM_CPUS_PER_TASK}

# ==============================================================================
# MODULE 6.1: QUAST ANALYSIS (Chromosomes-Only Isolation)
# ==============================================================================
echo ">>>> INITIALIZING QUAST ANALYSIS... <<<<"
QUAST_WORKSPACE="${OUTPUT_DIR}/QUAST_Analysis"
mkdir -p "${QUAST_WORKSPACE}"

# Extract chromosomes only to isolate the N-count calculations from UnChr padding walls
awk '/^>UnChr/{p=0;next} /^>/{p=1} p' "${HAP1_GAPCLOSED}" > "${QUAST_WORKSPACE}/Hap1_chroms_only.fasta"
awk '/^>UnChr/{p=0;next} /^>/{p=1} p' "${HAP2_GAPCLOSED}" > "${QUAST_WORKSPACE}/Hap2_chroms_only.fasta"

cd "${QUAST_WORKSPACE}"
quast.py \
    -t "${THREADS}" \
    -o final_quast_report \
    Hap1_chroms_only.fasta \
    Hap2_chroms_only.fasta

rm -f Hap1_chroms_only.fasta Hap2_chroms_only.fasta


# ==============================================================================
# MODULE 6.2: SEQUENTIAL BUSCO VALIDATION (SHARED PATH STORAGE FIX)
# ==============================================================================
echo ">>>> INITIALIZING BUSCO COMPLETENESS AUDIT... <<<<"

# --- Haplotype 1 Sequential Run ---
BUSCO_HAP1_DIR="${OUTPUT_DIR}/BUSCO_Hap1"
mkdir -p "${BUSCO_HAP1_DIR}" && cd "${BUSCO_HAP1_DIR}"

echo "  -> Running BUSCO on Haplotype 1 [Lineage: embryophyta]..."
# FIX: Removed '--offline' and added '--download_path' pointing to our central share folder
busco -m genome -i "${HAP1_GAPCLOSED}" -o Hap1_embryophyta -l embryophyta_odb10 -c "${THREADS}" --download_path "${BUSCO_DOWNLOAD_DIR}"

echo "  -> Running BUSCO on Haplotype 1 [Lineage: eudicots]..."
busco -m genome -i "${HAP1_GAPCLOSED}" -o Hap1_eudicots -l eudicots_odb10 -c "${THREADS}" --download_path "${BUSCO_DOWNLOAD_DIR}"

# --- Haplotype 2 Sequential Run ---
BUSCO_HAP2_DIR="${OUTPUT_DIR}/BUSCO_Hap2"
mkdir -p "${BUSCO_HAP2_DIR}" && cd "${BUSCO_HAP2_DIR}"

echo "  -> Running BUSCO on Haplotype 2 [Lineage: embryophyta]..."
busco -m genome -i "${HAP2_GAPCLOSED}" -o Hap2_embryophyta -l embryophyta_odb10 -c "${THREADS}" --download_path "${BUSCO_DOWNLOAD_DIR}"

echo "  -> Running BUSCO on Haplotype 2 [Lineage: eudicots]..."
busco -m genome -i "${HAP2_GAPCLOSED}" -o Hap2_eudicots -l eudicots_odb10 -c "${THREADS}" --download_path "${BUSCO_DOWNLOAD_DIR}"


# ==============================================================================
# MODULE 6.3: MERQURY QV & COMPLETENESS PROFILE
# ==============================================================================
echo ">>>> INITIALIZING MERQURY ASSESSMENT... <<<<"
MERQURY_WORKSPACE="${OUTPUT_DIR}/Merqury_Analysis"
mkdir -p "${MERQURY_WORKSPACE}"
cd "${MERQURY_WORKSPACE}"

echo "  -> Generating .meryl database from raw consensus HiFi dataset..."
meryl count k=21 memory=128G threads="${THREADS}" "${HIFI_READS}" output hypr01_hifi.meryl

echo "  -> Evaluating Consensus Quality Metrics via Merqury..."
if command -v merqury.sh &> /dev/null; then
    MERQURY_CMD="merqury.sh"
else
    MERQURY_CMD="merqury"
fi

${MERQURY_CMD} hypr01_hifi.meryl "${HAP1_GAPCLOSED}" "${HAP2_GAPCLOSED}" HyPR01_Final_QC


# ==============================================================================
# MODULE 6.4: RE-RUN SYNTENY STRUCTURAL AUDIT
# ==============================================================================
echo ">>>> INITIALIZING POST-GAP-CLOSING STRUCTURAL SYNTENY AUDIT... <<<<"
# (The inline Python code block below remains unchanged as it ran perfectly)
SYNTENY_WORKSPACE="${OUTPUT_DIR}/Synteny_PostGapClose"
mkdir -p "${SYNTENY_WORKSPACE}"
cd "${SYNTENY_WORKSPACE}"

samtools faidx "${HAP1_GAPCLOSED}"
samtools faidx "${HAP2_GAPCLOSED}"
cp "${HAP1_GAPCLOSED}.fai" ./Hap1_curated_chr_un_gapclosed.fasta.fai
cp "${HAP2_GAPCLOSED}.fai" ./Hap2_curated_chr_un_gapclosed.fasta.fai

minimap2 -x asm5 -t "${THREADS}" "${HAP2_GAPCLOSED}" "${HAP1_GAPCLOSED}" > PostGapClose_Hap1_vs_Hap2.paf

cat << 'EOF' > generate_post_gapclose_plots.py
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

hap1_lens = parse_fai("Hap1_curated_chr_un_gapclosed.fasta.fai")
hap2_lens = parse_fai("Hap2_curated_chr_un_gapclosed.fasta.fai")

chr_order = [f"Chr{i:02d}" for i in range(1, 9)]
all_seqs_order = chr_order + ["UnChr"]

alignment_data = []
pair_aligned_bases = {}
hap1_intervals = {}
hap2_intervals = {}

with open("PostGapClose_Hap1_vs_Hap2.paf", 'r') as f:
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

with open("curated_synteny_report_post_gapclose.txt", 'w') as out:
    out.write("========================================================================\n")
    out.write("Post-GapClosing Synteny Stability Audit Report (Blocks >= 50kb)\n")
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

build_dotplot(chr_order, "Hap1_vs_Hap2_Chromosomes_Only_post_gapclose.png", "Synteny Map Stability: Post-GapClosing Chromosomes")
build_dotplot(all_seqs_order, "Hap1_vs_Hap2_Whole_Genome_post_gapclose.png", "Synteny Map Stability: Post-GapClosing Whole Genome")
print("Post-gapclose plotting loops finished successfully.")
EOF

python3 generate_post_gapclose_plots.py
rm -f generate_post_gapclose_plots.py Hap1_curated_chr_un_gapclosed.fasta.fai Hap2_curated_chr_un_gapclosed.fasta.fai

echo "#########################################################"
echo "FINAL QC COMPILATION LOOP ENDED."
echo "Review individual subdirectories inside: ${OUTPUT_DIR}"
echo "#########################################################"