#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --partition=cpu        
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16     
#SBATCH --mem=32G              
#SBATCH --time=02:00:00        
#SBATCH --job-name=NUMT_NUPT_K19_control
#SBATCH --output=logs/numt_nupt_relaxed_k19_%j.out
#SBATCH --error=logs/numt_nupt_relaxed_k19_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  

# --- 1. ENVIRONMENT & MODULE INITIALIZATION ---
source ~/.bashrc
module load samtools/1.23.1
module load matplotlib/3.7.1

export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# ==============================================================================
#                      DIAGNOSTICS & LOGIC CHECK SECTION
# ==============================================================================
echo "[$(date)] Starting Pre-flight Diagnostic Checks..."

OUTPUT_DIR="/path/to/output/directory"
HAP1_GFF="${OUTPUT_DIR}/hap1_NUMT_NUPT.gff3"
HAP2_GFF="${OUTPUT_DIR}/hap2_NUMT_NUPT.gff3"

FASTA_DIR="/path/to/your/directory"
HAP1_FASTA="${FASTA_DIR}/hap1.softmasked.fa"
HAP2_FASTA="${FASTA_DIR}/hap2.softmasked.fa"

# Verify Essential Executables
for cmd in samtools python3; do
    if ! command -v $cmd &> /dev/null; then
        echo "ERROR: Required executable '$cmd' is missing from host environment." >&2
        exit 1
    fi
done

# Verify Python visualization toolkit is active
if ! python3 -c "import matplotlib" &> /dev/null; then
    echo "ERROR: 'matplotlib' library loading failed. Check module paths." >&2
    exit 1
fi
echo "✔ Toolchain and Python libraries verification passed."

# Verify Input Datasets
for file in "$HAP1_GFF" "$HAP2_GFF" "$HAP1_FASTA" "$HAP2_FASTA"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Target input asset missing: $file" >&2
        exit 1
    fi
done
echo "✔ Target input GFF3 and FASTA files verified."

if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi
cd "$OUTPUT_DIR" || exit 1
echo "[$(date)] Pre-flight diagnostics passed. Proceeding to execution payload."

# ==============================================================================
#                     DYNAMIC FASTA INDEX GENERATION
# ==============================================================================
echo "[$(date)] Generating assembly indices via containerized samtools faidx..."
samtools faidx "$HAP1_FASTA"
samtools faidx "$HAP2_FASTA"

# ==============================================================================
#               PYTHON3 TRIPLE-FORMAT HEADLESS VISUALIZER
# ==============================================================================
echo "[$(date)] Launching headless Python mapping engine..."

python3 - << 'EOF'
import os
import math
from collections import defaultdict
import matplotlib

# Force headless engine state to prevent X11 display server crashes
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

output_dir = "/path/to/output/directory"
hap1_gff_path = os.path.join(output_dir, "hap1_NUMT_NUPT.gff3")
hap2_gff_path = os.path.join(output_dir, "hap2_NUMT_NUPT.gff3")
hap1_fai_path = "/path/to/input/directory/hap1.softmasked.fa.fai"
hap2_fai_path = "/path/to/input/directory/hap2.softmasked.fa.fai"

chromosomes = [f"Chr0{i}" for i in range(1, 9)]

# Darker green optimization for NUPTs (#008F39) and rich orange for NUMTs
colors_brilliant = {"NUMT": "#FF6D00", "NUPT": "#008F39"}

def parse_individual_loci(gff_path, fai_path):
    chr_lens = {}
    with open(fai_path, 'r') as f:
        for line in f:
            if not line.strip(): continue
            parts = line.strip().split('\t')
            seqid = "Chr04" if parts[0] == "Ch4" else ("Chr05" if parts[0] == "Ch5" else parts[0])
            if seqid in chromosomes:
                chr_lens[seqid] = int(parts[1])
                
    loci_records = []
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            seqid = "Chr04" if parts[0] == "Ch4" else ("Chr05" if parts[0] == "Ch5" else parts[0])
            ftype = parts[2]
            if ftype not in ["NUMT", "NUPT"] or seqid == "UnChr": continue
            
            start, end = int(parts[3]), int(parts[4])
            length_kb = (end - start + 1) / 1000.0
            midpoint_mb = ((start + end) / 2.0) / 1000000.0
            
            loci_records.append({
                'seqid': seqid,
                'type': ftype,
                'midpoint_mb': midpoint_mb,
                'length_kb': length_kb
            })
                
    return chr_lens, loci_records

h1_lens, h1_loci = parse_individual_loci(hap1_gff_path, hap1_fai_path)
h2_lens, h2_loci = parse_individual_loci(hap2_gff_path, hap2_fai_path)

# Unified single column workspace canvas
fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(14, 8.5), sharex=False)
fig.patch.set_facecolor('white')

def render_structural_track(ax, chr_lens, loci_data, title_label):
    ax.set_facecolor('none')
    y_coords = {chromosomes[i]: 8 - i for i in range(8)}
    
    # Draw gray structural chromosome track baselines
    for chr_name, length in chr_lens.items():
        y = y_coords[chr_name]
        ax.hlines(y, xmin=0, xmax=length/1e6, colors='#E5E5E5', linewidth=11, zorder=1)
        
    # Plot true spatial coordinates scaled by insertion length
    for locus in loci_data:
        y = y_coords[locus['seqid']]
        marker_area = 20.0 + (locus['length_kb'] * 4.5) 
        
        ax.scatter(locus['midpoint_mb'], y, s=marker_area, c=colors_brilliant[locus['type']], 
                   alpha=0.85, edgecolors='none', zorder=2)

    # Clean layout limits and labeling architectures
    ax.set_ylim(0.2, 9.2)
    
    # Symmetrical maximum lock at exactly 60.0 Mb for horizontal coordinate tracks
    ax.set_xlim(-2, 60.0)
    
    ax.set_yticks(list(y_coords.values()))
    ax.set_yticklabels(list(y_coords.keys()), fontsize=11, fontweight='bold', color='#333333')
    
    ax.set_xlabel("Physical Position (Mb)", fontsize=12, fontweight='bold', labelpad=8, color='#222222')
    ax.xaxis.set_major_locator(plt.MultipleLocator(10))
    ax.xaxis.set_tick_params(labelsize=11, labelcolor='#333333', bottom=True)
        
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.grid(axis='x', linestyle=':', color='#E5E5E5', alpha=0.7)
    
    ax.text(-1.8, 8.7, title_label, fontsize=13, fontweight='bold', va='bottom', ha='left', color='#111111')

render_structural_track(ax0, h1_lens, h1_loci, "Haplotype 1 (Hap1)")
render_structural_track(ax1, h2_lens, h2_loci, "Haplotype 2 (Hap2)")

# ------------------------------------------------------------------------------
# PUBLICATION-GRADE LEGEND CONSTRUCTOR (DECONVOLUTED LOWER PLOT REGION)
# ------------------------------------------------------------------------------
# Color Legend Elements
color_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=colors_brilliant['NUPT'], markersize=10, label='NUPT (Chloroplast)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=colors_brilliant['NUMT'], markersize=10, label='NUMT (Mitochondrion)')
]

# Structural Size Benchmarks Legend Elements 
size_benchmarks = [2.0, 20.0, 100.0]
size_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#666666', alpha=0.7,
           markersize=math.sqrt(20.0 + b * 4.5), label=f"{int(b)} kb") for b in size_benchmarks
]

# Combined safe positioning lower canvas bounds to completely clear Chr01 lines
leg_color = ax1.legend(handles=color_elements, loc='upper left', bbox_to_anchor=(0.02, -0.28),
                       frameon=False, ncol=2, prop={'weight':'bold', 'size':11})
ax1.add_artist(leg_color)

ax1.legend(handles=size_elements, loc='upper right', bbox_to_anchor=(0.98, -0.28),
           frameon=False, ncol=3, title="Insertion Sequence Length", 
           prop={'weight':'bold', 'size':11}, title_fontproperties={'weight':'bold', 'size':11})

plt.subplots_adjust(bottom=0.18, hspace=0.35)

# Simultaneous triple format outputs
plt.savefig(os.path.join(output_dir, "HyPR01_organellar_density_map.png"), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(output_dir, "HyPR01_organellar_density_map.pdf"), format='pdf', bbox_inches='tight')
plt.savefig(os.path.join(output_dir, "HyPR01_organellar_density_map.svg"), format='svg', bbox_inches='tight')
plt.close()

print("✔ Redesigned triple-format visualization matrix exported successfully.")
EOF

echo "[$(date)] Pipeline workflow finished execution successfully."