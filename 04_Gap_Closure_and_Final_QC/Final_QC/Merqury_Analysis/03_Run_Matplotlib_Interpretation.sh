#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=00:30:00
#SBATCH --job-name=merqury_plot
#SBATCH --output=plot_%j.out
#SBATCH --error=plot_%j.err

# Load the required python plotting environment
module load matplotlib/3.7.1

# --- Write out the plotting script, then run it (merged from the former
# standalone Matplotlib_Merqury_Interpretation.py) ---
SCRIPT_PATH="Matplotlib_Merqury_Interpretation.py"

cat > "${SCRIPT_PATH}" << 'PYEOF'
import matplotlib.pyplot as plt
import numpy as np

def load_qv(path):
    """Parse a Merqury .qv file (tab-separated: seqname, kmers_only_in_asm,
    kmers_total, QV, error_rate) into parallel (names, qv_values) lists,
    preserving file order."""
    names, qvs = [], []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            names.append(fields[0])
            qvs.append(float(fields[3]))
    return names, qvs

HAP1_QV_FILE = "/path/to/input/directory/HyPR01_Final_QC.Hap1_curated_chr_un_gapclosed.qv"
HAP2_QV_FILE = "/path/to/input/directory/HyPR01_Final_QC.Hap2_curated_chr_un_gapclosed.qv"

chromosomes, hap1_qv = load_qv(HAP1_QV_FILE)
hap2_chromosomes, hap2_qv = load_qv(HAP2_QV_FILE)
assert chromosomes == hap2_chromosomes, (
    f"Hap1/Hap2 .qv files list different sequences or order: "
    f"{chromosomes} vs {hap2_chromosomes}"
)

x = np.arange(len(chromosomes))
width = 0.35

fig, ax = plt.subplots(figsize=(12, 7))

# Plot bars with classic scientific palettes
bar1 = ax.bar(x - width/2, hap1_qv, width, label='Haplotype 1', color='#1f77b4', edgecolor='black', linewidth=0.8)
bar2 = ax.bar(x + width/2, hap2_qv, width, label='Haplotype 2', color='#ff7f0e', edgecolor='black', linewidth=0.8)

# Add chart aesthetics matching plant biology journal specs
ax.set_ylabel('Consensus Quality Score (QV)', fontsize=12, fontweight='bold', labelpad=10)
ax.set_xlabel('Assembled Genomic Molecule Identifier', fontsize=12, fontweight='bold', labelpad=10)
ax.set_title('Merqury Chromosome-Level Base Accuracy Resolution (HyPR01)', fontsize=14, fontweight='bold', pad=20)
ax.set_xticks(x)
ax.set_xticklabels(chromosomes, fontsize=11)
ax.set_ylim(0, 80)

# Reference indicators for standard accuracy boundaries
ax.axhline(y=60, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
ax.text(7.7, 61, 'QV 60 Threshold', color='gray', fontsize=9, style='italic', ha='right')

ax.legend(loc='upper right', fontsize=11, frameon=True, shadow=False)
ax.grid(axis='y', linestyle=':', alpha=0.5)

# Dynamic value labels over each bar column
def autolabel(bars):
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.1f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=9)

autolabel(bar1)
autolabel(bar2)

plt.tight_layout()
plt.savefig('HyPR01_Chromosome_QV_Comparison.png', dpi=300)
print("Publication figure successfully rendered as 'HyPR01_Chromosome_QV_Comparison.png'.")
PYEOF

python3 "${SCRIPT_PATH}"
