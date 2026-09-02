import os
import matplotlib.pyplot as plt
import numpy as np

# Determine current script directory to ensure outputs land in the same folder
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Data matrices from Merqury outputs
chromosomes = [f"Chr{i:02d}" for i in range(1, 9)] + ["UnChr"]
hap1_qv = [65.081, 67.5256, 65.9478, 65.4518, 62.998, 67.1788, 65.9226, 68.984, 49.2526]
hap2_qv = [66.3876, 68.6721, 70.6444, 70.6299, 65.6121, 64.0452, 67.9882, 62.197, 47.2028]

x = np.arange(len(chromosomes))
width = 0.35

fig, ax = plt.subplots(figsize=(12, 7), facecolor='white')
ax.set_facecolor('white')

# Plot bars with exact color scheme
bar1 = ax.bar(x - width/2, hap1_qv, width, label='Haplotype 1', color='#1f77b4', edgecolor='black', linewidth=0.8)
bar2 = ax.bar(x + width/2, hap2_qv, width, label='Haplotype 2', color='#ff7f0e', edgecolor='black', linewidth=0.8)

# Aesthetics matching publication specs
ax.set_ylabel('Consensus Quality Score (QV)', fontsize=12, fontweight='bold', labelpad=10)
ax.set_xlabel('Assembled Chrosomosome ID', fontsize=12, fontweight='bold', labelpad=10)
ax.set_title('Merqury Chromosome-Level Base Accuracy Resolution (HyPR01)', fontsize=14, fontweight='bold', pad=20)
ax.set_xticks(x)
ax.set_xticklabels(chromosomes, fontsize=11)
ax.set_ylim(0, 80)

# Reference line for QV 60
ax.axhline(y=60, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)

# FIXED: 'va=bottom' sits directly ON TOP of the y=60 dotted line
ax.text(7.7, 60.5, 'QV 60 Threshold', color='#555555', fontsize=9.5,
        style='italic', ha='right', va='bottom',
        bbox=dict(boxstyle='square,pad=0.15', facecolor='white', edgecolor='none', alpha=0.85))

ax.legend(loc='upper right', fontsize=11, frameon=True, shadow=False)
ax.grid(axis='y', linestyle=':', alpha=0.5)

# Value annotations over bar columns
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

# Save PNG and PDF locally in the same directory as the script
out_png = os.path.join(SCRIPT_DIR, 'HyPR01_Chromosome_QV_Comparison.png')
out_pdf = os.path.join(SCRIPT_DIR, 'HyPR01_Chromosome_QV_Comparison.pdf')

plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(out_pdf, dpi=300, bbox_inches='tight', facecolor='white')

print(f"Saved PNG: {out_png}")
print(f"Saved PDF: {out_pdf}")