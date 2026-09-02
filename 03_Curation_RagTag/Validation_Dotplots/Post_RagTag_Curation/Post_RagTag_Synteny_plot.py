import sys
import re
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

# --- CONFIGURATION ---
paf_file = "/path/to/input/directory/Hap1_vs_Hap2_renamed_hap2.paf"
out_pdf = "/path/to/output/directory/synteny_dotplot_clean.pdf"
out_png = "/path/to/output/directory/synteny_dotplot_clean.png"
min_block_len = 50000  # 50kb threshold for clean synteny

# --- DATA STRUCTURES ---
hap1_lens = {}
hap2_lens = {}
alignments = []

print(f"Parsing PAF file: {paf_file}...")
try:
    with open(paf_file, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 12: 
                continue
            
            qname = cols[0]
            qlen = int(cols[1])
            tname = cols[5]
            tlen = int(cols[6])
            
            hap1_lens[qname] = qlen
            hap2_lens[tname] = tlen
            
            block_len = int(cols[10])
            if block_len >= min_block_len:
                alignments.append({
                    'q': qname, 't': tname, 
                    'qs': int(cols[2]), 'qe': int(cols[3]), 
                    'ts': int(cols[7]), 'te': int(cols[8]), 
                    'len': block_len
                })
except FileNotFoundError:
    print(f"Error: Could not find file at {paf_file}")
    sys.exit(1)

if not alignments:
    print(f"Warning: No alignments exceeding {min_block_len} bp were found.")
    sys.exit(0)


# --- CHROMOSOME SORTING LOGIC ---
def chr_sort_key(name):
    """
    Sorts Chromosomes numerically ('Chr01' -> 1) 
    and places unplaced scaffolds ('UnChr' / 'chrUn') at the end.
    """
    clean_name = name.replace("_RagTag", "")
    match = re.search(r'(\d+)', clean_name)
    if match and not ('un' in clean_name.lower()):
        return (0, int(match.group(1)))
    return (1, clean_name)

hap1_scaffolds = sorted(hap1_lens.keys(), key=chr_sort_key)
hap2_scaffolds = sorted(hap2_lens.keys(), key=chr_sort_key)

# Calculate cumulative offsets along X and Y axes
hap1_offsets = {}
curr_y = 0
for scaf in hap1_scaffolds:
    hap1_offsets[scaf] = curr_y
    curr_y += hap1_lens[scaf]
total_y = curr_y

hap2_offsets = {}
curr_x = 0
for scaf in hap2_scaffolds:
    hap2_offsets[scaf] = curr_x
    curr_x += hap2_lens[scaf]
total_x = curr_x


# --- PLOTTING ---
print("Generating publication-ready dotplot...")
plt.rcParams.update({'font.sans-serif': 'Arial', 'font.family': 'sans-serif'})
fig, ax = plt.subplots(figsize=(10, 10), facecolor='white')
ax.set_facecolor('white')

# 1. Grid Boundaries & Grid Lines
for scaf in hap2_scaffolds:
    x_end = (hap2_offsets[scaf] + hap2_lens[scaf])
    ax.axvline(x=x_end, color='#999999', linestyle='--', linewidth=0.6, zorder=1)

for scaf in hap1_scaffolds:
    y_end = (hap1_offsets[scaf] + hap1_lens[scaf])
    ax.axhline(y=y_end, color='#999999', linestyle='--', linewidth=0.6, zorder=1)

# 2. Alignment Lines (Solid Black for High Contrast Verification)
for aln in alignments:
    if aln['q'] in hap1_offsets and aln['t'] in hap2_offsets:
        x_start = hap2_offsets[aln['t']] + aln['ts']
        x_end = hap2_offsets[aln['t']] + aln['te']
        y_start = hap1_offsets[aln['q']] + aln['qs']
        y_end = hap1_offsets[aln['q']] + aln['qe']
        
        ax.plot([x_start, x_end], [y_start, y_end], color='black', linewidth=1.0, zorder=2)

# 3. Label Formatting & Axis Titles
h2_midpoints = [hap2_offsets[s] + (hap2_lens[s] / 2) for s in hap2_scaffolds]
hap2_clean_labels = [s.replace("_RagTag", "") for s in hap2_scaffolds]
ax.set_xticks(h2_midpoints)
ax.set_xticklabels(hap2_clean_labels, rotation=45, ha='right', fontsize=10)

h1_midpoints = [hap1_offsets[s] + (hap1_lens[s] / 2) for s in hap1_scaffolds]
hap1_clean_labels = [s.replace("_RagTag", "") for s in hap1_scaffolds]
ax.set_yticks(h1_midpoints)
ax.set_yticklabels(hap1_clean_labels, fontsize=10)

# Exact Matching Titles
ax.set_xlabel("Haplotype 2", fontsize=11, fontweight='normal', labelpad=12)
ax.set_ylabel("Haplotype 1", fontsize=11, fontweight='normal', labelpad=12)
ax.set_title("Synteny Map Verification: Post RagTag Curation Hap1 vs Hap2", fontsize=14, fontweight='bold', pad=15)

# Axis Box Spines & Limits
ax.set_xlim(0, total_x)
ax.set_ylim(0, total_y)

for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color('black')
    spine.set_linewidth(0.8)

plt.tight_layout()

# Save PDF Vector & 300 DPI PNG Raster Outputs
plt.savefig(out_pdf, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')

print(f"Saved vector PDF plot: {out_pdf}")
print(f"Saved 300 DPI PNG plot: {out_png}")