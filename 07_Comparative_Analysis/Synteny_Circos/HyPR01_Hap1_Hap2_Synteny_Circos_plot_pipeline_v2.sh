#!/bin/bash
#SBATCH --job-name=Hypericum_Hap1_Hap2_Synteny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --output=logs/hap_pipeline_%j.out
#SBATCH --error=logs/hap_pipeline_%j.err

# Fail fast on any error, unset variable, or pipe failure
set -euo pipefail

# ==============================================================================
# 1. PATH & STANDALONE DIRECTORY TREE SETUP
# ==============================================================================
PARENT_DIR="/path/to/your/directory"
BASE_DIR="$PARENT_DIR/HyPR01_Hap1_Hap2_circos_plot"

CLEAN_DIR="$BASE_DIR/cleaned_data"
SYN_DIR="$BASE_DIR/mcscanx"
CIRC_DIR="$BASE_DIR/circos"
SCRIPTS_DIR="$BASE_DIR/scripts"
LOG_DIR="$BASE_DIR/logs"

# Ensure isolated subfolders exist
mkdir -p "$CLEAN_DIR" "$SYN_DIR" "$CIRC_DIR/links" "$SCRIPTS_DIR" "$LOG_DIR"

# --- INTEGRATED REAL-TIME LOGGER ---
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="$LOG_DIR/hap1_vs_hap2_$TIMESTAMP.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "==================================================="
echo " Starting Isolated Hap1 vs Hap2 Synteny Pipeline   "
echo " Standardized Prefix: HyPR01_Hap1/2_ChrXX          "
echo " Color Mapping: Chromosome Color Wheel (_a4)       "
echo " Optimization: Mirrored Hemispheres & Balanced Gaps"
echo "==================================================="

# ENVIRONMENT INJECTION: Bind cluster binaries natively
export PATH="/path/to/your/directory/micromamba/envs/synteny_circos/bin:/path/to/your/directory/micromamba/envs/last_env/bin:$PATH"
export SINGULARITY_BIND="/path/to/your/directory"

# Source Data References
RAW_HAP1_FASTA="$PARENT_DIR/hap1.masked_unchr.fa"
RAW_HAP2_FASTA="$PARENT_DIR/hap2.masked_unchr.fa"
RAW_HAP1_GFF="$PARENT_DIR/harmonized_consensus_Hap1.gff3"
RAW_HAP2_GFF="$PARENT_DIR/harmonized_consensus_Hap2.gff3"

echo "==================================================="
echo " PHASE 1: GENOTYPE EXTRACTION & QUALITY FILTRATION "
echo "==================================================="

cat << 'EOF' > "$SCRIPTS_DIR/extract_haplotype_tracks.py"
import sys, os, re

def rev_comp(dna):
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', 'R':'Y', 'Y':'R', 'S':'S', 'W':'W', 'K':'M', 'M':'K', 'B':'V', 'D':'H', 'H':'D', 'V':'B'}
    return "".join(comp.get(base, base) for base in reversed(dna.upper()))

def translate(dna):
    mapping = {
        'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
        'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*', 'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
        'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
        'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
        'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
        'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
        'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
        'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'
    }
    dna = dna.upper()
    pep = []
    for i in range(0, len(dna) - 2, 3):
        pep.append(mapping.get(dna[i:i+3], 'X'))
    return "".join(pep)

hap_label, fasta_path, gff_path, out_pep, out_bed = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
prefix = "HyPR01_Hap1" if hap_label == "hap1" else "HyPR01_Hap2"

genome = {}
cur_chrom = None
cur_seq = []
with open(fasta_path, 'r') as f:
    for line in f:
        if line.startswith('>'):
            if cur_chrom: genome[cur_chrom] = "".join(cur_seq)
            cur_chrom = line.strip().split()[0][1:]
            cur_seq = []
        else:
            cur_seq.append(line.strip())
    if cur_chrom: genome[cur_chrom] = "".join(cur_seq)

transcripts = {}
bed_coords = {}
with open(gff_path, 'r') as f:
    for line in f:
        if line.startswith('#') or not line.strip(): continue
        parts = line.strip().split('\t')
        if len(parts) != 9: continue
        chrom, feat_type, start, end, strand = parts[0], parts[2], int(parts[3]), int(parts[4]), parts[6]
        attrs = {item.split('=', 1)[0]: item.split('=', 1)[1] for item in parts[8].split(';') if '=' in item}
        
        if feat_type in ['mRNA', 'transcript', 'primary_transcript']:
            tx_id = attrs.get('ID')
            if tx_id:
                tx_id_clean = tx_id.replace("transcript:", "").replace("rna:", "").replace("rna-", "").replace("cds-", "")
                bed_coords[tx_id_clean] = (chrom, start, end, strand)
        elif feat_type == 'CDS':
            parent = attrs.get('Parent')
            if parent:
                parent_clean = parent.replace("transcript:", "").replace("rna:", "").replace("rna-", "").replace("cds-", "")
                if parent_clean not in transcripts:
                    transcripts[parent_clean] = {'chrom': chrom, 'strand': strand, 'cds': []}
                transcripts[parent_clean]['cds'].append((start, end))

valid_txs = set()
with open(out_pep, 'w') as out:
    for tx_id, info in transcripts.items():
        chrom, strand = info['chrom'], info['strand']
        if chrom not in genome: continue
        sorted_cds = sorted(info['cds'], key=lambda x: x[0], reverse=(strand == '-'))
        full_cds_dna = "".join(genome[chrom][start-1:end] for start, end in sorted_cds)
        if strand == '-': full_cds_dna = rev_comp(full_cds_dna)
        
        protein_seq = translate(full_cds_dna)
        stop_count = protein_seq.count('*')
        
        if stop_count > 1: continue
        if stop_count == 1 and not protein_seq.endswith('*'): continue
        
        if protein_seq.endswith('*'): protein_seq = protein_seq[:-1]
        out.write(f">{tx_id}\n{protein_seq}\n")
        valid_txs.add(tx_id)

with open(out_bed, 'w') as bed_out:
    for tx_id in sorted(valid_txs):
        if tx_id in bed_coords:
            chrom, start, end, strand = bed_coords[tx_id]
            num_match = re.search(r'\d+', chrom)
            chrom_clean = f"{prefix}_Chr{int(num_match.group()):02d}" if num_match else f"{prefix}_{chrom}"
            bed_out.write(f"{chrom_clean}\t{start-1}\t{end}\t{tx_id}\t0\t{strand}\n")
EOF

echo "--> Formatting Haplotype 1 coding sequence tracks..."
python3 "$SCRIPTS_DIR/extract_haplotype_tracks.py" "hap1" "$RAW_HAP1_FASTA" "$RAW_HAP1_GFF" "$CLEAN_DIR/hap1.pep" "$CLEAN_DIR/hap1.bed"

echo "--> Formatting Haplotype 2 coding sequence tracks..."
python3 "$SCRIPTS_DIR/extract_haplotype_tracks.py" "hap2" "$RAW_HAP2_FASTA" "$RAW_HAP2_GFF" "$CLEAN_DIR/hap2.pep" "$CLEAN_DIR/hap2.bed"


echo "==================================================="
echo " PHASE 2: DIRECT PAIRWISE SYNTENY ALIGNMENTS       "
echo "==================================================="
cd "$SYN_DIR"
ln -sf "$CLEAN_DIR"/hap1.pep .
ln -sf "$CLEAN_DIR"/hap2.pep .
ln -sf "$CLEAN_DIR"/hap1.bed .
ln -sf "$CLEAN_DIR"/hap2.bed .

/bin/rm -f hap1.hap2.last hap1.hap2.blast hap1.idx hap2.idx

echo "⚡ Launching LAST reciprocal sequence search..."
python3 -m jcvi.compara.catalog ortholog hap1 hap2 --dbtype=prot --cpus=16 --cscore=.7 --no_strip_names > "$LOG_DIR/jcvi_alignment.log" 2>&1

echo "⚡ Resolving macro-collinear syntax anchors..."
python3 -m jcvi.compara.synteny screen hap1.hap2.anchors hap1.hap2.anchors.simple --minspan=5 --simple >> "$LOG_DIR/jcvi_alignment.log" 2>&1


echo "==================================================="
echo " PHASE 3: SCALING COORDINATES & PARSING RIBBONS    "
echo "==================================================="
cd "$BASE_DIR"

cat << 'EOF' > "$SCRIPTS_DIR/generate_circos_assets.py"
import sys, os, re
from collections import defaultdict

clean_dir = sys.argv[1]
circ_dir = sys.argv[2]
syn_dir = sys.argv[3]
fa_hap1 = sys.argv[4]
fa_hap2 = sys.argv[5]

open(os.path.join(circ_dir, "links", "hap1_hap2.links"), "w").close()

# Parse chromosome lengths directly from true sequence basepair distributions
chrom_lengths = {}
for hap_tag, fa_path, prefix in [("hap1", fa_hap1, "HyPR01_Hap1"), ("hap2", fa_hap2, "HyPR01_Hap2")]:
    cur_chrom, cur_len = None, 0
    for line in open(fa_path):
        if line.startswith(">"):
            if cur_chrom:
                n = re.search(r"\d+", cur_chrom)
                chrom_name = f"{prefix}_Chr{int(n.group()):02d}" if n else f"{prefix}_{cur_chrom}"
                chrom_lengths[(hap_tag, chrom_name)] = cur_len
            cur_chrom = line.strip().split()[0][1:]
            cur_len = 0
        else:
            cur_len += len(line.strip())
    if cur_chrom:
        n = re.search(r"\d+", cur_chrom)
        chrom_name = f"{prefix}_Chr{int(n.group()):02d}" if n else f"{prefix}_{cur_chrom}"
        chrom_lengths[(hap_tag, chrom_name)] = cur_len

# Custom color wheel mapping tracking individual chromosomes visually
color_wheel = {
    "Chr01": "hyp_c1", "Chr02": "hyp_c2", "Chr03": "hyp_c3", "Chr04": "hyp_c4",
    "Chr05": "hyp_c5", "Chr06": "hyp_c6", "Chr07": "hyp_c7", "Chr08": "hyp_c8"
}

# Output Standardized Karyotype File with Mirrored (Reversed H2) Genomes Order
with open(os.path.join(circ_dir, "karyotype.txt"), "w") as out:
    # 1. Plot Haplotype 1 in ascending order (Chr01 -> Chr08)
    for (sp, chrom), length in sorted(chrom_lengths.items()):
        if sp == "hap1" and "_Chr" in chrom:
            chrom_suffix = chrom.split("_")[-1]
            c_val = color_wheel.get(chrom_suffix, "grey")
            out.write(f"chr - {sp}_{chrom} {chrom_suffix} 0 {length} {c_val}\n")
            
    # 2. Plot Haplotype 2 in descending order (Chr08 -> Chr01) to face each other perfectly
    h2_chroms = [(chrom, length) for (sp, chrom), length in chrom_lengths.items() if sp == "hap2" and "_Chr" in chrom]
    for chrom, length in sorted(h2_chroms, key=lambda x: x[0], reverse=True):
        chrom_suffix = chrom.split("_")[-1]
        out.write(f"chr - hap2_{chrom} {chrom_suffix} 0 {length} tgt_color\n")

# Load BED files for feature lookups
def load_bed(pfx):
    d = {}
    path = os.path.join(clean_dir, f"{pfx}.bed")
    for line in open(path):
        p = line.split()
        if len(p) >= 4: d[p[3]] = (p[0], int(p[1]), int(p[2]))
    return d

beds = {"hap1": load_bed("hap1"), "hap2": load_bed("hap2")}

# Compile Ribbon Links files injecting direct translucent line colors (_a4)
simple_file = os.path.join(syn_dir, "hap1.hap2.anchors.simple")
if os.path.exists(simple_file):
    with open(os.path.join(circ_dir, "links", "hap1_hap2.links"), "w") as out:
        for line in open(simple_file):
            if line.startswith("#"): continue
            p = line.split()
            if p[0] in beds["hap1"] and p[2] in beds["hap2"]:
                chrA, sA1, eA1 = beds["hap1"][p[0]]
                chrB, sB1, eB1 = beds["hap2"][p[2]]
                _, sA2, eA2 = beds["hap1"][p[1]]
                _, sB2, eB2 = beds["hap2"][p[3]]
                
                chrom_suffix = chrA.split("_")[-1]
                base_color = color_wheel.get(chrom_suffix, "grey")
                out.write(f"hap1_{chrA} {min(sA1,sA2,eA1,eA2)} {max(sA1,sA2,eA1,eA2)} hap2_{chrB} {min(sB1,sB2,eB1,eB2)} {max(sB1,sB2,eB1,eB2)} color={base_color}_a4\n")

# Position sector titles exactly at center midpoints
sp_chroms = defaultdict(list)
for line in open(os.path.join(circ_dir, "karyotype.txt")):
    p = line.split()
    if len(p) >= 3:
        sp_chroms[p[2].split("_")[0]].append(p)

titles = {"hap1": "Hypericum_Haplotype_1", "hap2": "Hypericum_Haplotype_2"}
with open(os.path.join(circ_dir, "species_labels.txt"), "w") as out:
    for code, chrs in sp_chroms.items():
        if code in titles:
            mid = chrs[len(chrs)//2]
            mid_pos = int(mid[5]) // 2
            out.write(f"{mid[2]} {mid_pos} {mid_pos} {titles[code]}\n")
EOF

python3 "$SCRIPTS_DIR/generate_circos_assets.py" "$CLEAN_DIR" "$CIRC_DIR" "$SYN_DIR" "$RAW_HAP1_FASTA" "$RAW_HAP2_FASTA"


# ==============================================================================
# PHASE 4: CIRCOS MAP COMPILATION
# ==============================================================================
cd "$CIRC_DIR"

cat << 'EOF' > ideogram.conf
<ideogram>
<spacing>
# Intra-chromosome gap padding
default = 0.015r

# MIRRORED BOUNDARY RULES FIXED: Perfectly maps poles symmetrically (01-vs-01 and 08-vs-08)
<pairwise hap1_HyPR01_Hap1_Chr08 hap2_HyPR01_Hap2_Chr08>
spacing = 0.15r
</pairwise>
<pairwise hap2_HyPR01_Hap2_Chr01 hap1_HyPR01_Hap1_Chr01>
spacing = 0.15r
</pairwise>
</spacing>
radius           = 0.72r
thickness        = 40p
fill             = yes
stroke_color     = dgrey
stroke_thickness = 2p
show_label       = yes
label_font       = bold
label_radius     = 1.10r
label_size       = 38
label_parallel   = no
label_snuggle    = yes
</ideogram>
EOF

cat << 'EOF' > ticks.conf
show_ticks          = yes
show_tick_labels    = yes
<ticks>
radius           = 1.01r
color            = black
thickness        = 2p
multiplier       = 1e-6
format           = %d
<tick>
spacing        = 10u
size           = 10p
show_label     = no
</tick>
<tick>
spacing        = 20u
size           = 14p
show_label     = yes
label_size     = 18p
label_offset   = 4p
format         = %d
</tick>
</ticks>
EOF

cat << 'EOF' > circos.conf
karyotype = karyotype.txt
chromosomes_units = 1000000
chromosomes_display_default = yes

<colors>
hyp_c1 = 220,80,80
hyp_c2 = 230,140,70
hyp_c3 = 220,190,80
hyp_c4 = 90,180,100
hyp_c5 = 70,170,160
hyp_c6 = 80,140,220
hyp_c7 = 120,100,200
hyp_c8 = 200,90,160
tgt_color = 235,235,225
</colors>

<<include ideogram.conf>>
<<include ticks.conf>>

<links>
radius = 0.99r
crest  = 1
ribbon = yes
flat   = yes
bezier_radius = 0r
bezier_radius_purity = 0.5
<rules>
<rule>
condition = var(size1) < 100000 || var(size2) < 100000
show = no
</rule>
</rules>

<link hap1_hap2>
file = links/hap1_hap2.links
</link>
</links>

<plots>
<plot>
type             = text
color            = black
file             = species_labels.txt
r0               = 1.14r
r1               = 1.34r
show_links       = no
label_size       = 62p
label_font       = bold
</plot>
</plots>

<image>
<<include etc/image.conf>>
dir* = .
file* = Hypericum_Hap1_Hap2_Synteny
image_format* = png,svg
radius* = 1500p
background* = white
# OVERRIDE FIX: Applied standard explicit override syntax to clear core templates natively
angle_offset* = -90
</image>

<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
EOF

echo "[$(date)] Launching main vector map construction..."
circos -conf circos.conf

echo "==================================================="
echo " RUN SUCCESSFUL! GENOTYPE COMPARISON COMPLETE     "
echo " PNG Graphic: $CIRC_DIR/Hypericum_Hap1_Hap2_Synteny.png"
echo " SVG Graphic: $CIRC_DIR/Hypericum_Hap1_Hap2_Synteny.svg"
echo "==================================================="