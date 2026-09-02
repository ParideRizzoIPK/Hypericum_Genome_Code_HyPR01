#!/bin/bash
#SBATCH --job-name=HyPR01_Pairwise_Master
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=45G
#SBATCH --time=24:00:00
#SBATCH --output=logs/pairwise_master_%j.out
#SBATCH --error=logs/pairwise_master_%j.err

# Fail fast on any error, unset variable, or pipe failure
set -euo pipefail

# ==============================================================================
# 1. PATH & SANDBOX DIRECTORY TREE SETUP
# ==============================================================================
OLD_BASE_DIR="/path/to/your/directory"
NEW_BASE_DIR="$OLD_BASE_DIR/HyPR01_Hap2_all_vs_all"

CLEAN_DIR="$NEW_BASE_DIR/cleaned_data"
SYN_DIR="$NEW_BASE_DIR/mcscanx"
CIRC_DIR="$NEW_BASE_DIR/circos"
SCRIPTS_DIR="$NEW_BASE_DIR/scripts"
LOG_DIR="$NEW_BASE_DIR/logs"

# Construct a completely isolated sandbox subdirectory tree
mkdir -p "$CLEAN_DIR" "$SYN_DIR" "$CIRC_DIR/links" "$CIRC_DIR/configs" "$SCRIPTS_DIR" "$LOG_DIR"

# --- INTEGRATED REAL-TIME LOGGER ---
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="$LOG_DIR/pairwise_pipeline_$TIMESTAMP.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "==================================================="
echo " Starting Pairwise 1-vs-1 Master Synteny Pipeline  "
echo " Sandbox Path: $NEW_BASE_DIR"
echo "==================================================="

# ENVIRONMENT INJECTION: Bind cluster binaries natively
export PATH="/path/to/your/directory/micromamba/envs/synteny_circos/bin:/path/to/your/directory/micromamba/envs/last_env/bin:$PATH"
export SINGULARITY_BIND="/path/to/your/directory"

# ==============================================================================
# 2. ASSET STAGING & MANIHOT CORRECTION
# ==============================================================================
echo "📋 Staging local reference tracks into sandbox workspace..."
cp "$OLD_BASE_DIR/cleaned_data/"*.bed "$CLEAN_DIR/"
cp "$OLD_BASE_DIR/cleaned_data/"*.pep "$CLEAN_DIR/"

# AUTOMATED MANIHOT ID HARMONIZATION FILTER: Strips any leftover .CDS strings to ensure a clean match
python3 -c '
import os
clean_dir = "'$CLEAN_DIR'"
for ext in ["bed", "pep"]:
    path = os.path.join(clean_dir, f"man.{ext}")
    if os.path.exists(path):
        with open(path, "r") as f: lines = f.readlines()
        with open(path, "w") as f:
            for line in lines:
                if ext == "bed":
                    p = line.strip().split("\t")
                    if len(p) >= 4: p[3] = p[3].split(".CDS")[0]
                    f.write("\t".join(p) + "\n")
                else:
                    if line.startswith(">"): line = line.strip().split()[0].split(".CDS")[0] + "\n"
                    f.write(line)
'

# ==============================================================================
# 3. HIGH-SPEED SEQUENTIAL SYNTENY GENERATION
# ==============================================================================
echo "==================================================="
echo " PHASE 3: RUNNING JCVI ORTHOLOG DETECTION ENGINE   "
echo "==================================================="
cd "$SYN_DIR"
ln -sf "$CLEAN_DIR"/*.pep .
ln -sf "$CLEAN_DIR"/*.bed .

SUBJECTS=("ara" "pop" "man" "hev" "ery")
for SBJ in "${SUBJECTS[@]}"; do
    echo "⚡ Processing alignment track sequentially: hyp vs $SBJ"
    /bin/rm -f hyp.${SBJ}.last hyp.${SBJ}.blast hyp.idx ${SBJ}.idx
    
    # FIXED: Reallocated full 16-core power sequentially to prevent I/O thrashing deadlocks
    python3 -m jcvi.compara.catalog ortholog hyp $SBJ --dbtype=prot --cpus=16 --cscore=.7 --no_strip_names > "$LOG_DIR/jcvi_${SBJ}.log" 2>&1
    python3 -m jcvi.compara.synteny screen hyp.${SBJ}.anchors hyp.${SBJ}.anchors.simple --minspan=5 --simple >> "$LOG_DIR/jcvi_${SBJ}.log" 2>&1
done

# ==============================================================================
# 4. CONFIGURATION GENERATOR FOR INDIVIDUAL SPECIES
# ==============================================================================
echo "==================================================="
echo " PHASE 4: GENERATING SECTOR ASSETS & CONFIGS       "
echo "==================================================="
cd "$NEW_BASE_DIR"

cat << 'EOF' > "$SCRIPTS_DIR/build_pairwise_configs.py"
import sys, os

clean_dir, circ_dir, syn_dir = sys.argv[1], sys.argv[2], sys.argv[3]
species_codes = ["ara", "pop", "man", "hev", "ery"]
species_names = ["Arabidopsis", "Populus", "Manihot", "Hevea", "Erythroxylum"]
prefix_map = {"ara": "At", "pop": "Pt", "man": "Me", "hev": "Hb", "ery": "En", "hyp": "Hp"}

# Distinct 8-color palette wheel for Hypericum chromosomes
color_wheel = {
    "HpChr01": "hyp_c1", "HpChr02": "hyp_c2", "HpChr03": "hyp_c3", "HpChr04": "hyp_c4",
    "HpChr05": "hyp_c5", "HpChr06": "hyp_c6", "HpChr07": "hyp_c7", "HpChr08": "hyp_c8"
}

for idx, target in enumerate(species_codes):
    target_name = species_names[idx]
    
    beds = {"hyp": {}, target: {}}
    lengths = {"hyp": {}, target: {}}
    min_size = 10000000 
    
    for sp in ["hyp", target]:
        path = os.path.join(clean_dir, f"{sp}.bed")
        if not os.path.exists(path): continue
        for line in open(path):
            p = line.split()
            if len(p) >= 3:
                beds[sp][p[3]] = (p[0], int(p[1]), int(p[2]))
                lengths[sp][p[0]] = max(lengths[sp].get(p[0], 0), int(p[2]))
                
    # 1. Output dedicated karyotype tracking file
    kary_path = os.path.join(circ_dir, "configs", f"karyotype_{target}.txt")
    with open(kary_path, "w") as out:
        for sp in ["hyp", target]:
            sp_prefix = prefix_map[sp]
            for chrom, length in sorted(lengths[sp].items()):
                if length >= min_size and chrom.startswith(f"{sp_prefix}Chr"):
                    c_val = color_wheel.get(chrom, "black") if sp == "hyp" else "tgt_color"
                    out.write(f"chr - {sp}_{chrom} {chrom} 0 {length} {c_val}\n")
                    
    # 2. Extract ribbon track link endpoints with matching transparency rules (_a4)
    links_path = os.path.join(circ_dir, "links", f"hyp_{target}.links")
    simple_file = os.path.join(syn_dir, f"hyp.{target}.anchors.simple")
    with open(links_path, "w") as out:
        if os.path.exists(simple_file):
            for line in open(simple_file):
                if line.startswith("#"): continue
                p = line.split()
                if p[0] in beds["hyp"] and p[2] in beds[target]:
                    chrA, sA1, eA1 = beds["hyp"][p[0]]; chrB, sB1, eB1 = beds[target][p[2]]
                    _, sA2, eA2 = beds["hyp"][p[1]]; _, sB2, eB2 = beds[target][p[3]]
                    
                    if chrA.startswith("HpChr") and chrB.startswith(prefix_map[target] + "Chr"):
                        base_color = color_wheel.get(chrA, "black")
                        out.write(f"hyp_{chrA} {min(sA1,sA2)} {max(eA1,eA2)} {target}_{chrB} {min(sB1,sB2)} {max(eB1,eB2)} color={base_color}_a4\n")

    # 3. Compute exact centered positions for titles
    labels_path = os.path.join(circ_dir, "configs", f"labels_{target}.txt")
    with open(labels_path, "w") as out:
        h_chroms = sorted([c for c in lengths["hyp"] if c.startswith("HpChr") and lengths["hyp"][c] >= min_size])
        if h_chroms:
            h_mid = h_chroms[len(h_chroms)//2]
            out.write(f"hyp_{h_mid} {lengths['hyp'][h_mid]//2} {lengths['hyp'][h_mid]//2} Hypericum_Hap2\n")
        t_chroms = sorted([c for c in lengths[target] if c.startswith(prefix_map[target] + "Chr") and lengths[target][c] >= min_size])
        if t_chroms:
            t_mid = t_chroms[len(t_chroms)//2]
            out.write(f"{target}_{t_mid} {lengths[target][t_mid]//2} {lengths[target][t_mid]//2} {target_name}\n")
EOF

python3 "$SCRIPTS_DIR/build_pairwise_configs.py" "$CLEAN_DIR" "$CIRC_DIR" "$SYN_DIR"

# ==============================================================================
# 5. CIRCOS INDEPENDENT VISUAL LAYOUT COMPILATION
# ==============================================================================
echo "==================================================="
echo " PHASE 5: CIRCOS GRAPHIC ENGINE COMPILATION         "
echo "==================================================="
cd "$CIRC_DIR"

for CODE in "ara" "pop" "man" "hev" "ery"; do
    echo "🎨 Building standalone config for outgroup: $CODE"
    
    cat << EOF > "configs/circos_${CODE}.conf"
karyotype = configs/karyotype_${CODE}.txt
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
tgt_color = 225,225,215
</colors>

<ideogram>
<spacing>
default = 0.015r
<pairwise hyp* ${CODE}*>
spacing = 0.15r
</pairwise>
<pairwise ${CODE}* hyp*>
spacing = 0.15r
</pairwise>
</spacing>
radius           = 0.68r
thickness        = 40p
fill             = yes
stroke_color     = dgrey
stroke_thickness = 2p
show_label       = yes
label_font       = bold
label_radius     = 1.08r
label_size       = 52p
label_parallel   = no
label_snuggle    = yes
</ideogram>

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
size           = 15p
show_label     = yes
label_size     = 18p
label_offset   = 4p
format         = %d
</tick>
</ticks>

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

<link pairwise_track>
file = links/hyp_${CODE}.links
</link>
</links>

<plots>
<plot>
type             = text
color            = black
file             = configs/labels_${CODE}.txt
r0               = 1.16r
r1               = 1.36r
show_links       = no
label_size       = 60p
label_font       = bold
</plot>
</plots>

<image>
<<include etc/image.conf>>
dir* = .
file* = Hypericum_vs_${CODE}_Final_large_labels
image_format* = png,svg
radius* = 1500p
background* = white
angle_offset* = -90
</image>

<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
EOF

    echo "Compiling vector files for: Hypericum vs $CODE"
    circos -conf "configs/circos_${CODE}.conf"
done

echo "==================================================="
echo " MASTER PAIRWISE PIPELINE EXECUTED SUCCESSFULLY   "
echo " All final images saved to directory: $CIRC_DIR"
echo "==================================================="