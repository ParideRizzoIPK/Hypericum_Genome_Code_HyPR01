# SUPERSEDED / UNUSED -- kept for the audit trail only.
#
# An earlier, unified draft of the karyotype/links Circos config logic for
# the Malpighiales comparison (Figure 16). The kept pipeline script,
# ../Malpighiales_HyPR01_Hap2_Synteny_Circos_plot_pipeline_v1.sh, does not
# call this file -- it builds the same config via two separate inline
# `python3 -c` blocks with different arguments (3 positional args across
# two calls, vs. this file's single 4-arg entry point) and different code,
# not a byte-for-byte match the way this repo's other embedded/static-copy
# script pairs are. Do not treat this file as the source for Figure 16 --
# see PIPELINE_OVERVIEW.md.

import sys, os, re
from collections import defaultdict

clean_dir = sys.argv[1]
circ_dir = sys.argv[2]
syn_dir = sys.argv[3]
hyp_fasta = sys.argv[4]

species_order = ["hyp", "pop", "man", "hev", "ara", "ery"]
species_info = {"hyp": "102,0,153", "pop": "255,127,80", "man": "241,196,15", "hev": "46,204,113", "ara": "112,128,144", "ery": "0,191,255"}
PREFIX_MAP = {"ara":"At", "pop":"Pt", "man":"Me", "hev":"Hb", "ery":"En", "hyp":"Hp"}
species_names = {"hyp": "Hypericum_Hap2", "pop": "Populus", "man": "Manihot", "hev": "Hevea", "ara": "Arabidopsis", "ery": "Erythroxylum"}
MIN_SIZE = 10000000 

# 1. Pre-create/touch all expected links tracks to satisfy Circos file constraints
for sbj in ["ara", "pop", "man", "hev", "ery"]:
    open(os.path.join(circ_dir, "links", f"hyp_{sbj}.links"), "w").close()

# 2. Extract accurate coordinates from BED tracking metrics
chrom_lengths = {}
for pfx in species_order:
    bed_path = os.path.join(clean_dir, f"{pfx}.bed")
    if os.path.exists(bed_path):
        for line in open(bed_path):
            parts = line.split()
            if len(parts) >= 3:
                chrom_lengths[(pfx, parts[0])] = max(chrom_lengths.get((pfx, parts[0]), 0), int(parts[2]))

# Update Hypericum coordinates from true assembly basepair counts
cur_chrom, cur_len = None, 0
for line in open(hyp_fasta):
    if line.startswith(">"):
        if cur_chrom:
            n = re.search(r"\d+", cur_chrom)
            chrom_lengths[("hyp", f"HpChr{int(n.group()):02d}" if n else cur_chrom)] = cur_len
        cur_chrom = line.strip().split()[0][1:]
        cur_len = 0
    else: cur_len += len(line.strip())
if cur_chrom:
    n = re.search(r"\d+", cur_chrom)
    chrom_lengths[("hyp", f"HpChr{int(n.group()):02d}" if n else cur_chrom)] = cur_len

# 3. Output structural Karyotype index 
with open(os.path.join(circ_dir, "karyotype.txt"), "w") as out:
    for pfx in species_order:
        prefix_chars = PREFIX_MAP.get(pfx, "")
        for (sp, chrom), length in sorted(chrom_lengths.items()):
            if sp == pfx and length >= MIN_SIZE and chrom.startswith(prefix_chars + "Chr"):
                out.write(f"chr - {sp}_{chrom} {chrom} 0 {length} {species_info[sp]}\n")

# 4. Generate Link Ribbon Intersections (BUGFIX: removed destructive ID replacements)
def load_bed_lookup(pfx):
    d = {}
    path = os.path.join(clean_dir, f"{pfx}.bed")
    if not os.path.exists(path): return d
    for line in open(path):
        p = line.split()
        if len(p) >= 4: d[p[3]] = (p[0], int(p[1]), int(p[2]))
    return d

beds = {s: load_bed_lookup(s) for s in species_order}

for f in os.listdir(syn_dir):
    if not f.endswith(".simple"): continue
    ref_sp, qry_sp = f.split(".")[0], f.split(".")[1]
    with open(os.path.join(circ_dir, "links", f"{ref_sp}_{qry_sp}.links"), "w") as out:
        for line in open(os.path.join(syn_dir, f)):
            if line.startswith("#"): continue
            p = line.split()
            if p[0] in beds[ref_sp] and p[2] in beds[qry_sp]:
                chrA, sA1, eA1 = beds[ref_sp][p[0]]; chrB, sB1, eB1 = beds[qry_sp][p[2]]
                _, sA2, eA2 = beds[ref_sp][p[1]]; _, sB2, eB2 = beds[qry_sp][p[3]]
                if chrA.startswith(PREFIX_MAP[ref_sp]) and chrB.startswith(PREFIX_MAP[qry_sp]):
                    out.write(f"{ref_sp}_{chrA} {min(sA1,sA2,eA1,eA2)} {max(sA1,sA2,eA1,eA2)} {qry_sp}_{chrB} {min(sB1,sB2,eB1,eB2)} {max(sB1,sB2,eB1,eB2)}\n")

# 5. Format centered visual text titles safely
sp_chroms = defaultdict(list)
kary_path = os.path.join(circ_dir, "karyotype.txt")
if os.path.exists(kary_path):
    for line in open(kary_path):
        p = line.split()
        if len(p) >= 3: sp_chroms[p[2].split("_")[0]].append(p)

with open(os.path.join(circ_dir, "species_labels.txt"), "w") as out:
    for s_code, chrs in sp_chroms.items():
        if s_code in species_names:
            mid = chrs[len(chrs)//2]
            mid_pos = int(mid[5]) // 2
            out.write(f"{mid[2]} {mid_pos} {mid_pos} {species_names[s_code]}\n")
