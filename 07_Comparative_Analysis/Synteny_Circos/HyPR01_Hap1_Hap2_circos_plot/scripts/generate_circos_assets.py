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
