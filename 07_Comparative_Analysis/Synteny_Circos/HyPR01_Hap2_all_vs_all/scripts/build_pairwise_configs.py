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
