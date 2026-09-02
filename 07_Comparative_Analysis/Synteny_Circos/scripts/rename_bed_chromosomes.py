import sys, os, re
species, gff_file, bed_file, tracking_file = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
mapping = {}
prefix_map = {'ara': 'At', 'pop': 'Pt', 'man': 'Me', 'hev': 'Hb', 'ery': 'En', 'hyp': 'Hp'}
sp_prefix = prefix_map.get(species, species.capitalize())

with open(gff_file, "r") as f:
    for line in f:
        if line.startswith("#"): continue
        cols = line.strip().split("\t")
        if len(cols) != 9: continue
        if cols[2] in ["region", "chromosome"]:
            raw_id, attrs_str = cols[0], cols[8]
            attrs = {item.split('=', 1)[0]: item.split('=', 1)[1] for item in attrs_str.split(";") if "=" in item}
            clean_num = None
            if "chromosome" in attrs: clean_num = attrs["chromosome"]
            elif "Name" in attrs and "LG" in attrs["Name"]: clean_num = attrs["Name"].replace("LG", "")
            elif "Alias" in attrs:
                for a in attrs["Alias"].split(","):
                    if "Chromosome" in a or "chr" in a.lower():
                        clean_num = re.sub(r'[^\d]', '', a); break
            elif "ID" in attrs and "chromosome:" in attrs["ID"]: clean_num = attrs["ID"].replace("chromosome:", "")
            elif raw_id.isdigit() or raw_id.lower().startswith("chr"): clean_num = re.sub(r'[^\d]', '', raw_id)
            
            if clean_num:
                num_only = re.sub(r'[^\d]', '', str(clean_num))
                if num_only.isdigit():
                    harmonized_id = f"{sp_prefix}Chr{int(num_only):02d}"
                    if raw_id not in mapping:
                        mapping[raw_id] = harmonized_id
                        with open(tracking_file, "a") as tf: tf.write(f"{species}\t{raw_id}\t{harmonized_id}\n")

if os.path.exists(bed_file):
    new_bed_lines = []
    with open(bed_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                raw_chrom = parts[0]
                if raw_chrom not in mapping:
                    num_only = re.sub(r'[^\d]', '', raw_chrom)
                    if num_only.isdigit() and len(num_only) <= 2:
                        mapping[raw_chrom] = f"{sp_prefix}Chr{int(num_only):02d}"
                parts[0] = mapping.get(raw_chrom, raw_chrom)
                new_bed_lines.append("\t".join(parts) + "\n")
    with open(bed_file, "w") as f: f.writelines(new_bed_lines)
