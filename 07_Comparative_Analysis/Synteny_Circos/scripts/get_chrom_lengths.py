import sys, os, re

prefix = sys.argv[1]
fasta_path = sys.argv[2]
out_path = sys.argv[3]

with open(out_path, "a") as out:
    cur_chrom, cur_len = None, 0
    for line in open(fasta_path):
        if line.startswith(">"):
            if cur_chrom:
                if prefix == "hyp":
                    n = re.search(r"\d+", cur_chrom)
                    chrom_name = f"HpChr{int(n.group()):02d}" if n else cur_chrom
                else:
                    chrom_name = cur_chrom
                out.write(f"{prefix}\t{chrom_name}\t{cur_len}\n")
            cur_chrom = line.strip().split()[0][1:]
            cur_len = 0
        else:
            cur_len += len(line.strip())
    if cur_chrom:
        if prefix == "hyp":
            n = re.search(r"\d+", cur_chrom)
            chrom_name = f"HpChr{int(n.group()):02d}" if n else cur_chrom
        else:
            chrom_name = cur_chrom
        out.write(f"{prefix}\t{chrom_name}\t{cur_len}\n")
