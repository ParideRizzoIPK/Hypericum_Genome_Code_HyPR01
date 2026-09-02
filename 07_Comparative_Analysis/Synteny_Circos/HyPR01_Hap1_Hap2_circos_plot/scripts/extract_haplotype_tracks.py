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
