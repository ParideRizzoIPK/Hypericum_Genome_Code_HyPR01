#!/bin/bash
#SBATCH --job-name=Malp_Hap2_Master
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=45G
#SBATCH --time=24:00:00
#SBATCH --output=logs/master_pipeline_%j.out
#SBATCH --error=logs/master_pipeline_%j.err

# Fail fast on any error, unset variable, or pipe failure
set -euo pipefail

# ==============================================================================
# 1. PATH & ARCHITECTURE DEFINITIONS
# ==============================================================================
BASE_DIR="/path/to/your/directory"
GENOMES_DIR="$BASE_DIR/genomes"
CLEAN_DIR="$BASE_DIR/cleaned_data"
SYN_DIR="$BASE_DIR/mcscanx"
CIRC_DIR="$BASE_DIR/circos"
SCRIPTS_DIR="$BASE_DIR/scripts"
LOG_DIR="$BASE_DIR/logs"

TRACKING_FILE="$BASE_DIR/Chromosome_Mapping_Table.tsv"

# Enforce clean sandbox subdirectory tree creation
mkdir -p "$GENOMES_DIR" "$CLEAN_DIR" "$SYN_DIR" "$CIRC_DIR/links" "$SCRIPTS_DIR" "$LOG_DIR"

# ENVIRONMENT INJECTION: Prepend environment binaries directly to PATH.
export PATH="/path/to/your/directory/micromamba/envs/synteny_circos/bin:/path/to/your/directory/micromamba/envs/last_env/bin:$PATH"
export SINGULARITY_BIND="/path/to/your/directory"

# Define raw input target references
NEW_HYP_FASTA="$BASE_DIR/hap2.masked_unchr.fa"
NEW_HYP_GFF="$BASE_DIR/harmonized_consensus_Hap2.gff3"
NEW_HYP_PEP="$BASE_DIR/hyp.pep"

echo "==================================================="
echo " PHASE 1: RECONSTRUCTING GENOMES & ASSET STAGING   "
echo "==================================================="
if [ ! -f "$NEW_HYP_GFF" ] || [ ! -f "$NEW_HYP_FASTA" ]; then
    echo "ERROR: Mandatory Hypericum input files missing from $BASE_DIR!" >&2
    exit 1
fi

cd "$GENOMES_DIR"
if [ ! -f "datasets" ]; then
    echo "Staging NCBI CLI datasets tool..."
    curl -sS -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets'
    chmod +x datasets
fi

# -- Ensembl Plants Bulk Engine (Theobroma cacao removed) --
CNSEMB_SPECIES=("arabidopsis_thaliana" "populus_trichocarpa" "manihot_esculenta")
FOLDER_NAMES=("arabidopsis_col0" "populus_trichocarpa" "manihot_esculenta")

for i in "${!CNSEMB_SPECIES[@]}"; do
    sp="${CNSEMB_SPECIES[$i]}"
    folder="${FOLDER_NAMES[$i]}"
    echo "--> Downloading Ensembl core package: $folder"
    mkdir -p "$GENOMES_DIR/$folder" && cd "$GENOMES_DIR/$folder"
    
    wget -q ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/$sp/dna/*.dna.toplevel.fa.gz
    wget -q ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/$sp/pep/*.pep.all.fa.gz
    wget -q ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/gff3/$sp/*.gff3.gz
    
    rm -f *chromosome*.gff3.gz *primary_assembly*.gff3.gz *chr.gff3.gz *v6.dna.toplevel.fa.gz
    gunzip -f *.gz
    rm -f *chromosome*.gff3 *primary_assembly*.gff3 *chr.gff3
done

# -- NCBI RefSeq / GenBank Assembly Fetching --
echo "--> Pulling Hevea brasiliensis from NCBI datasets..."
mkdir -p "$GENOMES_DIR/hevea_brasiliensis" && cd "$GENOMES_DIR/hevea_brasiliensis"
"$GENOMES_DIR/datasets" download genome accession GCF_030052815.1 --include gff3,protein,genome --filename hevea_dataset.zip
unzip -qo hevea_dataset.zip
mv ncbi_dataset/data/GCF_030052815.1/*.fna Hevea_brasiliensis.genome.fasta
mv ncbi_dataset/data/GCF_030052815.1/*.gff Hevea_brasiliensis.gff3
mv ncbi_dataset/data/GCF_030052815.1/protein.faa Hevea_brasiliensis.proteome.fasta
rm -rf ncbi_dataset hevea_dataset.zip README.md md5sum.txt

echo "--> Pulling Erythroxylum novogranatense from NCBI datasets..."
mkdir -p "$GENOMES_DIR/erythroxylum_novogranatense" && cd "$GENOMES_DIR/erythroxylum_novogranatense"
"$GENOMES_DIR/datasets" download genome accession PRJNA749480 --include gff3,protein,genome --filename erythroxylum_dataset.zip
unzip -qo erythroxylum_dataset.zip
mv ncbi_dataset/data/GCA_029891385.1/*.fna Erythroxylum_novogranatense.genome.fasta
mv ncbi_dataset/data/GCA_029891385.1/*.gff Erythroxylum_novogranatense.gff3
mv ncbi_dataset/data/GCA_029891385.1/protein.faa Erythroxylum_novogranatense.proteome.fasta
rm -rf ncbi_dataset erythroxylum_dataset.zip README.md md5sum.txt

# -- Processing Local Hypericum & Quality-Filtering Proteins --
echo "--> Compiling high-quality complete CDS protein tracks for Hypericum..."

cat << 'EOF' > "$SCRIPTS_DIR/extract_clean_proteins.py"
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

fasta_path, gff_path, out_pep, out_bed = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

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
            chrom_clean = f"HpChr{int(num_match.group()):02d}" if num_match else chrom
            bed_out.write(f"{chrom_clean}\t{start-1}\t{end}\t{tx_id}\t0\t{strand}\n")
EOF

python3 "$SCRIPTS_DIR/extract_clean_proteins.py" "$NEW_HYP_FASTA" "$NEW_HYP_GFF" "$CLEAN_DIR/hyp.pep" "$CLEAN_DIR/hyp.bed"
cp "$CLEAN_DIR/hyp.pep" "$NEW_HYP_PEP"

# ==============================================================================
# PHASE 2: CHROMOSOME HARMONIZATION ENGINE         
# ==============================================================================
echo "==================================================="
echo " PHASE 2: CHROMOSOME HARMONIZATION ENGINE         "
echo "==================================================="
cd "$BASE_DIR"

cat << 'EOF' > "$SCRIPTS_DIR/rename_bed_chromosomes.py"
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
EOF

echo -e "Species\tRaw_ID\tHarmonized_ID" > "$TRACKING_FILE"
SPECIES_ARRAY=("ara:arabidopsis_col0" "pop:populus_trichocarpa" "man:manihot_esculenta" "hev:hevea_brasiliensis" "ery:erythroxylum_novogranatense")

for ENTRY in "${SPECIES_ARRAY[@]}"; do
    IFS=":" read -r PREFIX FOLDER <<< "$ENTRY"
    
    cp "/path/to/input/directory/$PREFIX.bed" "$CLEAN_DIR/" 2>/dev/null || true
    cp "/path/to/your/directory" "$CLEAN_DIR/" 2>/dev/null || true
    
    GFF_FILE=$(find "$GENOMES_DIR/$FOLDER" -maxdepth 1 -name "*.gff3" 2>/dev/null | head -n 1) || true
    if [ -n "$GFF_FILE" ]; then
        python3 "$SCRIPTS_DIR/rename_bed_chromosomes.py" "$PREFIX" "$GFF_FILE" "$CLEAN_DIR/$PREFIX.bed" "$TRACKING_FILE"
    fi
done

# AUTOMATED MANIHOT ID HARMONIZATION FILTER: Strips '.CDS' extensions from BOTH files to force a byte-for-byte structural match
python3 -c '
import os
clean_dir = "'$CLEAN_DIR'"

# 1. Standardize man.bed 4th column entries
bed_path = os.path.join(clean_dir, "man.bed")
if os.path.exists(bed_path):
    lines = []
    for line in open(bed_path):
        parts = line.strip().split("\t")
        if len(parts) >= 4:
            parts[3] = parts[3].split(".CDS")[0]
        lines.append("\t".join(parts) + "\n")
    with open(bed_path, "w") as f: f.writelines(lines)

# 2. Standardize man.pep sequence headers
pep_path = os.path.join(clean_dir, "man.pep")
if os.path.exists(pep_path):
    lines = []
    for line in open(pep_path):
        if line.startswith(">"):
            token = line.strip().split()[0]
            line = token.split(".CDS")[0] + "\n"
        lines.append(line)
    with open(pep_path, "w") as f: f.writelines(lines)
'

echo "==================================================="
echo " PHASE 3: CONCURRENT GEOMETRIC SYNTENY ALIGNMENTS  "
echo "==================================================="
cd "$SYN_DIR"
ln -sf "$CLEAN_DIR"/*.pep .
ln -sf "$CLEAN_DIR"/*.bed .

# CLEAR ALIGNMENT CACHE: Deletes the old index traps to guarantee a fresh alignment calculation
/bin/rm -f hyp.man.last hyp.man.blast man.idx hyp.idx

SUBJECTS=("ara" "pop" "man" "hev" "ery")
for SBJ in "${SUBJECTS[@]}"; do
    echo "⚡ Spawning async background threads: hyp vs $SBJ"
    python3 -m jcvi.compara.catalog ortholog hyp $SBJ --dbtype=prot --cpus=3 --cscore=.7 --no_strip_names > "../logs/jcvi_${SBJ}.log" 2>&1 && \
    python3 -m jcvi.compara.synteny screen hyp.$SBJ.anchors hyp.$SBJ.anchors.simple --minspan=5 --simple >> "../logs/jcvi_${SBJ}.log" 2>&1 &
done

echo "⏳ Holding main core till concurrent tasks complete..."
wait

echo "==================================================="
echo " PHASE 4: EXACT CHROMOSOME SCALING & CIRCOS MAPS  "
echo "==================================================="
cd "$BASE_DIR"

/bin/rm -f "$CLEAN_DIR/chrom_sizes.tsv"
touch "$CLEAN_DIR/chrom_sizes.tsv"

cat << 'EOF' > "$SCRIPTS_DIR/get_chrom_lengths.py"
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
EOF

python3 "$SCRIPTS_DIR/get_chrom_lengths.py" "hyp" "$NEW_HYP_FASTA" "$CLEAN_DIR/chrom_sizes.tsv"

SPECIES_GENOMES=(
    "ara:arabidopsis_col0/*.dna.toplevel.fa"
    "pop:populus_trichocarpa/*.dna.toplevel.fa"
    "man:manihot_esculenta/*.dna.toplevel.fa"
    "hev:hevea_brasiliensis/*genome.fasta"
    "ery:erythroxylum_novogranatense/*genome.fasta"
)
for ENTRY in "${SPECIES_GENOMES[@]}"; do
    IFS=":" read -r PFX GLOB <<< "$ENTRY"
    FASTA_FILE=$(find "$GENOMES_DIR" -maxdepth 2 -path "$GENOMES_DIR/$GLOB" | head -n 1) || true
    if [ -n "$FASTA_FILE" ]; then
        python3 "$SCRIPTS_DIR/get_chrom_lengths.py" "$PFX" "$FASTA_FILE" "$CLEAN_DIR/chrom_sizes.tsv"
    fi
done

python3 -c '
import os, sys, re
clean_dir, circ_dir, tracking_file = sys.argv[1], sys.argv[2], sys.argv[3]
mapping = {}
if os.path.exists(tracking_file):
    for line in open(tracking_file):
        p = line.strip().split("\t")
        if len(p) >= 3: mapping[(p[0], p[1])] = p[2]

species_order = ["hyp", "pop", "man", "hev", "ara", "ery"]
species_info = {"hyp": "102,0,153", "pop": "255,127,80", "man": "241,196,15", "hev": "46,204,113", "ara": "112,128,144", "ery": "0,191,255"}
PREFIX_MAP = {"ara":"At", "pop":"Pt", "man":"Me", "hev":"Hb", "ery":"En", "hyp":"Hp"}
MIN_SIZE = 10000000 

with open(f"{circ_dir}/karyotype.txt", "w") as out:
    for line in open(f"{clean_dir}/chrom_sizes.tsv", "r"):
        p = line.strip().split("\t")
        if len(p) != 3: continue
        sp, raw_chrom, length = p[0], p[1], int(p[2])
        harm = raw_chrom if sp == "hyp" else mapping.get((sp, raw_chrom))
        if not harm:
            n = re.sub(r"[^\d]", "", raw_chrom)
            if n.isdigit() and len(n) <= 2: harm = f"{PREFIX_MAP.get(sp, sp.capitalize())}Chr{int(n):02d}"
        if harm:
            prefix_chars = PREFIX_MAP.get(sp, "")
            if length >= MIN_SIZE and harm.startswith(f"{prefix_chars}Chr"):
                out.write(f"chr - {sp}_{harm} {harm} 0 {length} {species_info[sp]}\n")
' "$CLEAN_DIR" "$CIRC_DIR" "$TRACKING_FILE"

python3 -c '
import os, sys
from collections import defaultdict

clean_dir = sys.argv[1]
syn_dir = sys.argv[2]
circ_dir = sys.argv[3]

for sbj in ["ara", "pop", "man", "hev", "ery"]:
    open(f"{circ_dir}/links/hyp_{sbj}.links", "w").close()

def load_bed(pfx):
    d = {}
    if not os.path.exists(f"{clean_dir}/{pfx}.bed"): return d
    for line in open(f"{clean_dir}/{pfx}.bed"):
        p = line.split()
        if len(p) >= 4: d[p[3]] = (p[0], int(p[1]), int(p[2]))
    return d

species = ["hyp", "ara", "pop", "man", "hev", "ery"]
species_names = {"hyp": "Hypericum_Hap2", "pop": "Populus", "man": "Manihot", "hev": "Hevea", "ara": "Arabidopsis", "ery": "Erythroxylum"}
PREFIX_MAP = {"ara":"At", "pop":"Pt", "man":"Me", "hev":"Hb", "ery":"En", "hyp":"Hp"}
beds = {s: load_bed(s) for s in species}

for f in os.listdir(syn_dir):
    if not f.endswith(".simple"): continue
    ref_sp, qry_sp = f.split(".")[0], f.split(".")[1]
    with open(f"{circ_dir}/links/{ref_sp}_{qry_sp}.links", "w") as out:
        for line in open(f"{syn_dir}/{f}"):
            if line.startswith("#"): continue
            p = line.split()
            if p[0] in beds[ref_sp] and p[2] in beds[qry_sp]:
                chrA, sA1, eA1 = beds[ref_sp][p[0]]; chrB, sB1, eB1 = beds[qry_sp][p[2]]
                _, sA2, eA2 = beds[ref_sp][p[1]]; _, sB2, eB2 = beds[qry_sp][p[3]]
                if chrA.startswith(PREFIX_MAP[ref_sp]) and chrB.startswith(PREFIX_MAP[qry_sp]):
                    out.write(f"{ref_sp}_{chrA} {min(sA1,sA2,eA1,eA2)} {max(sA1,sA2,eA1,eA2)} {qry_sp}_{chrB} {min(sB1,sB2,eB1,eB2)} {max(sB1,sB2,eB1,eB2)}\n")

sp_chroms = defaultdict(list)
for line in open(f"{circ_dir}/karyotype.txt"):
    p = line.split()
    if len(p) >= 3:
        sp_chroms[p[2].split("_")[0]].append(p)

with open(f"{circ_dir}/species_labels.txt", "w") as out:
    for s_code, chrs in sp_chroms.items():
        if s_code in species_names:
            mid = chrs[len(chrs)//2]
            mid_pos = int(mid[5]) // 2
            out.write(f"{mid[2]} {mid_pos} {mid_pos} {species_names[s_code]}\n")
' "$CLEAN_DIR" "$SYN_DIR" "$CIRC_DIR"

cd "$CIRC_DIR"

cat << 'EOF' > ideogram.conf
<ideogram>
<spacing>
default = 0.005r
<pairwise hyp pop>
spacing = 0.10r
</pairwise>
<pairwise pop man>
spacing = 0.10r
</pairwise>
<pairwise man hev>
spacing = 0.10r
</pairwise>
<pairwise hev ara>
spacing = 0.10r
</pairwise>
<pairwise ara ery>
spacing = 0.10r
</pairwise>
<pairwise ery hyp>
spacing = 0.10r
</pairwise>
</spacing>
radius           = 0.70r
thickness        = 40p
fill             = yes
stroke_color     = dgrey
stroke_thickness = 2p
show_label       = yes
label_font       = bold
label_radius     = 1.16r
label_size       = 40
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
size           = 15p
show_label     = yes
label_size     = 20p
label_offset   = 5p
format         = %d
</tick>
</ticks>
EOF

cat << 'EOF' > circos.conf
karyotype = karyotype.txt
chromosomes_units = 1000000
chromosomes_display_default = yes

<colors>
hyp_color = 102,0,153
ara_color = 112,128,144
pop_color = 255,127,80
man_color = 241,196,15
hev_color = 46,204,113
ery_color = 0,191,255
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
condition = var(size1) < 250000 || var(size2) < 250000
show = no
</rule>
</rules>

<link ara>
file = links/hyp_ara.links
color = ara_color_a4
</link>
<link pop>
file = links/hyp_pop.links
color = pop_color_a4
</link>
<link man>
file = links/hyp_man.links
color = man_color_a4
</link>
<link hev>
file = links/hyp_hev.links
color = hev_color_a4
</link>
<link ery>
file = links/hyp_ery.links
color = ery_color_a4
</link>
</links>

<plots>
<plot>
type             = text
color            = black
file             = species_labels.txt
r0               = 1.10r
r1               = 1.40r
show_links       = no
label_size       = 70p
label_font       = bold
</plot>
</plots>

<image>
<<include etc/image.conf>>
dir* = .
file* = Malpighiales_Synteny_Final
image_format* = png,svg
radius* = 1500p
background* = white
</image>

<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
EOF

echo "[$(date)] Launching main multi-format vector compilation..."
circos -conf circos.conf

echo "==================================================="
echo " MASTER ARCHITECTURE RUN COMPLETED SUCCESSFULLY    "
echo " PNG Render: $CIRC_DIR/Malpighiales_Synteny_Final.png"
echo " SVG Render: $CIRC_DIR/Malpighiales_Synteny_Final.svg"
echo "==================================================="