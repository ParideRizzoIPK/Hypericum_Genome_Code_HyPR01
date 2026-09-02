#!/bin/bash
#SBATCH --job-name=HyPR01_DualHap_EggNOG
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/EggNOG_DualHap_%j.out
#SBATCH --error=logs/EggNOG_DualHap_%j.err

set -e

# ────────────────────────────────────────────────────────
# 1. MODULE LOADING & ENVIRONMENT
# ────────────────────────────────────────────────────────
. /etc/profile.d/modules.sh
module load gffread/0.12.6
module load eggnog-mapper/2.1.12
module load python/3.12.0

EGGNOG_DATA_DIR="/path/to/your/directory/eggnog_data"
export EGGNOG_DATA_DIR

# --- MAPPED BIOINFORMATICS ASSEMBLY PATHS ---
GENOME_HAP1="/path/to/input/directory/hap1.masked_unchr.fa"
GENOME_HAP2="/path/to/input/directory/hap2.masked_unchr.fa"

BASE_OUT_DIR="/path/to/output/directory"

# ────────────────────────────────────────────────────────
# 2. EMBEDDED ATTR MERGER SCRIPT
# ────────────────────────────────────────────────────────
PY_SCRIPT="annotate_gff_embedded.py"
cat << 'EOF' > "$PY_SCRIPT"
import csv
import sys
import re

if len(sys.argv) != 4:
    print("Usage: python script.py <input_gff> <eggnog_file> <output_gff>")
    sys.exit(1)

INPUT_GFF, EGGNOG_FILE, OUTPUT_GFF = sys.argv[1:4]

def clean_gff_value(value):
    if not value or value == "-": return None
    return value.replace(';', ',').replace('\n', ' ').replace('\r', '').replace('=', '%3D')

def load_eggnog_data(filepath):
    data = {}
    with open(filepath, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row or row[0].startswith("#"): continue
            try:
                query_id, desc, gene_name, go, ec = row[0], row[7], row[8], row[9], row[10]
                data[query_id] = {
                    'note': clean_gff_value(desc),
                    'gene': clean_gff_value(gene_name),
                    'Ontology_term': clean_gff_value(go),
                    'EC_number': clean_gff_value(ec)
                }
            except IndexError: continue
    return data

def process_gff(input_gff, output_gff, annotation_map):
    match_count = 0
    with open(input_gff, 'r') as fin, open(output_gff, 'w') as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            parts = line.strip().split('\t')
            if len(parts) != 9:
                fout.write(line)
                continue
            
            match = re.search(r'ID=([^;]+)', parts[8])
            if match:
                current_id = match.group(1)
                if current_id in annotation_map:
                    match_count += 1
                    ann = annotation_map[current_id]
                    new_attrs = []
                    if ann['note']: new_attrs.extend([f";Note={ann['note']}", f";product={ann['note']}"])
                    if ann['gene']: new_attrs.append(f";gene={ann['gene']}")
                    if ann['Ontology_term']: new_attrs.append(f";Ontology_term={ann['Ontology_term']}")
                    if ann['EC_number']: new_attrs.append(f";EC_number={ann['EC_number']}")
                    parts[8] += "".join(new_attrs)
            fout.write('\t'.join(parts) + '\n')
    print(f"    -> Successfully appended entries to {match_count} transcript features.")

if __name__ == "__main__":
    annot_map = load_eggnog_data(EGGNOG_FILE)
    process_gff(INPUT_GFF, OUTPUT_GFF, annot_map)
EOF

# ────────────────────────────────────────────────────────
# 3. RUN LOOP PROCESSING STRATEGY
# ────────────────────────────────────────────────────────
declare -A HAP_CONFIGS
HAP_CONFIGS[Hap1]="$GENOME_HAP1"
HAP_CONFIGS[Hap2]="$GENOME_HAP2"

for HAP in Hap1 Hap2; do
    GENOME_FASTA=${HAP_CONFIGS[$HAP]}
    TARGET_DIR="${BASE_OUT_DIR}/HyPR01_${HAP}"
    INPUT_GFF="${TARGET_DIR}/HyPR01_${HAP}_final_annotations.gff3"
    EXTRACTED_FAA="${TARGET_DIR}/HyPR01_${HAP}_extracted_proteins.faa"
    CLEANED_FAA="${TARGET_DIR}/HyPR01_${HAP}_extracted_proteins_cleaned.faa"
    OUT_PREFIX="HyPR01_${HAP}_emapper"
    FINAL_ANNOTATED_GFF="${TARGET_DIR}/HyPR01_${HAP}_final_annotations_functional.gff3"
    
    # Temporary renamed fasta for gffread compatibility
    TEMP_FASTA="${TARGET_DIR}/temp_genome_sorted_headers.fa"

    echo "========================================================================"
    echo "[$(date)] STARTING WORKFLOW FOR RUN TRACK: $HAP"
    echo "========================================================================"
    echo "Workspace Target:  $TARGET_DIR"
    echo "Genome Reference:  $GENOME_FASTA"
    echo "Coordinate Source: $INPUT_GFF"
    
    if [ ! -f "$INPUT_GFF" ]; then
        echo "[ERROR] Missing input coordinates: $INPUT_GFF. Skipping Haplotype."
        continue
    fi

    # Phase A: Fix Fasta Headers and Run gffread structural peptide extraction
    echo "[-] Step 1: Normalizing fasta headers and running gffread translation..."
    sed 's/^>UnChr/>ChrUn/' "$GENOME_FASTA" > "$TEMP_FASTA"
    
    # Adding -V flag drops trailing stop codons automatically from transcripts
    gffread "$INPUT_GFF" -g "$TEMP_FASTA" -V -y "$EXTRACTED_FAA"
    rm -f "$TEMP_FASTA"

    # Strict secondary FASTA cleanup pass: Strip out any remaining '.' or '*' from the sequence lines
    # This leaves the sequence text intact but sanitizes characters that break DIAMOND
    sed '/^[^>]/s/[\.\*]//g' "$EXTRACTED_FAA" > "$CLEANED_FAA"
    rm -f "$EXTRACTED_FAA"

    # Phase B: Quantify calculated sequence counts
    PROTEIN_COUNT=$(grep -c "^>" "$CLEANED_FAA" || true)
    echo "[+] Extraction complete. Total proteins available for analysis: $PROTEIN_COUNT"

    if [ "$PROTEIN_COUNT" -eq 0 ]; then
        echo "[ERROR] Fasta file is empty. Check assembly alignment compatibility. Skipping Haplotype."
        continue
    fi

    # Phase C: eggNOG-mapper search run using Viridiplantae (tax_scope 33090)
    echo "[-] Step 2: Running eggNOG-mapper orthology mapping loop..."
    emapper.py \
        -i "$CLEANED_FAA" \
        --itype proteins \
        -o "$OUT_PREFIX" \
        --output_dir "$TARGET_DIR" \
        --cpu "$SLURM_CPUS_PER_TASK" \
        -m diamond \
        --sensmode more-sensitive \
        --tax_scope 33090 \
        --override

    # Phase D: Execute embedded python structural column append
    echo "[-] Step 3: Injecting functional metadata terms into target file..."
    EGG_OUTPUT_MAP="${TARGET_DIR}/${OUT_PREFIX}.emapper.annotations"

    if [ -f "$EGG_OUTPUT_MAP" ]; then
        python3 "$PY_SCRIPT" "$INPUT_GFF" "$EGG_OUTPUT_MAP" "$FINAL_ANNOTATED_GFF"
        echo "[+] Production functional model saved to: $FINAL_ANNOTATED_GFF"
    else
        echo "[ERROR] eggNOG annotations table output missed for path: $EGG_OUTPUT_MAP"
    fi
    
    # Final cleanup of haplotype-specific intermediate fasta files
    rm -f "$CLEANED_FAA"
    echo ""
done

# Clear cache script file
rm -f "$PY_SCRIPT"
echo "[*] Dual Haplotype Annotation Run finished completely at $(date)"