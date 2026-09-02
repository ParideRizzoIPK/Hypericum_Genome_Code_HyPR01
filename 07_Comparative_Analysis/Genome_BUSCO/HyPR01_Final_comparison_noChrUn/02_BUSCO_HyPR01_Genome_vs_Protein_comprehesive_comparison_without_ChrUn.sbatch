#!/bin/bash
#
#SBATCH --job-name=HyPR01_BUSCO_noChrUn
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=logs/busco_final_%j.out
#SBATCH --error=logs/busco_final_%j.err

# Strict error handling: exit on command failure (-e), unset variables (-u), and pipeline failure (-o pipefail)
set -euo pipefail

# ==============================================================================
# 1. ENVIRONMENT & PATH SETUP
# ==============================================================================
export SINGULARITY_BIND="/path/to/your/directory"

# Load modules from the cluster environment
module load gffread/0.12.6
module load BUSCO/5.8.2
module load python/3.12.0
module load samtools/1.23.1

# Master Output Directory (noChrUn)
OUT_DIR="/path/to/output/directory"
GENOME_DIR="${OUT_DIR}/filtered_genomes"
GFF_DIR="${OUT_DIR}/filtered_gffs"
PROTEIN_DIR="${OUT_DIR}/extracted_proteins"
DOWNLOAD_DIR="/path/to/your/directory/busco_downloads"

# Raw Source Paths (Genome FASTAs)
RAW_HAP1_GENOME="/path/to/input/directory/hap1.masked_ChrUn.fa"
RAW_HAP2_GENOME="/path/to/input/directory/hap2.masked_ChrUn.fa"

# Raw Source Paths (Updated Functional Annotation GFF3s)
RAW_HAP1_GFF="/path/to/input/directory/HyPR01_Hap1_final_annotations_functional.gff3"
RAW_HAP2_GFF="/path/to/input/directory/HyPR01_Hap2_final_annotations_functional.gff3"

# Filtered Output Paths (Nuclear-Only: Chr01 - Chr08)
HAP1_GENOME="${GENOME_DIR}/hap1_nuclear_Chr01-Chr08.fa"
HAP2_GENOME="${GENOME_DIR}/hap2_nuclear_Chr01-Chr08.fa"

HAP1_GFF="${GFF_DIR}/hap1_nuclear_Chr01-Chr08.gff3"
HAP2_GFF="${GFF_DIR}/hap2_nuclear_Chr01-Chr08.gff3"

# Initialize directories
mkdir -p "${OUT_DIR}"
mkdir -p "${GENOME_DIR}"
mkdir -p "${GFF_DIR}"
mkdir -p "${PROTEIN_DIR}"

# Input validation step
echo "[$(date)] ==== INPUT VALIDATION ===="
for file in "${RAW_HAP1_GENOME}" "${RAW_HAP2_GENOME}" "${RAW_HAP1_GFF}" "${RAW_HAP2_GFF}"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required input file does not exist: $file" >&2
        exit 1
    fi
done
echo "All input raw genome and functional annotation GFF3 files verified."

# ==============================================================================
# 2. FASTA & GFF3 FILTERING: REMOVE ChrUn (KEEP ONLY Chr01 - Chr08)
# ==============================================================================
echo -e "\n[$(date)] ==== STEP 0: FILTERING GENOMES & GFF3 FILES ===="

# Function to extract nuclear chromosomes from FASTA
extract_nuclear_chromosomes () {
    local input_fa="$1"
    local output_fa="$2"
    
    echo "Indexing ${input_fa}..."
    samtools faidx "${input_fa}"
    
    echo "Extracting Chr01 to Chr08 into ${output_fa}..."
    samtools faidx "${input_fa}" Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 > "${output_fa}"
    samtools faidx "${output_fa}"
}

# Function to remove ChrUn records and pragma headers from GFF3
filter_gff_no_chrun () {
    local input_gff="$1"
    local output_gff="$2"
    
    echo "Filtering ChrUn out of ${input_gff}..."
    grep -v -E '^##sequence-region[[:space:]]+ChrUn' "${input_gff}" | grep -v -E '^ChrUn[[:space:]]' > "${output_gff}"
}

echo "[$(date)] Filtering Haplotype 1 Genome & GFF3..."
extract_nuclear_chromosomes "${RAW_HAP1_GENOME}" "${HAP1_GENOME}"
filter_gff_no_chrun "${RAW_HAP1_GFF}" "${HAP1_GFF}"

echo "[$(date)] Filtering Haplotype 2 Genome & GFF3..."
extract_nuclear_chromosomes "${RAW_HAP2_GENOME}" "${HAP2_GENOME}"
filter_gff_no_chrun "${RAW_HAP2_GFF}" "${HAP2_GFF}"

echo "[$(date)] Verification - Chromosomes remaining in Hap1 GFF3:"
awk -F'\t' '!/^#/ {print $1}' "${HAP1_GFF}" | sort | uniq

echo "[$(date)] Verification - Chromosomes remaining in Hap2 GFF3:"
awk -F'\t' '!/^#/ {print $1}' "${HAP2_GFF}" | sort | uniq

# ==============================================================================
# 3. STEP 1: GFF3-TO-PROTEIN TRANSLATION (gffread)
# ==============================================================================
echo -e "\n[$(date)] ==== STEP 1: EXTRACT PROTEINS WITH GFFREAD ===="

HAP1_FAA="${PROTEIN_DIR}/hap1_proteins.faa"
HAP2_FAA="${PROTEIN_DIR}/hap2_proteins.faa"

# Haplotype 1 Protein Extraction using filtered FASTA and functional GFF3
echo "[$(date)] Extracting and translating Haplotype 1 proteins..."
gffread -y "${HAP1_FAA}" -g "${HAP1_GENOME}" "${HAP1_GFF}"
HAP1_COUNT=$(grep -c '^>' "${HAP1_FAA}")
echo "[$(date)] Verification: Extracted ${HAP1_COUNT} nuclear protein sequences for Haplotype 1."

# Haplotype 2 Protein Extraction using filtered FASTA and functional GFF3
echo "[$(date)] Extracting and translating Haplotype 2 proteins..."
gffread -y "${HAP2_FAA}" -g "${HAP2_GENOME}" "${HAP2_GFF}"
HAP2_COUNT=$(grep -c '^>' "${HAP2_FAA}")
echo "[$(date)] Verification: Extracted ${HAP2_COUNT} nuclear protein sequences for Haplotype 2."

if [ "${HAP1_COUNT}" -eq 0 ] || [ "${HAP2_COUNT}" -eq 0 ]; then
    echo "ERROR: Protein extraction resulted in an empty FASTA file. Aborting pipeline." >&2
    exit 1
fi

# ==============================================================================
# 4. STEP 2: MULTI-LINEAGE BUSCO ANALYSIS
# ==============================================================================
echo -e "\n[$(date)] ==== STEP 2: RUNNING BUSCO COMPARISONS ===="

DB_KEYS=("eukaryota" "embryophyta" "eudicots" "malpighiales")
DB_PATHS=(
  "/path/to/your/directory/busco_downloads/lineages/eukaryota_odb12"
  "/path/to/your/directory/busco_downloads/lineages/embryophyta_odb10"
  "/path/to/your/directory/busco_downloads/lineages/eudicots_odb10"
  "/path/to/your/directory/busco_downloads/lineages/malpighiales_odb12"
)

RUN_NAMES=("Hap1_Genome" "Hap1_Protein" "Hap2_Genome" "Hap2_Protein")
RUN_FILES=("${HAP1_GENOME}" "${HAP1_FAA}" "${HAP2_GENOME}" "${HAP2_FAA}")
RUN_MODES=("genome" "protein" "genome" "protein")

for i in "${!DB_KEYS[@]}"; do
    DB_KEY="${DB_KEYS[$i]}"
    DB_PATH="${DB_PATHS[$i]}"
    
    echo "=========================================================="
    echo "[$(date)] Processing Database: ${DB_KEY}"
    echo "Lineage Path: ${DB_PATH}"
    echo "=========================================================="
    
    DB_OUT_DIR="${OUT_DIR}/${DB_KEY}"
    mkdir -p "${DB_OUT_DIR}"
    cd "${DB_OUT_DIR}"
    
    for j in "${!RUN_NAMES[@]}"; do
        R_NAME="${RUN_NAMES[$j]}"
        R_FILE="${RUN_FILES[$j]}"
        R_MODE="${RUN_MODES[$j]}"
        
        BUSCO_RUN_NAME="${R_NAME}_${DB_KEY}"
        
        echo "[$(date)] [RUN] ${BUSCO_RUN_NAME} | Mode: ${R_MODE}"
        
        busco \
          -i "${R_FILE}" \
          -o "${BUSCO_RUN_NAME}" \
          -l "${DB_PATH}" \
          -m "${R_MODE}" \
          -c 16 \
          --offline \
          --download_path "${DOWNLOAD_DIR}" \
          -f
    done
done

# ==============================================================================
# 5. STEP 3: RE-ACTIVATE PYTHON ENVIRONMENT & PLOT DATA
# ==============================================================================
echo -e "\n[$(date)] ==== STEP 3: STARTING PER-DATABASE PLOT GENERATION ===="

VENV_DIR="/path/to/your/directory/busco_plot_env"

if [ ! -d "${VENV_DIR}" ]; then
    python3 -m venv "${VENV_DIR}"
fi

source "${VENV_DIR}/bin/activate"
pip install --upgrade pip --quiet
pip install matplotlib numpy --quiet

cat << 'EOF' > "${OUT_DIR}/generate_final_plots.py"
import os
import sys
import glob
import re
import matplotlib.pyplot as plt
import numpy as np

def generate_db_plots(base_dir):
    db_keys = ["eukaryota", "embryophyta", "eudicots", "malpighiales"]
    run_types = ["Hap1_Genome", "Hap1_Protein", "Hap2_Genome", "Hap2_Protein"]
    
    c_single = '#34b3e4'
    c_dup = '#047495'
    c_frag = '#f7d838'
    c_miss = '#e93c26'

    for db in db_keys:
        print(f"Aggregating data for database: {db}")
        db_dir = os.path.join(base_dir, db)
        if not os.path.isdir(db_dir):
            print(f"Skipping database {db}: Folder not found.")
            continue
            
        data = {}
        for r_type in run_types:
            run_folder = f"{r_type}_{db}"
            search_path = os.path.join(db_dir, run_folder, "short_summary.specific.*.txt")
            summary_files = glob.glob(search_path)
            
            if not summary_files:
                print(f"  Warning: No summary files detected in {run_folder}")
                continue
                
            filepath = summary_files[0]
            
            with open(filepath, 'r') as f:
                content = f.read()
                
            try:
                s = int(re.search(r'\s+(\d+)\s+Complete and single-copy', content).group(1))
                d = int(re.search(r'\s+(\d+)\s+Complete and duplicated', content).group(1))
                frag = int(re.search(r'\s+(\d+)\s+Fragmented', content).group(1))
                miss = int(re.search(r'\s+(\d+)\s+Missing', content).group(1))
                total = int(re.search(r'\s+(\d+)\s+Total BUSCO groups searched', content).group(1))
                
                pretty_label = r_type.replace("_", " ")
                
                data[pretty_label] = {
                    'Single': (s / total) * 100,
                    'Duplicated': (d / total) * 100,
                    'Fragmented': (frag / total) * 100,
                    'Missing': (miss / total) * 100
                }
            except AttributeError as e:
                print(f"  Error: Failed to regex-parse BUSCO values from {filepath}: {e}")
                continue
        
        if not data:
            print(f"  Skipping plot generation for {db} - empty data.")
            continue
            
        labels = list(data.keys())
        single = [data[k]['Single'] for k in labels]
        dup = [data[k]['Duplicated'] for k in labels]
        frag = [data[k]['Fragmented'] for k in labels]
        miss = [data[k]['Missing'] for k in labels]
        
        fig, ax = plt.subplots(figsize=(11, 4.5))
        y_pos = np.arange(len(labels))
        
        ax.barh(y_pos, single, color=c_single, edgecolor='none', height=0.55, label='Complete (single-copy)')
        ax.barh(y_pos, dup, left=single, color=c_dup, edgecolor='none', height=0.55, label='Complete (duplicated)')
        ax.barh(y_pos, frag, left=np.array(single)+np.array(dup), color=c_frag, edgecolor='none', height=0.55, label='Fragmented')
        ax.barh(y_pos, miss, left=np.array(single)+np.array(dup)+np.array(frag), color=c_miss, edgecolor='none', height=0.55, label='Missing')
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=11, fontweight='bold')
        ax.set_xlabel('% BUSCOs', fontsize=12, fontweight='semibold')
        ax.set_title(f'BUSCO Comparison (Nuclear Chromosomes) - {db.upper()} Lineage', fontsize=13, fontweight='bold', pad=15)
        ax.set_xlim(0, 100)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#cccccc')
        ax.spines['bottom'].set_color('#cccccc')
        ax.xaxis.grid(True, linestyle='--', alpha=0.5, color='#cccccc')
        ax.set_axisbelow(True)
        
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.18), ncol=4, frameon=False, fontsize=9)
        
        plt.tight_layout()
        out_plot_path = os.path.join(base_dir, f"{db}_busco_comparison.png")
        plt.savefig(out_plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Successfully saved comparison plot to: {out_plot_path}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python generate_final_plots.py <base_dir>")
        sys.exit(1)
    generate_db_plots(sys.argv[1])
EOF

python3 "${OUT_DIR}/generate_final_plots.py" "${OUT_DIR}"

deactivate

echo -e "\n[$(date)] ==== ALL COMPLETED SUCCESSFULLY ===="