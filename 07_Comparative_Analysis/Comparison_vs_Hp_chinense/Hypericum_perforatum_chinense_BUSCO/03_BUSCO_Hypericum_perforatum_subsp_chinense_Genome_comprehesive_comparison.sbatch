#!/bin/bash
#
#SBATCH --job-name=Hp_Subgenomes_BUSCO
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=logs/busco_subgenomes_%j.out
#SBATCH --error=logs/busco_subgenomes_%j.err

# Strict error handling
set -euo pipefail

# ==============================================================================
# 1. ENVIRONMENT & PATH SETUP
# ==============================================================================
export SINGULARITY_BIND="/path/to/your/directory"

# Load modules from cluster
module load BUSCO/5.8.2
module load python/3.12.0

# Master paths
BASE_DIR="/path/to/your/directory"
SPLIT_DIR="${BASE_DIR}/Hypericum_perforatum_subsp_chinense/HyPR01_synteny_v2/split"

# New Target Output Directory (Updated per your request)
OUT_DIR="${BASE_DIR}/Hypericum_perforatum_subsp_chinense/Hypericum_perforatum_chinense_BUSCO"
CONCAT_DIR="${OUT_DIR}/concatenated_subgenomes"
DOWNLOAD_DIR="/path/to/your/directory/busco_downloads"

# Initialize workspaces dynamically prior to logging or concatenation
mkdir -p "${OUT_DIR}"
mkdir -p "${CONCAT_DIR}"

# Input verification guard rail
echo "[$(date)] ==== STEP 1: VERIFYING SPLIT CHROMOSOMES ===="
for S in A B C D; do
    for CHR in {1..8}; do
        FILE="${SPLIT_DIR}/ref_Chr${CHR}${S}.fa"
        if [ ! -f "$FILE" ]; then
            echo "ERROR: Missing split file required for concatenation: $FILE" >&2
            exit 1
        fi
    done
done
echo "All 32 split chromosome files detected and successfully validated."

# ==============================================================================
# 2. STEP 2: RECONSTRUCTING SUBGENOMES (CONCATENATION)
# ==============================================================================
echo -e "\n[$(date)] ==== STEP 2: CONCATENATING SUBGENOMES A, B, C, AND D ===="

for S in A B C D; do
    TARGET_FA="${CONCAT_DIR}/subgenome_${S}.fa"
    echo "[$(date)] Merging chromosomes 1 to 8 into subgenome_${S}.fa..."
    
    # Wipe target file if it already exists to avoid double-appending
    > "${TARGET_FA}"
    
    # Append chromosome FASTA sequences chronologically
    for CHR in {1..8}; do
        cat "${SPLIT_DIR}/ref_Chr${CHR}${S}.fa" >> "${TARGET_FA}"
    done
    
    # Validate final file size and header structure
    HEADER_COUNT=$(grep -c '^>' "${TARGET_FA}")
    echo "  -> Created ${TARGET_FA} (Contains ${HEADER_COUNT} sequences)."
done

# ==============================================================================
# 3. STEP 3: RUNNING MULTI-LINEAGE GENOME MODE BUSCO
# ==============================================================================
echo -e "\n[$(date)] ==== STEP 3: RUNNING BUSCO COMPARISONS ===="

# Define target lineage paths
DB_KEYS=("eukaryota" "embryophyta" "eudicots" "malpighiales")
DB_PATHS=(
  "${DOWNLOAD_DIR}/lineages/eukaryota_odb12"
  "${DOWNLOAD_DIR}/lineages/embryophyta_odb10"
  "${DOWNLOAD_DIR}/lineages/eudicots_odb10"
  "${DOWNLOAD_DIR}/lineages/malpighiales_odb12"
)

# Run definitions targeting your newly concatenated files
RUN_NAMES=("Subgenome_A" "Subgenome_B" "Subgenome_C" "Subgenome_D")
RUN_FILES=(
  "${CONCAT_DIR}/subgenome_A.fa"
  "${CONCAT_DIR}/subgenome_B.fa"
  "${CONCAT_DIR}/subgenome_C.fa"
  "${CONCAT_DIR}/subgenome_D.fa"
)

# Nested execution over each database and subgenome FASTA
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
        
        BUSCO_RUN_NAME="${R_NAME}_${DB_KEY}"
        
        echo "[$(date)] [RUNNING] ${BUSCO_RUN_NAME} | Mode: genome"
        
        busco \
          -i "${R_FILE}" \
          -o "${BUSCO_RUN_NAME}" \
          -l "${DB_PATH}" \
          -m genome \
          -c 16 \
          --offline \
          --download_path "${DOWNLOAD_DIR}" \
          -f
    done
done

# ==============================================================================
# 4. STEP 4: PLOTTING FINAL COMPARATIVE CHARTS
# ==============================================================================
echo -e "\n[$(date)] ==== STEP 4: STARTING VISUALIZATION PASS ===="

VENV_DIR="/path/to/your/directory/busco_plot_env"
if [ ! -d "${VENV_DIR}" ]; then
    python3 -m venv "${VENV_DIR}"
fi

source "${VENV_DIR}/bin/activate"
pip install --upgrade pip --quiet
pip install matplotlib numpy --quiet

# Dynamically generate the updated Python plotting engine
cat << 'EOF' > "${OUT_DIR}/generate_final_plots.py"
import os
import sys
import glob
import re
import matplotlib.pyplot as plt
import numpy as np

def generate_db_plots(base_dir):
    db_keys = ["eukaryota", "embryophyta", "eudicots", "malpighiales"]
    run_types = ["Subgenome_A", "Subgenome_B", "Subgenome_C", "Subgenome_D"]
    
    # Publication-quality color scheme
    c_single = '#34b3e4'  # Light Blue
    c_dup = '#047495'     # Deep Slate Blue
    c_frag = '#f7d838'    # Warm Yellow
    c_miss = '#e93c26'    # Crimson Red

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
                
                # Format label (e.g., Subgenome_A -> Subgenome A)
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
            
        # Parse datasets
        labels = list(data.keys())
        single = [data[k]['Single'] for k in labels]
        dup = [data[k]['Duplicated'] for k in labels]
        frag = [data[k]['Fragmented'] for k in labels]
        miss = [data[k]['Missing'] for k in labels]
        
        fig, ax = plt.subplots(figsize=(11, 4.5))
        y_pos = np.arange(len(labels))
        
        # Plot bars
        ax.barh(y_pos, single, color=c_single, edgecolor='none', height=0.55, label='Complete (single-copy)')
        ax.barh(y_pos, dup, left=single, color=c_dup, edgecolor='none', height=0.55, label='Complete (duplicated)')
        ax.barh(y_pos, frag, left=np.array(single)+np.array(dup), color=c_frag, edgecolor='none', height=0.55, label='Fragmented')
        ax.barh(y_pos, miss, left=np.array(single)+np.array(dup)+np.array(frag), color=c_miss, edgecolor='none', height=0.55, label='Missing')
        
        # Plot Styling
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=11, fontweight='bold')
        ax.set_xlabel('% BUSCOs', fontsize=12, fontweight='semibold')
        ax.set_title(f'Tetraploid Subgenome BUSCO Comparison - {db.upper()} Lineage', fontsize=13, fontweight='bold', pad=15)
        ax.set_xlim(0, 100)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#cccccc')
        ax.spines['bottom'].set_color('#cccccc')
        ax.xaxis.grid(True, linestyle='--', alpha=0.5, color='#cccccc')
        ax.set_axisbelow(True)
        
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.18), ncol=4, frameon=False, fontsize=9)
        
        plt.tight_layout()
        out_plot_path = os.path.join(base_dir, f"{db}_subgenome_busco_comparison.png")
        plt.savefig(out_plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Successfully saved comparison plot to: {out_plot_path}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python generate_final_plots.py <base_dir>")
        sys.exit(1)
    generate_db_plots(sys.argv[1])
EOF

# Run the plotting wrapper script
python3 "${OUT_DIR}/generate_final_plots.py" "${OUT_DIR}"

# Terminate Python environment safely
deactivate

echo -e "\n[$(date)] ==== SUBGENOME PIPELINE EXECUTED SUCCESSFULLY ===="