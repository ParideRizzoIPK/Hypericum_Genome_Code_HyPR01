#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --partition=cpu        
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16     
#SBATCH --mem=35G              
#SBATCH --time=02:00:00        
#SBATCH --job-name=STAR_genome_indexer
#SBATCH --output=logs/star_indexer_%j.out
#SBATCH --error=logs/star_indexer_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  # Defer strict -u until after custom user profiles are loaded

# ==============================================================================
#                       ENVIRONMENT MODULE LOADING
# ==============================================================================
echo "[$(date)] Loading environment modules..."
source ~/.bashrc

# Safely enable strict variable evaluation now that user profile layers are loaded
set -u  

module purge 2>/dev/null || true
module load STAR/2.7.9a

export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# ==============================================================================
#                      DIAGNOSTICS & PATH VALIDATION
# ==============================================================================
echo "[$(date)] Verifying reference genome assets..."

OUTPUT_DIR="/path/to/output/directory"
INDEX_DIR="${OUTPUT_DIR}/STAR_combined_index"

HAP1_FA="/path/to/input/directory/hap1.masked_unchr.fa"
HAP2_FA="/path/to/input/directory/hap2.masked_unchr.fa"
CP_FA="/path/to/input/directory/HyPR01_Chloroplast_Genome.fasta"
MT_FA="/path/to/input/directory/HyPR01_Mitochondria_genome_collapsed.fasta"

GM_KEY_GZ="/path/to/your/directory/GeneMark_ETP_Software_Data/gm_key.gz"
GM_TAR_GZ="/path/to/your/directory/GeneMark_ETP_Software_Data/gmetp_linux_64.tar.gz"

for file in "$HAP1_FA" "$HAP2_FA" "$CP_FA" "$MT_FA" "$GM_KEY_GZ" "$GM_TAR_GZ"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required genomic asset missing: $file" >&2
        exit 1
    fi
done

mkdir -p "$INDEX_DIR"

# ==============================================================================
#           STEP 1: PRE-STAGING GENEMARK ENGINE AND LICENSES (FOR PHASE 3)
# ==============================================================================
echo "[$(date)] Unpacking GeneMark-ETP cluster software environments..."
SOFTWARE_OUT_DIR="${OUTPUT_DIR}/software"
mkdir -p "$SOFTWARE_OUT_DIR"

if [ ! -d "${SOFTWARE_OUT_DIR}/gmetp_linux_64" ]; then
    tar -xzf "$GM_TAR_GZ" -C "$SOFTWARE_OUT_DIR"
    echo "✔ GeneMark suite extracted locally."
fi

if [ ! -f ~/.gm_key ]; then
    zcat "$GM_KEY_GZ" > ~/.gm_key
    echo "✔ Global GeneMark runtime license key staged at ~/.gm_key."
fi

# ==============================================================================
#           STEP 2: FASTA SUFFIX INJECTION & RE-HEADER CONCATENATION
# ==============================================================================
COMBINED_GENOME="${OUTPUT_DIR}/HyPR01_diploid_competitor_ref.fa"

if [ ! -f "$COMBINED_GENOME" ]; then
    echo "[$(date)] Injecting suffix tags and compiling joint reference space..."
    
    python3 - "$HAP1_FA" "$HAP2_FA" "$CP_FA" "$MT_FA" "$COMBINED_GENOME" << 'EOF'
import sys

h1_in, h2_in, cp_in, mt_in, combined_out = sys.argv[1:6]

def stream_suffixed_records(input_path, output_handle, suffix):
    with open(input_path, 'r') as fin:
        for line in fin:
            if line.startswith('>'):
                parts = line.strip().split()
                suffixed_header = f"{parts[0]}_{suffix}"
                if len(parts) > 1:
                    suffixed_header += " " + " ".join(parts[1:])
                output_handle.write(suffixed_header + '\n')
            else:
                output_handle.write(line)

with open(combined_out, 'w') as fout:
    print("  -> Appending _Hap1 tags...")
    stream_suffixed_records(h1_in, fout, "Hap1")
    print("  -> Appending _Hap2 tags...")
    stream_suffixed_records(h2_in, fout, "Hap2")
    print("  -> Appending Organellar Sponge tracks...")
    stream_suffixed_records(cp_in, fout, "Chloroplast_Sponge")
    stream_suffixed_records(mt_in, fout, "Mitochondrion_Sponge")
EOF
    echo "✔ Unified competitive genome compiled successfully."
else
    echo "✔ Unified competitive reference genome already exists."
fi

# ==============================================================================
#               STEP 3: CONSTRUCT UNIFIED STAR GENOME INDEX
# ==============================================================================
THREADS=${SLURM_CPUS_PER_TASK}

echo "[$(date)] Commencing unified STAR genome index generation..."
STAR --runThreadN "$THREADS" \
     --runMode genomeGenerate \
     --genomeDir "$INDEX_DIR" \
     --genomeFastaFiles "$COMBINED_GENOME" \
     --genomeSAindexNbases 13

echo "[$(date)] Genome indexing complete. You can now submit Script 2."