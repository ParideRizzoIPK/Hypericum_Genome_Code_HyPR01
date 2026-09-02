#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=48:00:00       
#SBATCH --job-name=pre_ragtag_eudicots
#SBATCH --output=logs/qc_%j.out
#SBATCH --error=logs/qc_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  

echo "#########################################################"
echo "Pipeline: Step 1 - Pre-RagTag QC (Eudicots Database & All Scaffolds)"
echo "Date: $(date)"
echo "#########################################################"

# --- 1. ENVIRONMENT & MODULE INITIALIZATION ---
source ~/.bashrc
set -u            

# Load the cluster-native modules
module load samtools/1.23.1
module load BUSCO/5.8.2

# Container Bind Fixes for Apptainer module architecture
export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# --- 2. DIRECTORIES & FILES ---
# New designated output folder
OUTPUT_DIR="/path/to/output/directory"
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

HAP1_FA="/path/to/input/directory/HyPR-01_v10_telo.hap1_scaffolds_final.fa"
HAP2_FA="/path/to/input/directory/HyPR-01_v10_telo.hap2_scaffolds_final.fa"

THREADS=${SLURM_CPUS_PER_TASK}

# CHANGE 2: Targeting the narrower eudicots database
LINEAGE="eudicots_odb10"

# ==============================================================================
# STEP 1.1: GENERATE SAMTOOLS INDEXES
# ==============================================================================
echo ">>>> Indexing Assemblies via Samtools faidx <<<<"
samtools faidx "${HAP1_FA}"
samtools faidx "${HAP2_FA}"

cp "${HAP1_FA}.fai" ./hap1.fai
cp "${HAP2_FA}.fai" ./hap2.fai

# ==============================================================================
# STEP 1.2: COMPUTE ASSEMBLY STATISTICS VIA INLINE PYTHON
# ==============================================================================
echo ">>>> Running Python Stats Parser <<<<"

cat << 'EOF' > parse_fai.py
import sys

fai_file = sys.argv[1]
out_file = sys.argv[2]
label = sys.argv[3]

lengths = []
with open(fai_file, 'r') as f:
    for line in f:
        if line.strip():
            parts = line.split('\t')
            lengths.append(int(parts[1]))

lengths.sort(reverse=True)
total_len = sum(lengths)
num_scaffs = len(lengths)

cum_sum = 0
n50 = None
l50 = None
n90 = None

for idx, l in enumerate(lengths):
    cum_sum += l
    if n50 is None and cum_sum >= total_len * 0.5:
        n50 = l
        l50 = idx + 1
    if n90 is None and cum_sum >= total_len * 0.9:
        n90 = l

with open(out_file, 'w') as out:
    out.write(f"========================================\n")
    out.write(f"Assembly Statistics Report: {label}\n")
    out.write(f"========================================\n")
    out.write(f"Total Assembly Length: {total_len:,} bp\n")
    out.write(f"Number of Scaffolds:   {num_scaffs}\n")
    out.write(f"N50 Length:            {n50:,} bp\n")
    out.write(f"L50 Count:             {l50}\n")
    out.write(f"N90 Length:            {n90:,} bp\n\n")
    
    # CHANGE 1: Printing sizes for ALL scaffolds instead of restricting to top 10
    out.write(f"Sizes of All Scaffolds (Ordered by Size):\n")
    for i, l in enumerate(lengths):
        out.write(f"  Scaffold {i+1:03d}: {l:,} bp\n")
EOF

python3 parse_fai.py hap1.fai hap1_contiguity_stats.txt "Haplotype 1 (Pre-RagTag Eudicots Run)"
python3 parse_fai.py hap2.fai hap2_contiguity_stats.txt "Haplotype 2 (Pre-RagTag Eudicots Run)"

rm -f parse_fai.py hap1.fai hap2.fai

# ==============================================================================
# STEP 1.3: RUN BUSCO ON BOTH HAPLOTYPES (EUDICOTS)
# ==============================================================================
echo ">>>> Starting BUSCO Eudicots Analysis for Haplotype 1 <<<<"
busco \
    -i "${HAP1_FA}" \
    -o busco_hap1_eudicots \
    -m genome \
    -l "${LINEAGE}" \
    -c ${THREADS} \
    --quiet

echo ">>>> Starting BUSCO Eudicots Analysis for Haplotype 2 <<<<"
busco \
    -i "${HAP2_FA}" \
    -o busco_hap2_eudicots \
    -m genome \
    -l "${LINEAGE}" \
    -c ${THREADS} \
    --quiet

echo "#########################################################"
echo "PIPELINE COMPLETE."
echo "Outputs generated in: ${OUTPUT_DIR}"
echo "#########################################################"