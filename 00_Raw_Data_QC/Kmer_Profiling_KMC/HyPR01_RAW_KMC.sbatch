#!/bin/bash
#SBATCH --job-name=KMC_genome_estimation
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8 # Adjust based on your KMC version and needs
#SBATCH --mem=512G          # Adjust based on your genome size
#SBATCH --time=2-00:00:00    # 2 days, adjust as needed
#SBATCH --output=logs/HyPR01RAW_KMC.out
#SBATCH --error=logs/HyPR01RAW_KMC.err
#SBATCH --partition=cpu # Specify your Slurm partition

# --- Load necessary modules ---
module load kmc # Or the specific module name for KMC on your cluster

# --- Define variables ---
FASTQ_FILE="/path/to/input/directory/hypr01_hifi_reads.fastq" # Your input PacBio CCS fastq file
KMER_SIZE=21                     # A common k-mer size for estimation
THREADS=8
OUTPUT_PREFIX="HyPR01_kmc"

# --- Create a temporary directory for KMC database ---
TMP_DIR=$(mktemp -d)

# --- Run KMC ---
echo "Starting KMC..."
kmc -k${KMER_SIZE} -t${THREADS} -m32 -ci1 -cs10000 @<(echo ${FASTQ_FILE}) ${OUTPUT_PREFIX} ${TMP_DIR}

# --- Generate k-mer histogram ---
kmc_tools transform ${OUTPUT_PREFIX} histogram ${OUTPUT_PREFIX}.hist -cx10000

echo "KMC finished. Histogram in ${OUTPUT_PREFIX}.hist"
echo "You can now use a tool like GenomeScope to estimate genome size from the .hist file."

# --- Clean up temporary directory ---
rm -r ${TMP_DIR}