#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G            
#SBATCH --time=24:00:00
#SBATCH --job-name=hic_verification_check
#SBATCH --output=logs/hic_chk_%j.out
#SBATCH --error=logs/hic_chk_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  

echo "#########################################################"
echo "Pipeline: Step 4A - Hi-C Contact Map Validation (Path Typo Fixed)"
echo "Date: $(date)"
echo "#########################################################"

# --- 1. ENVIRONMENT & MODULE INITIALIZATION ---
source ~/.bashrc
set -u            

# Load required cluster modules
module load bwa-mem2/2.2.1
module load pairtools/1.0.3
module load samtools/1.23.1

# Container Bind Safety variables
export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# --- 2. CONFIGURATION AND DEFINITIONS ---
OUTPUT_DIR="/path/to/output/directory"
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

# Cluster pointer to Juicer Tools Java archive
JUICER_TOOLS_JAR="/path/to/your/directory/juicer_tools_1.22.01.jar"

# --- PATH FIX: Corrected directory folder from HiC_RAW to HyPR01_RAW ---
HIC_R1="/path/to/input/directory/2037511_Hypericum_PR1_HiC166_S68_L002_R1_001.fastq.gz"
HIC_R2="/path/to/input/directory/2037511_Hypericum_PR1_HiC166_S68_L002_R2_001.fastq.gz"

# Input Curated Assemblies from Step 4
HAP1_FA="/path/to/input/directory/Hap1_curated_chr_un.fasta"
HAP2_FA="/path/to/input/directory/Hap2_curated_chr_un_RENAMED.fasta"

THREADS=${SLURM_CPUS_PER_TASK}

# ==============================================================================
# PIPELINE RUN 1: HAPLOTYPE 1 VALIDATION MATRIX
# ==============================================================================
echo ">>>> STARTING VALIDATION FOR HAPLOTYPE 1 <<<<"
mkdir -p "${OUTPUT_DIR}/Hap1_Workspace"
cd "${OUTPUT_DIR}/Hap1_Workspace"

echo "  -> Indexing Haplotype 1 Assembly..."
ln -sf "${HAP1_FA}" ./Hap1_ref.fasta
samtools faidx Hap1_ref.fasta
cut -f1,2 Hap1_ref.fasta.fai > hap1.chrom.sizes
bwa-mem2 index Hap1_ref.fasta

echo "  -> Aligning Hi-C reads and streaming into Pairtools..."
bwa-mem2 mem -t "${THREADS}" -5SP Hap1_ref.fasta "${HIC_R1}" "${HIC_R2}" | \
    pairtools parse --chroms-path hap1.chrom.sizes | \
    pairtools sort --nproc "${THREADS}" --memory 64G -o hap1_sorted.pairs

echo "  -> Compiling final binary .hic matrix via Juicer Tools..."
java -Xmx160G -jar "${JUICER_TOOLS_JAR}" pre \
    -j "${THREADS}" \
    hap1_sorted.pairs \
    "${OUTPUT_DIR}/Hap1_verification.hic" \
    hap1.chrom.sizes

echo "  -> Cleaning up Haplotype 1 intermediate workspace scratch files..."
rm -rf ./Hap1_ref.fasta* hap1_sorted.pairs hap1.chrom.sizes


# ==============================================================================
# PIPELINE RUN 2: HAPLOTYPE 2 VALIDATION MATRIX
# ==============================================================================
echo ">>>> STARTING VALIDATION FOR HAPLOTYPE 2 <<<<"
mkdir -p "${OUTPUT_DIR}/Hap2_Workspace"
cd "${OUTPUT_DIR}/Hap2_Workspace"

echo "  -> Indexing Haplotype 2 Assembly..."
ln -sf "${HAP2_FA}" ./Hap2_ref.fasta
samtools faidx Hap2_ref.fasta
cut -f1,2 Hap2_ref.fasta.fai > hap2.chrom.sizes
bwa-mem2 index Hap2_ref.fasta

echo "  -> Aligning Hi-C reads and streaming into Pairtools..."
bwa-mem2 mem -t "${THREADS}" -5SP Hap2_ref.fasta "${HIC_R1}" "${HIC_R2}" | \
    pairtools parse --chroms-path hap2.chrom.sizes | \
    pairtools sort --nproc "${THREADS}" --memory 64G -o hap2_sorted.pairs

echo "  -> Compiling final binary .hic matrix via Juicer Tools..."
java -Xmx160G -jar "${JUICER_TOOLS_JAR}" pre \
    -j "${THREADS}" \
    hap2_sorted.pairs \
    "${OUTPUT_DIR}/Hap2_verification.hic" \
    hap2.chrom.sizes

echo "  -> Cleaning up Haplotype 2 intermediate workspace scratch files..."
rm -rf ./Hap2_ref.fasta* hap2_sorted.pairs hap2.chrom.sizes


# ==============================================================================
# PIPELINE WRAPUP
# ==============================================================================
cd "${OUTPUT_DIR}"
rm -rf "${OUTPUT_DIR}/Hap1_Workspace" "${OUTPUT_DIR}/Hap2_Workspace"

echo "#########################################################"
echo "HI-C VERIFICATION RUN SUCCESFULLY ENDED."
echo "Final Juicebox Matrix Targets:"
echo "  - Haplotype 1 Matrix: ${OUTPUT_DIR}/Hap1_verification.hic"
echo "  - Haplotype 2 Matrix: ${OUTPUT_DIR}/Hap2_verification.hic"
echo "#########################################################"