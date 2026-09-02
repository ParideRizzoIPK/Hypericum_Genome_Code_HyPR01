#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=24:00:00       
#SBATCH --job-name=ragtag_bidirectional
#SBATCH --output=logs/ragtag_%j.out
#SBATCH --error=logs/ragtag_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  

echo "#########################################################"
echo "Pipeline: Steps 2.2 & 2.3 - Labeled Bidirectional RagTag Curation"
echo "Date: $(date)"
echo "#########################################################"

# ==============================================================================
# USER CONFIGURATION CONTROLS (SANDBOX TOGGLE)
# ==============================================================================
# Leave as "false" for your primary pass. Toggle to "true" only if downstream 
# synteny plots suggest local structural changes are needed for Hap2.
RUN_OPTIONAL_HAP2_CURATION="false"

# --- 1. ENVIRONMENT & MODULE INITIALIZATION ---
source ~/.bashrc
set -u            

# Load the cluster-native software tools 
module load minimap2/2.24
module load RagTag/2.1.0

# --- CRITICAL CONTAINER BIND FIX ---
# Forces the Apptainer/Singularity runtime to mount your network storage
export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# --- 2. DIRECTORIES & INPUT FILE PATHS ---
OUTPUT_DIR="/path/to/output/directory"
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

# References and Query Assemblies
HAP2_REFF_ORIG="/path/to/input/directory/HyPR-01_v10_telo.hap2_scaffolds_final.fa"
HAP1_QRY_100KB="${OUTPUT_DIR}/Hap1_YAHS_min100kb.fa"
HAP2_QRY_100KB="${OUTPUT_DIR}/Hap2_YAHS_min100kb.fa"

THREADS=${SLURM_CPUS_PER_TASK}

# ==============================================================================
# STEP 2.2: PRIMARY PASS (Hap2 as Reference, Hap1 as Query)
# ==============================================================================
echo ">>>> STARTING STEP 2.2: PRIMARY PASS (Hap1 Curation) <<<<"

# --- Script B2: RagTag Correct ---
echo "  -> Running RagTag Correct on Haplotype 1..."
ragtag.py correct \
    --aligner minimap2 \
    -u \
    -t "${THREADS}" \
    -o ragtag_hap1_correct \
    "${HAP2_REFF_ORIG}" \
    "${HAP1_QRY_100KB}"

echo "  -> Renaming Hap1 Correct Fasta to explicitly include Haplotype ID..."
mv ragtag_hap1_correct/ragtag.correct.fasta ragtag_hap1_correct/HyPR-01_Hap1_ragtag_correct.fasta


# --- Script B3: RagTag Scaffold ---
echo "  -> Running RagTag Scaffold on Haplotype 1..."
ragtag.py scaffold \
    --aligner minimap2 \
    -C \
    -t "${THREADS}" \
    -o ragtag_hap1_scaffold \
    "${HAP2_REFF_ORIG}" \
    "ragtag_hap1_correct/HyPR-01_Hap1_ragtag_correct.fasta"

echo "  -> Renaming Hap1 Scaffold Fasta to explicitly include Haplotype ID..."
mv ragtag_hap1_scaffold/ragtag.scaffold.fasta ragtag_hap1_scaffold/HyPR-01_Hap1_ragtag_scaffold.fasta

echo ">>>> STEP 2.2 COMPLETE. Curated, Labeled Hap1 generated successfully. <<<<"
echo "---------------------------------------------------------"

# ==============================================================================
# STEP 2.3: OPTIONAL REVERSE PASS (Curated Hap1 as Reference, Hap2 as Query)
# ==============================================================================
if [ "${RUN_OPTIONAL_HAP2_CURATION}" = "true" ]; then
    echo ">>>> STARTING STEP 2.3: OPTIONAL REVERSE PASS (Hap2 Curation) <<<<"
    
    # Updated reference pointer to track the new explicitly labeled filename
    CURATED_HAP1_REF="ragtag_hap1_scaffold/HyPR-01_Hap1_ragtag_scaffold.fasta"
    
    # Run Correction on Hap2 guided by Curated Hap1
    echo "  -> Running Strict RagTag Correct on Haplotype 2..."
    ragtag.py correct \
        --aligner minimap2 \
        -u \
        -t "${THREADS}" \
        -o ragtag_hap2_correct \
        "${CURATED_HAP1_REF}" \
        "${HAP2_QRY_100KB}"
        
    echo "  -> Renaming Hap2 Correct Fasta to explicitly include Haplotype ID..."
    mv ragtag_hap2_correct/ragtag.correct.fasta ragtag_hap2_correct/HyPR-01_Hap2_ragtag_correct.fasta
        
    # Run Scaffolding on Hap2 using highly conservative thresholds
    echo "  -> Running Conservative Strict RagTag Scaffold on Haplotype 2..."
    ragtag.py scaffold \
        --aligner minimap2 \
        -C \
        -i 0.90 \
        -s 0.90 \
        --remove-small \
        -f 50000 \
        -t "${THREADS}" \
        -o ragtag_hap2_scaffold \
        "${CURATED_HAP1_REF}" \
        "ragtag_hap2_correct/HyPR-01_Hap2_ragtag_correct.fasta"
        
    echo "  -> Renaming Hap2 Scaffold Fasta to explicitly include Haplotype ID..."
    mv ragtag_hap2_scaffold/ragtag.scaffold.fasta ragtag_hap2_scaffold/HyPR-01_Hap2_ragtag_scaffold.fasta
        
    echo ">>>> STEP 2.3 COMPLETE. Curated, Labeled Hap2 generated successfully. <<<<"
else
    echo ">>>> STEP 2.3 NOTICE: Optional Reverse Pass skipped as requested. <<<<"
    echo "     Haplotype 2 remains your original high-quality YAHS-telo structure."
fi

echo "#########################################################"
echo "PIPELINE JOB RUN ENDED."
echo "Outputs stored in: ${OUTPUT_DIR}"
echo "#########################################################"