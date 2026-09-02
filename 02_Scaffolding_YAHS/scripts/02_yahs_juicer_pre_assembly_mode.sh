#!/bin/bash
 
# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G
#SBATCH --time=12:00:00
#SBATCH --job-name=yahs_jbat_fix_v3
#SBATCH --output=logs/jbat_fix_v3_%j.out
#SBATCH --error=logs/jbat_fix_v3_%j.err
 
# --- SCRIPT SETUP ---
set -e
set -o pipefail
 
# --- 1. DEFINE DIRECTORIES & PATHS ---
OUT_DIR="/path/to/output/directory"
JUICER_TOOLS_JAR="/path/to/your/directory/juicer_tools_1.22.01.jar"
cd "${OUT_DIR}"

# --- 2. LOAD MODULES ---
module load bwa-mem2/2.2.1
module load samtools/1.23.1
module load yahs/1.2
 
# --- 3. AUTO-DETECT THE JUICER CONVERSION TOOL ---
YAHS_BIN_DIR=$(dirname $(which yahs))
if [ -f "${YAHS_BIN_DIR}/juicer" ]; then
    JUICER_PRE_CMD="${YAHS_BIN_DIR}/juicer"
elif [ -f "${YAHS_BIN_DIR}/juicer_pre" ]; then
    JUICER_PRE_CMD="${YAHS_BIN_DIR}/juicer_pre"
else
    JUICER_PRE_CMD="/opt/Bio/yahs/1.2/bin/juicer"
fi

echo "Using Juicer conversion tool: ${JUICER_PRE_CMD}"
 
# --- FUNCTION TO PROCESS HAPLOTYPES ---
generate_jbat_files() {
    local HAP_NUM=$1
    local PREFIX="HyPR-01_v3.hap${HAP_NUM}"
    
    echo "========================================================"
    echo ">>>> GENERATING JBAT FILES FOR HAPLOTYPE ${HAP_NUM} <<<<"
    echo "========================================================"
    
    # Input Files
    local BIN_FILE="${PREFIX}.bin"
    local AGP_FILE="${PREFIX}_scaffolds_final.agp"
    local CTG_FA_INDEX="HyPR01_hap${HAP_NUM}.fa.fai"
    
    # Check inputs
    if [[ ! -f "${BIN_FILE}" || ! -f "${AGP_FILE}" || ! -f "${CTG_FA_INDEX}" ]]; then
        echo "ERROR: Missing required input files for HAP${HAP_NUM}."
        return
    fi

    local JBAT_PREFIX="${PREFIX}_JBAT"
    local TEMP_LOG="${JBAT_PREFIX}_pre_log.txt"
    
    # --- STEP 1: Run Juicer Pre (-a) and CAPTURE LOGS ---
    echo "[HAP${HAP_NUM}] 1. Running 'juicer pre -a'..."
    
    # We pipe stderr to a temp file so we can grep the size later
    if ${JUICER_PRE_CMD} pre -a -o "${JBAT_PREFIX}" "${BIN_FILE}" "${AGP_FILE}" "${CTG_FA_INDEX}" 2> "${TEMP_LOG}"; then
        echo "[HAP${HAP_NUM}] 'juicer pre' completed."
        cat "${TEMP_LOG}" # Print log to main stdout for SLURM record
    else
        echo "ERROR: 'juicer pre' failed."
        cat "${TEMP_LOG}"
        exit 1
    fi
    
    # Verify outputs
    local JBAT_TXT="${JBAT_PREFIX}.txt"
    if [ ! -f "${JBAT_TXT}" ]; then
        echo "ERROR: ${JBAT_TXT} was not created."
        exit 1
    fi

    # Handling liftover renaming
    if [ -f "${JBAT_PREFIX}.liftover" ]; then
        mv "${JBAT_PREFIX}.liftover" "${JBAT_PREFIX}.liftover.agp"
    fi

    # --- STEP 2: Extract EXACT Assembly Size ---
    # We grep the value from the log: "PRE_C_SIZE: assembly 123456789"
    echo "[HAP${HAP_NUM}] 2. Extracting exact assembly size from logs..."
    
    # grep the line, get the last column (size)
    local CORRECT_SIZE=$(grep "PRE_C_SIZE: assembly" "${TEMP_LOG}" | awk '{print $NF}')
    
    if [[ -z "${CORRECT_SIZE}" ]]; then
        echo "ERROR: Could not find 'PRE_C_SIZE: assembly' in logs. Cannot determine correct size."
        exit 1
    fi
    
    echo "[HAP${HAP_NUM}]    Detected Exact Size: ${CORRECT_SIZE}"
    
    local JBAT_CHROM_SIZES="${JBAT_PREFIX}.chrom.sizes"
    echo "assembly ${CORRECT_SIZE}" > "${JBAT_CHROM_SIZES}"

    # --- STEP 3: Generate the .hic file (Java) ---
    local JBAT_HIC="${JBAT_PREFIX}.hic"
    
    echo "[HAP${HAP_NUM}] 3. Building .hic matrix (Java)..."
    
    java -Xmx160G -jar ${JUICER_TOOLS_JAR} pre \
        -j ${SLURM_CPUS_PER_TASK} \
        "${JBAT_TXT}" \
        "${JBAT_HIC}" \
        "${JBAT_CHROM_SIZES}"
        
    if [ -f "${JBAT_HIC}" ]; then
        echo "[HAP${HAP_NUM}] SUCCESS: Created ${JBAT_HIC}"
        
        # Cleanup
        rm "${JBAT_TXT}"
        rm "${JBAT_CHROM_SIZES}"
        rm "${TEMP_LOG}"
    else
        echo "ERROR: .hic generation failed."
        exit 1
    fi
}
 
# --- EXECUTE ---
generate_jbat_files "1"
generate_jbat_files "2"
 
echo "========================================================"
echo "FULL JBAT GENERATION COMPLETE"
echo "========================================================"