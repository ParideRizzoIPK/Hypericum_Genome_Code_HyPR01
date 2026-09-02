#!/bin/bash
 
# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G
#SBATCH --time=96:00:00
#SBATCH --job-name=yahs_telo_v22
#SBATCH --output=logs/yahs_integrated_%j.out
#SBATCH --error=logs/yahs_integrated_%j.err
 
# --- DEBUG START ---
echo "DEBUG: Script started at $(date)"
 
# --- SCRIPT SETUP ---
set -e
set -o pipefail
 
# --- 1. DIRECTORIES ---
OUT_DIR="/path/to/output/directory"
mkdir -p "${OUT_DIR}"
 
# Tool Paths
YAHS_ROOT="/path/to/your/directory"
YAHS_BIN="${YAHS_ROOT}/yahs"
JUICER_PRE_CMD="${YAHS_ROOT}/juicer"
YAHS_SCRIPTS="${YAHS_ROOT}/scripts"
JUICER_TOOLS_JAR="/path/to/your/directory/juicer_tools_1.22.01.jar"

# --- 2. PERMISSIONS & PATH ---
export PATH="${YAHS_SCRIPTS}:${YAHS_ROOT}:${PATH}"

chmod +x "${YAHS_BIN}" "${JUICER_PRE_CMD}"
if [ -f "${YAHS_SCRIPTS}/statistics" ]; then chmod +x "${YAHS_SCRIPTS}/statistics"; fi

# --- 3. LOAD MODULES ---
# FIX: Updated to the correct version available on your cluster
module load bwa-mem2/2.2.1
module load samtools/1.23.1
module load openjdk/11 || echo "Warning: module openjdk/11 not found, using system Java"

# --- 4. VARIABLES ---
TOTAL_CPUS=${SLURM_CPUS_PER_TASK}
BWA_THREADS=40
SAM_THREADS=8
Q_THRESHOLD=20
TELO_MOTIF="CCCTAAA" 
 
# Input Files
SRC_HAP1="/path/to/input/directory/HyPR01_genome_hic_assembly.hic.hap1.p_ctg.fa"
SRC_HAP2="/path/to/input/directory/HyPR01_genome_hic_assembly.hic.hap2.p_ctg.fa"
HIC_READ1="/path/to/input/directory/2037511_Hypericum_PR1_HiC166_S68_L002_R1_001_val_1.fq.gz"
HIC_READ2="/path/to/input/directory/2037511_Hypericum_PR1_HiC166_S68_L002_R2_001_val_2.fq.gz"
 
# --- FUNCTION ---
process_haplotype() {
    local HAP_NUM=$1
    local SRC_FA=$2
    local LOCAL_FA="${OUT_DIR}/HyPR01_hap${HAP_NUM}.fa"
    local PREFIX="${OUT_DIR}/HyPR-01_v10_telo.hap${HAP_NUM}"
    
    echo "========================================================"
    echo ">>>> PROCESSING HAPLOTYPE ${HAP_NUM} <<<<"
    echo "========================================================"
    cd "${OUT_DIR}"
 
    # 5.1 Indexing
    if [ ! -f "${LOCAL_FA}" ]; then cp "${SRC_FA}" "${LOCAL_FA}"; fi
    if [ ! -f "${LOCAL_FA}.fai" ]; then samtools faidx "${LOCAL_FA}"; fi
    if [ ! -f "${LOCAL_FA}.bwt.2bit.64" ]; then bwa-mem2.avx2 index "${LOCAL_FA}"; fi

    # 5.2 Alignment (With Integrity Check)
    BAM_FILE="${PREFIX}.hic_aligned.bam"
    NEED_ALIGNMENT=true

    if [ -s "${BAM_FILE}" ]; then
        echo "[HAP${HAP_NUM}] Checking existing BAM integrity..."
        if samtools quickcheck "${BAM_FILE}"; then
            echo "[HAP${HAP_NUM}] BAM is valid. Skipping alignment."
            NEED_ALIGNMENT=false
        else
            echo "[HAP${HAP_NUM}] BAM is corrupt. Re-aligning..."
            rm -f "${BAM_FILE}"
        fi
    fi

    if [ "$NEED_ALIGNMENT" = true ]; then
        echo "[HAP${HAP_NUM}] Aligning Hi-C reads..."
        bwa-mem2.avx2 mem -5SP -t ${BWA_THREADS} "${LOCAL_FA}" "${HIC_READ1}" "${HIC_READ2}" | \
        samtools view -Sb -@ ${SAM_THREADS} - > "${BAM_FILE}"
    fi
 
    # 5.3 YAHS Scaffolding
    # Cleanup previous partial run outputs
    rm -f "${PREFIX}.bin" "${PREFIX}_scaffolds_final.agp" "${PREFIX}_scaffolds_final.fa"

    echo "[HAP${HAP_NUM}] Running YAHS with --telo-motif ${TELO_MOTIF}..."
    "${YAHS_BIN}" -q ${Q_THRESHOLD} --telo-motif ${TELO_MOTIF} -o "${PREFIX}" "${LOCAL_FA}" "${BAM_FILE}"
 
    # Verification
    if [ ! -f "${PREFIX}.bin" ]; then
        echo "ERROR: YAHS failed to create output files."
        exit 1
    fi

    # 5.4 Generate JBAT Files
    echo "[HAP${HAP_NUM}] Running Juicer Pre (Assembly Mode)..."
    "${JUICER_PRE_CMD}" pre -a -o "${PREFIX}_JBAT" \
        "${PREFIX}.bin" \
        "${PREFIX}_scaffolds_final.agp" \
        "${LOCAL_FA}.fai"
    
    [ -f "${PREFIX}_JBAT.liftover" ] && mv "${PREFIX}_JBAT.liftover" "${PREFIX}_JBAT.liftover.agp"

    # 5.5 Java Matrix
    echo "[HAP${HAP_NUM}] Indexing final scaffolds to calculate size..."
    
    # Generate the index for the NEW scaffolds (Fixed in v20)
    samtools faidx "${PREFIX}_scaffolds_final.fa"

    echo "[HAP${HAP_NUM}] Building .hic matrix (Java)..."
    local TOTAL_SIZE=$(awk '{sum+=$2} END {print sum}' "${PREFIX}_scaffolds_final.fa.fai")
    
    echo "assembly ${TOTAL_SIZE}" > "${PREFIX}_JBAT.chrom.sizes"

    java -Xmx160G -jar ${JUICER_TOOLS_JAR} pre \
        -j ${SLURM_CPUS_PER_TASK} \
        "${PREFIX}_JBAT.txt" \
        "${PREFIX}_JBAT.hic" \
        "${PREFIX}_JBAT.chrom.sizes"
    
    echo "[HAP${HAP_NUM}] COMPLETE: ${PREFIX}_JBAT.hic created."
}
 
# --- EXECUTE ---
process_haplotype "1" "${SRC_HAP1}"
process_haplotype "2" "${SRC_HAP2}"
 
echo "========================================================"
echo "ALL PROCESSES COMPLETE!"
echo "========================================================"