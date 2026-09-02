#!/bin/bash

# ==============================================================================
# SBATCH Directives for SLURM Cluster
# ==============================================================================
#SBATCH --job-name=HiMT_HyPR01_organelle
#SBATCH --auks=yes
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --output=logs/himt_organelle_%j.out
#SBATCH --error=logs/himt_organelle_%j.err
#SBATCH --time=1-00:00:00

# ==============================================================================
# Strict Bash Error Handling
# ==============================================================================
set -e          # Exit immediately if a command exits with a non-zero status
set -u          # Treat unset variables as an error
set -o pipefail # Return exit status of the last command in a pipeline to fail

# ==============================================================================
# Directory and Input Path Definitions
# ==============================================================================
BASE_OUTPUT_DIR="/path/to/output/directory"
LOG_DIR="${BASE_OUTPUT_DIR}/slurm_logs"

# Input raw PacBio HiFi reads (full CCS read set, not downsampled -- per project
# decision to remove "why only a subsample" as a reviewer question)
RAW_READS="/path/to/input/directory/hypr01_hifi_reads.fastq"

# Published H. perforatum mitochondrial reference genome, for direct
# collinearity comparison via `himt compare` (future_workflow_plan.md, Stage 6).
# Downloaded automatically below (Step: "Reference genome acquisition").
#
# ------------------------------------------------------------------------
# REFERENCE METADATA (verified by direct download + header/length check,
# not taken from the paper's text alone -- see notes below):
#
#   Source database : GenBase, National Genomics Data Center (NGDC),
#                      China National Center for Bioinformation (CNCB)
#                      https://ngdc.cncb.ac.cn/genbase/
#   Retrieval API    : https://ngdc.cncb.ac.cn/genbase/api/file/fasta?acc=<ACCESSION>
#   Associated paper : Hypericum perforatum mitochondrial genome study,
#                      Frontiers in Plant Science (2026),
#                      https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2026.1763082/full
#   Raw reads (SRA/GSA equiv.) : CRA035721 (Genome Sequence Archive, NGDC) -- not downloaded here
#
#   Accession   Chromosome  GenBase label              Verified length
#   ----------  ----------  -------------------------  ---------------
#   C_AA133241.1  Chr1      "partial genome"*          204,492 bp
#   C_AA133242.1  Chr2      "complete genome"          280,636 bp
#                                              Combined total: 485,128 bp
#
#   * GenBase's "partial genome" label on Chr1 is a labeling convention for
#     a companion/second chromosome record, not truncated sequence -- its
#     length matches the paper's reported Chr1 size exactly, confirmed here
#     by direct download rather than assumed from the label text.
#
#   IMPORTANT -- accession specificity: the paper's Data Availability
#   statement lists TWO accession ranges covering TWO different Hypericum
#   species (H. perforatum and H. ascyron). C_AA133240.1, adjacent to the
#   correct accessions, is H. ascyron mitochondrion, NOT H. perforatum --
#   confirmed by downloading and checking its header directly. Do not
#   substitute nearby accession numbers without re-verifying the header.
# ------------------------------------------------------------------------
REFERENCE_MITO_ACCESSIONS=("C_AA133241.1" "C_AA133242.1")
REFERENCE_DIR="${BASE_OUTPUT_DIR}/reference"
REFERENCE_MITO_GENOME="${REFERENCE_DIR}/Hperforatum_mitochondrion_reference.fasta"

# Bind paths for Singularity/Apptainer subroutines
export SINGULARITY_BIND="/path/to/your/directory,/path/to/your/other/directory"

# ==============================================================================
# Pre-Run Setup & Diagnostics
# ==============================================================================
echo "======================================================================"
echo "INFO: Starting HiMT Organelle Genome Assembly (chloroplast + mitochondrion)"
echo "Date / Time      : $(date)"
echo "Execution Host   : $(hostname)"
echo "Allocated CPUs   : ${SLURM_CPUS_PER_TASK}"
echo "Allocated Memory : 32G"
echo "Target Directory : ${BASE_OUTPUT_DIR}"
echo "Log Directory    : ${LOG_DIR}"
echo "Input HiFi Reads : ${RAW_READS}"
echo "======================================================================"

# Ensure directories exist
mkdir -p "${BASE_OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"

# Verify input read file presence
if [ ! -s "${RAW_READS}" ]; then
    echo "ERROR: Input FASTQ file not found or is empty: ${RAW_READS}" >&2
    exit 1
fi

# Change to target directory
cd "${BASE_OUTPUT_DIR}"

# ==============================================================================
# Direct Environment Injection
# ==============================================================================
echo "INFO: Injecting HiMT environment directly into PATH..."

# Directly prepend the environment bin to the system PATH to avoid hook errors
export PATH="/path/to/your/directory/micromamba/envs/himt/bin:${PATH}"

# Validate tool presence. These are diagnostics only (non-fatal, `|| true`) --
# the real test is the `himt assemble` call itself below. HiMT's documented
# dependencies are python>=3.10, blast>=2.14.0, flye>=2.9.4, miniprot>=0.13,
# minimap2>=2.24-r1122, plotly, numpy -- checking all the external binaries,
# not just flye/blastn, so a missing one is caught here rather than mid-run.
echo "INFO: Checking dependency tool versions (non-fatal -- for diagnostics only):"
himt -h    > /dev/null 2>&1 && echo "  himt: OK"     || echo "  WARNING: 'himt -h' failed -- check the himt conda env is active"
flye --version                 || echo "  WARNING: flye executable not found in PATH"
blastn -version | head -n 1    || echo "  WARNING: blastn executable not found in PATH"
miniprot --version             || echo "  WARNING: miniprot executable not found in PATH"
minimap2 --version             || echo "  WARNING: minimap2 executable not found in PATH"

# ==============================================================================
# Run Assembly
# ==============================================================================
echo "======================================================================"
echo "INFO: Executing 'himt assemble'..."
echo "INFO: No organelle-selection flag is passed -- HiMT assembles both the"
echo "      chloroplast and mitochondrial genomes in a single run by default."
echo "      (An earlier draft of this script included '--organelle both',"
echo "      which is not a real HiMT option and would have caused an"
echo "      argument-parsing failure; removed.)"
echo "======================================================================"

himt assemble \
  -i "${RAW_READS}" \
  -o "${BASE_OUTPUT_DIR}" \
  -t "${SLURM_CPUS_PER_TASK}" \
  -s plant

echo "INFO: 'himt assemble' completed successfully."

# ==============================================================================
# Rename outputs to HyPR01/H. perforatum-specific names.
#
# HiMT has no --prefix/--name flag (confirmed against the full documented
# option list for `himt assemble`) -- its output filenames are hardcoded
# (himt_mitochondrial.gfa, himt_chloroplast.gfa, etc.) regardless of input
# filename or sample name. Copying (not moving, so the originals remain in
# case anything downstream still expects the default names) to sample-
# specific names here so every later step in this script, and anything
# referencing these files afterward (the audit docs, the comparison step
# below), uses an unambiguous HyPR01-specific name rather than a generic
# one that would collide if this pipeline is ever re-run on another sample.
# ==============================================================================
echo "======================================================================"
echo "INFO: Renaming HiMT outputs to HyPR01-specific filenames"
echo "======================================================================"

declare -A HIMT_RENAME_MAP=(
    ["himt_mitochondrial.gfa"]="HyPR01_HiMT_mitochondrion.gfa"
    ["himt_mitochondrial.html"]="HyPR01_HiMT_mitochondrion_QC.html"
    ["himt_chloroplast.gfa"]="HyPR01_HiMT_chloroplast.gfa"
    ["himt_chloroplast.html"]="HyPR01_HiMT_chloroplast_QC.html"
    ["extract.fa"]="HyPR01_HiMT_extracted_organelle_reads.fa"
)

for default_name in "${!HIMT_RENAME_MAP[@]}"; do
    src="${BASE_OUTPUT_DIR}/${default_name}"
    dst="${BASE_OUTPUT_DIR}/${HIMT_RENAME_MAP[$default_name]}"
    if [ -f "${src}" ]; then
        cp "${src}" "${dst}"
        echo "  ${default_name}  ->  ${HIMT_RENAME_MAP[$default_name]}"
    else
        echo "  NOT FOUND (check actual output layout, names may differ from README): ${default_name}"
    fi
done

# ==============================================================================
# Built-in QC reports (produced automatically by `himt assemble`, not a
# separate step -- surfacing their locations here so they aren't missed)
# ==============================================================================
echo "======================================================================"
echo "INFO: Built-in assembly QC reports (open these first):"
echo "======================================================================"
for report in "${BASE_OUTPUT_DIR}/HyPR01_HiMT_mitochondrion_QC.html" "${BASE_OUTPUT_DIR}/HyPR01_HiMT_chloroplast_QC.html"; do
    if [ -f "${report}" ]; then
        echo "  Found: ${report}"
    else
        echo "  NOT FOUND (check actual output layout): ${report}"
    fi
done

# ==============================================================================
# Reference genome acquisition: download the published H. perforatum
# mitochondrial genome (GenBase C_AA133241.1 + C_AA133242.1) for direct
# collinearity comparison via `himt compare` (future_workflow_plan.md, Stage 6).
#
# Source, version, and accession detail: see the REFERENCE_MITO_ACCESSIONS /
# REFERENCE_MITO_GENOME block near the top of this script for the full
# metadata (database, API endpoint, paper citation, verified lengths, and
# the H. ascyron mix-up this was checked against).
#
# Cached: only downloads if the combined reference file doesn't already
# exist, so re-running this script doesn't re-fetch it every time.
# ==============================================================================
echo "======================================================================"
echo "INFO: Reference genome acquisition (GenBase, NGDC/CNCB)"
echo "======================================================================"

mkdir -p "${REFERENCE_DIR}"

# This whole block is best-effort / non-fatal: it downloads a *bonus*
# comparison reference, not a required input for the assembly itself (which
# has already completed by this point in the script). REFERENCE_OK tracks
# whether it's safe to run the comparison step below; any failure here sets
# it to "no" and moves on rather than aborting the whole job. Every network
# call is guarded against `set -e` killing the script on a connection-level
# curl failure (as opposed to just a bad HTTP status) -- `|| echo "000"` /
# `|| true` ensure a failing command substitution doesn't trigger `set -e`.
REFERENCE_OK="yes"

if [ ! -s "${REFERENCE_MITO_GENOME}" ]; then
    echo "INFO: Reference not yet present -- downloading from GenBase."
    : > "${REFERENCE_MITO_GENOME}.tmp"

    for acc in "${REFERENCE_MITO_ACCESSIONS[@]}"; do
        acc_file="${REFERENCE_DIR}/${acc}.fasta"
        echo "  Fetching ${acc} from https://ngdc.cncb.ac.cn/genbase/api/file/fasta?acc=${acc}"

        http_code=$(curl -sL --max-time 60 \
            "https://ngdc.cncb.ac.cn/genbase/api/file/fasta?acc=${acc}" \
            -o "${acc_file}" -w "%{http_code}" || echo "000")

        if [ "${http_code}" != "200" ] || [ ! -s "${acc_file}" ]; then
            echo "WARNING: Download failed for ${acc} (HTTP ${http_code}). This compute node" >&2
            echo "         may not have outbound HTTPS access -- if so, download the" >&2
            echo "         reference on a login node instead and place the combined" >&2
            echo "         FASTA at ${REFERENCE_MITO_GENOME} before re-running just this section." >&2
            echo "         Skipping reference comparison for this run; assembly results above are unaffected." >&2
            REFERENCE_OK="no"
            break
        fi

        # Verify it's actually the expected species/molecule before trusting it --
        # this exact check caught a wrong-species accession (H. ascyron, not
        # H. perforatum) during script development. Never skip this.
        header=$(head -1 "${acc_file}" || echo "")
        if ! echo "${header}" | grep -qi "Hypericum perforatum"; then
            echo "WARNING: ${acc} header does not mention 'Hypericum perforatum':" >&2
            echo "         ${header}" >&2
            echo "         Refusing to use this file as the reference -- check the accession." >&2
            REFERENCE_OK="no"
            break
        fi
        echo "    Header : ${header}"

        seqlen=$(awk '/^>/{next}{n+=length($0)}END{print n}' "${acc_file}" || echo "?")
        echo "    Length : ${seqlen} bp"

        cat "${acc_file}" >> "${REFERENCE_MITO_GENOME}.tmp"
    done

    if [ "${REFERENCE_OK}" = "yes" ]; then
        mv "${REFERENCE_MITO_GENOME}.tmp" "${REFERENCE_MITO_GENOME}"
        total_len=$(awk '/^>/{next}{n+=length($0)}END{print n}' "${REFERENCE_MITO_GENOME}" || echo "?")
        echo "INFO: Combined reference written to ${REFERENCE_MITO_GENOME}"
        echo "INFO: Combined length: ${total_len} bp (expected 485,128 bp -- Chr1 204,492 + Chr2 280,636)"
        if [ "${total_len}" != "485128" ]; then
            echo "WARNING: Combined length does not match the expected 485,128 bp." >&2
            echo "         Proceeding anyway, but check ${REFERENCE_MITO_GENOME} manually." >&2
        fi

        # Record keeping: checksum + retrieval date, for the methods text
        { sha256sum "${REFERENCE_MITO_GENOME}" > "${REFERENCE_MITO_GENOME}.sha256"; } 2>/dev/null || \
            { shasum -a 256 "${REFERENCE_MITO_GENOME}" > "${REFERENCE_MITO_GENOME}.sha256"; } || true
        echo "Retrieved: $(date -u +'%Y-%m-%d %H:%M:%S UTC')" > "${REFERENCE_MITO_GENOME}.provenance.txt"
        echo "Source: GenBase (NGDC/CNCB), https://ngdc.cncb.ac.cn/genbase/" >> "${REFERENCE_MITO_GENOME}.provenance.txt"
        echo "Accessions: ${REFERENCE_MITO_ACCESSIONS[*]}" >> "${REFERENCE_MITO_GENOME}.provenance.txt"
        echo "Paper: https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2026.1763082/full" >> "${REFERENCE_MITO_GENOME}.provenance.txt"
        echo "INFO: Provenance recorded in ${REFERENCE_MITO_GENOME}.provenance.txt and .sha256"
    else
        rm -f "${REFERENCE_MITO_GENOME}.tmp"
    fi
else
    echo "INFO: Reference already present at ${REFERENCE_MITO_GENOME} -- skipping download."
fi

# ==============================================================================
# Direct collinearity comparison against the published reference
# (future_workflow_plan.md, Stage 6)
# ==============================================================================
echo "======================================================================"
echo "INFO: Running reference comparison..."
echo "======================================================================"

HIMT_MITO_RESULT="${BASE_OUTPUT_DIR}/HyPR01_HiMT_mitochondrion.gfa"

if [ -f "${REFERENCE_MITO_GENOME}" ] && [ -f "${HIMT_MITO_RESULT}" ]; then
    echo "INFO: Running 'himt compare' against published H. perforatum mitogenome."
    COMPARE_OUT_DIR="${BASE_OUTPUT_DIR}/reference_comparison"
    mkdir -p "${COMPARE_OUT_DIR}"
    himt compare \
      -r "${REFERENCE_MITO_GENOME}" \
      -q "${HIMT_MITO_RESULT}" \
      -o "${COMPARE_OUT_DIR}" \
      -c mitochondrial
    echo "INFO: Comparison output written to ${COMPARE_OUT_DIR}"
else
    echo "WARNING: Skipping comparison -- expected HiMT output not found at:"
    echo "         ${HIMT_MITO_RESULT}"
    echo "         (Check the actual output filename once 'himt assemble' has run once --"
    echo "          this name is taken from the README's documented output table and"
    echo "          hasn't been confirmed against a real run yet.)"
fi

# ==============================================================================
# Output Inspection
# ==============================================================================
echo "======================================================================"
echo "INFO: Summary of generated assembly files:"
echo "======================================================================"

find "${BASE_OUTPUT_DIR}" -maxdepth 2 \( -name "*.fasta" -o -name "*.fa" -o -name "*.gfa" -o -name "*.html" \) -exec ls -lh {} + || true

echo "======================================================================"
echo "INFO: Job completed successfully at $(date)"
echo "======================================================================"
