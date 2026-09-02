#!/bin/bash
#
# ==============================================================================
# HyPR01 Chr08 NOR Confirmation Pipeline
# ==============================================================================
# Tests the hypothesis that the ~9x baseline read-depth anomaly near the
# terminus of Chr08 (both HyPR01 haplotypes; see HyPR01_Coverage_Pipeline_Report.md
# section 8) is a 45S ribosomal DNA array / nucleolar organizer region (NOR).
#
# Four independent, orthogonal lines of evidence are generated:
#   1. barrnap-HGV  - direct eukaryotic 18S/5.8S/28S rRNA gene prediction
#   2. BLASTn       - identity against a real, same-species (H. perforatum) rDNA
#                     GenBank reference
#   3. nucmer       - reference-free, model-free self-alignment dotplot
#   4. TRF          - independent short-tandem-repeat confirmation
#
# Full write-up of results and interpretation:
#   HyPR01_Chr08_NOR_Confirmation_Report.md
#
# Tested on: macOS, Apple Silicon (arm64). Requires network access (bioconda,
# GitHub, NCBI E-utilities). Idempotent - safe to re-run; steps whose expected
# output already exists are skipped.
#
# Deliverables (all written to ${ANALYSIS_DIR}):
#   Hap{1,2}_Chr08.fa                         Full Chr08, extracted per haplotype
#   Hap1_Chr08_NORcandidate_1-120kb.fa        Hap1 NOR-candidate region
#   Hap2_Chr08_NORcandidate_1-1050kb.fa       Hap2 NOR-candidate region (broad)
#   Hap2_Chr08_1-200kb.fa                     Hap2 NOR-candidate region (BLAST/dotplot)
#   Hap{1,2}_Chr08_barrnap_euk.gff3           barrnap-HGV eukaryotic rRNA gene calls
#   reference_rDNA/ON685357.1.fa              H. perforatum rDNA GenBank reference
#   reference_rDNA/HpITS_db.*                 BLAST database built from the above
#   Hap1_vs_HpITS.blast.tsv, Hap2_vs_HpITS.blast.tsv   BLASTn results
#   Hap{1,2}_selfdot.delta / .coords          nucmer self-alignment output
#   Chr08_NOR_selfdotplot.png                 Rendered self-dotplot figure
#   trf_out/*.dat                             TRF tandem-repeat calls
# ==============================================================================

set -euo pipefail
trap 'echo "[FATAL] line ${LINENO} exited with status $?" >&2' ERR

# ------------------------------------------------------------------------------
# STEP 0: CONFIGURATION
# ------------------------------------------------------------------------------
BASE_DIR="/path/to/your/directory"
ANALYSIS_DIR="${BASE_DIR}/analysis"
ENV_NAME="nor_test"

HAP1_GENOME="${BASE_DIR}/hap1.masked_ChrUn.fa"
HAP2_GENOME="${BASE_DIR}/hap2.masked_ChrUn.fa"

# Same-species (H. perforatum) rDNA reference used for BLAST identity confirmation
# (18S partial - ITS1 - 5.8S - ITS2 - 28S partial, 3730 bp)
NCBI_REF_ACC="ON685357.1"

PHYLOFLASH_DIR="${ANALYSIS_DIR}/phyloFlash"
BARRNAP="${PHYLOFLASH_DIR}/barrnap-HGV/bin/barrnap_HGV"
NHMMER_BIN="${PHYLOFLASH_DIR}/barrnap-HGV/binaries/darwin/nhmmer"

mkdir -p "${ANALYSIS_DIR}" "${ANALYSIS_DIR}/reference_rDNA" "${ANALYSIS_DIR}/trf_out"

# ------------------------------------------------------------------------------
# STEP 1: MAMBA ENVIRONMENT
# ------------------------------------------------------------------------------
echo "=========================================================================="
echo "[$(date)] STEP 1: MAMBA ENVIRONMENT (${ENV_NAME})"
echo "=========================================================================="

if ! mamba env list | grep -qE "^${ENV_NAME}\s"; then
    echo "[$(date)] Creating ${ENV_NAME}..."
    mamba create -y -n "${ENV_NAME}" -c conda-forge -c bioconda \
        barrnap hmmer mummer4 trf blast samtools seqkit python=3.11 matplotlib numpy
else
    echo "[$(date)] Environment ${ENV_NAME} already exists. Skipping creation."
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

echo "[$(date)] Tool versions:"
echo "  hmmer:    $(nhmmer -h 2>&1 | grep -m1 HMMER)"
echo "  blastn:   $(blastn -version | head -1)"
echo "  mummer4:  $(nucmer --version 2>&1 | head -1)"
echo "  seqkit:   $(seqkit version 2>&1 | head -1)"

# ------------------------------------------------------------------------------
# STEP 2: barrnap-HGV (eukaryotic-capable rRNA predictor)
# ------------------------------------------------------------------------------
echo "=========================================================================="
echo "[$(date)] STEP 2: RETRIEVE barrnap-HGV (bioconda barrnap has no euk HMMs;"
echo "         full phyloflash package fails to solve on osx-arm64 due to emirge)"
echo "=========================================================================="

if [ ! -f "${BARRNAP}" ]; then
    echo "[$(date)] Cloning phyloFlash (shallow) for its bundled barrnap-HGV fork..."
    git clone --depth 1 https://github.com/HRGV/phyloFlash.git "${PHYLOFLASH_DIR}"
else
    echo "[$(date)] barrnap-HGV already present. Skipping clone."
fi

# The bundled binaries/darwin/nhmmer is an old x86_64 build. On Apple Silicon,
# swap in our conda environment's native arm64 nhmmer, fix its rpath (the
# copied binary's embedded @loader_path/../lib/ no longer resolves from its
# new location), and re-sign (install_name_tool invalidates the code
# signature, and arm64 binaries must be signed to execute at all).
if [[ "$(uname -s)" == "Darwin" && "$(uname -m)" == "arm64" ]]; then
    if ! file "${NHMMER_BIN}" 2>/dev/null | grep -q "arm64"; then
        echo "[$(date)] Patching barrnap-HGV's bundled nhmmer for Apple Silicon..."
        cp "$(command -v nhmmer)" "${NHMMER_BIN}"
        chmod +x "${NHMMER_BIN}"
        install_name_tool -add_rpath "${CONDA_PREFIX}/lib" "${NHMMER_BIN}"
        codesign --force --sign - "${NHMMER_BIN}"
    else
        echo "[$(date)] nhmmer binary already patched (arm64). Skipping."
    fi
fi

# ------------------------------------------------------------------------------
# STEP 3: EXTRACT Chr08 AND NOR-CANDIDATE SUBREGIONS
# ------------------------------------------------------------------------------
echo "=========================================================================="
echo "[$(date)] STEP 3: EXTRACTING Chr08 AND NOR-CANDIDATE SUBREGIONS"
echo "=========================================================================="

cd "${ANALYSIS_DIR}"

seqkit faidx "${HAP1_GENOME}"
seqkit faidx "${HAP2_GENOME}"

if [ ! -s Hap1_Chr08.fa ]; then
    seqkit faidx "${HAP1_GENOME}" Chr08 > Hap1_Chr08.fa
fi
if [ ! -s Hap2_Chr08.fa ]; then
    seqkit faidx "${HAP2_GENOME}" Chr08 > Hap2_Chr08.fa
fi

if [ ! -s Hap1_Chr08_NORcandidate_1-120kb.fa ]; then
    seqkit subseq --chr Chr08 -r 1:120000 "${HAP1_GENOME}" > Hap1_Chr08_NORcandidate_1-120kb.fa
fi
if [ ! -s Hap2_Chr08_NORcandidate_1-1050kb.fa ]; then
    seqkit subseq --chr Chr08 -r 1:1050000 "${HAP2_GENOME}" > Hap2_Chr08_NORcandidate_1-1050kb.fa
fi
if [ ! -s Hap2_Chr08_1-200kb.fa ]; then
    seqkit subseq --chr Chr08 -r 1:200000 "${HAP2_GENOME}" > Hap2_Chr08_1-200kb.fa
fi

seqkit stats Hap1_Chr08.fa Hap2_Chr08.fa Hap1_Chr08_NORcandidate_1-120kb.fa \
    Hap2_Chr08_NORcandidate_1-1050kb.fa Hap2_Chr08_1-200kb.fa

# ------------------------------------------------------------------------------
# STEP 4: barrnap-HGV EUKARYOTIC rRNA GENE PREDICTION (evidence line 1)
# ------------------------------------------------------------------------------
echo "=========================================================================="
echo "[$(date)] STEP 4: barrnap-HGV --kingdom euk ON FULL Chr08 (BOTH HAPLOTYPES)"
echo "=========================================================================="

if [ ! -s Hap1_Chr08_barrnap_euk.gff3 ]; then
    perl "${BARRNAP}" --kingdom euk --threads 8 --quiet Hap1_Chr08.fa \
        > Hap1_Chr08_barrnap_euk.gff3 2> Hap1_Chr08_barrnap_euk.log
fi
if [ ! -s Hap2_Chr08_barrnap_euk.gff3 ]; then
    perl "${BARRNAP}" --kingdom euk --threads 8 --quiet Hap2_Chr08.fa \
        > Hap2_Chr08_barrnap_euk.gff3 2> Hap2_Chr08_barrnap_euk.log
fi

echo "[$(date)] Hap1 rRNA gene calls:"
grep -c $'\trRNA\t' Hap1_Chr08_barrnap_euk.gff3 || true
echo "[$(date)] Hap2 rRNA gene calls:"
grep -c $'\trRNA\t' Hap2_Chr08_barrnap_euk.gff3 || true

# ------------------------------------------------------------------------------
# STEP 5: SAME-SPECIES rDNA REFERENCE + BLASTn (evidence line 2)
# ------------------------------------------------------------------------------
echo "=========================================================================="
echo "[$(date)] STEP 5: FETCH H. perforatum rDNA REFERENCE (${NCBI_REF_ACC}) + BLASTn"
echo "=========================================================================="

REF_FA="reference_rDNA/${NCBI_REF_ACC}.fa"
if [ ! -s "${REF_FA}" ]; then
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${NCBI_REF_ACC}&rettype=fasta&retmode=text" \
        > "${REF_FA}"
fi

if [ ! -s reference_rDNA/HpITS_db.nsq ]; then
    makeblastdb -in "${REF_FA}" -dbtype nucl -out reference_rDNA/HpITS_db
fi

blastn -query Hap1_Chr08_NORcandidate_1-120kb.fa -db reference_rDNA/HpITS_db \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -max_target_seqs 20 | sort -k4,4 -rn > Hap1_vs_HpITS.blast.tsv

blastn -query Hap2_Chr08_1-200kb.fa -db reference_rDNA/HpITS_db \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -max_target_seqs 20 | sort -k4,4 -rn > Hap2_vs_HpITS.blast.tsv

echo "[$(date)] Hap1 BLAST hits >=3000bp (near-full-length rDNA reference matches):"
awk '$4>=3000' Hap1_vs_HpITS.blast.tsv | wc -l
echo "[$(date)] Hap2 BLAST hits >=3000bp:"
awk '$4>=3000' Hap2_vs_HpITS.blast.tsv | wc -l

# ------------------------------------------------------------------------------
# STEP 6: nucmer SELF-ALIGNMENT / DOTPLOT (evidence line 3, model-free)
# ------------------------------------------------------------------------------
echo "=========================================================================="
echo "[$(date)] STEP 6: nucmer SELF-ALIGNMENT (MODEL-FREE TANDEM REPEAT CHECK)"
echo "=========================================================================="

if [ ! -s Hap1_selfdot.coords ]; then
    nucmer --maxmatch --nosimplify -p Hap1_selfdot \
        Hap1_Chr08_NORcandidate_1-120kb.fa Hap1_Chr08_NORcandidate_1-120kb.fa
    show-coords -rclT Hap1_selfdot.delta > Hap1_selfdot.coords
fi
if [ ! -s Hap2_selfdot.coords ]; then
    nucmer --maxmatch --nosimplify -p Hap2_selfdot \
        Hap2_Chr08_1-200kb.fa Hap2_Chr08_1-200kb.fa
    show-coords -rclT Hap2_selfdot.delta > Hap2_selfdot.coords
fi

echo "[$(date)] Dominant repeat-offset periodicity (Hap1, long/high-identity hits only):"
awk 'NR>4 && $1<$3 && $7>90 && ($2-$1)>2000 {print $3-$1}' Hap1_selfdot.coords | sort -n

# Render the dotplot figure
python3 << 'PYEOF'
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load_coords(path):
    pts = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i < 4:
                continue
            parts = line.split('\t') if '\t' in line else line.split()
            try:
                s1, e1, s2, e2 = int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3])
                ident = float(parts[6])
            except (ValueError, IndexError):
                continue
            pts.append((s1, e1, s2, e2, ident))
    return pts

fig, axes = plt.subplots(1, 2, figsize=(13, 6.2))

for ax, (label, coordfile, xmax) in zip(
    axes,
    [("Hap1 Chr08:1-120,000 (self vs self)", "Hap1_selfdot.coords", 120000),
     ("Hap2 Chr08:1-200,000 (self vs self)", "Hap2_selfdot.coords", 200000)]
):
    pts = load_coords(coordfile)
    for s1, e1, s2, e2, ident in pts:
        if s1 == s2 and e1 == e2:
            continue  # skip the trivial identity diagonal
        alpha = min(1.0, max(0.15, (ident - 75) / 25))
        ax.plot([s1, e1], [s2, e2], color="#1f77b4", lw=1.1, alpha=alpha)
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, xmax)
    ax.set_aspect("equal")
    ax.set_xlabel("Position (bp)")
    ax.set_ylabel("Position (bp)")
    ax.set_title(label, fontsize=10, fontweight="bold")
    ax.grid(alpha=0.25, lw=0.4)

fig.suptitle("HyPR01 Chr08 NOR-candidate region - nucmer self-alignment dotplot\n"
             "(off-diagonal ladder = tandem repeat; rung spacing = monomer length ~10,347 bp)",
             fontsize=11, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.90))
fig.savefig("Chr08_NOR_selfdotplot.png", dpi=200)
print("[+] saved Chr08_NOR_selfdotplot.png")
PYEOF

# ------------------------------------------------------------------------------
# STEP 7: TRF - COMPLEMENTARY INTERNAL SUB-REPEAT CHECK (evidence line 4)
# ------------------------------------------------------------------------------
echo "=========================================================================="
echo "[$(date)] STEP 7: TRF (complementary check; MaxPeriod capped at 2000bp, so"
echo "         this cannot see the primary ~10.3kb monomer directly - it finds"
echo "         shorter internal sub-repeats within it instead)"
echo "=========================================================================="

cd "${ANALYSIS_DIR}/trf_out"
trf ../Hap1_Chr08_NORcandidate_1-120kb.fa 2 7 7 80 10 50 2000 -d -h -l 1 || true
trf ../Hap2_Chr08_1-200kb.fa 2 7 7 80 10 50 2000 -d -h -l 1 || true

echo "[$(date)] Top TRF hits by score (Hap1):"
DATFILE=$(ls Hap1_Chr08_NORcandidate_1-120kb.fa.*.dat 2>/dev/null | head -1)
if [ -n "${DATFILE}" ]; then
    awk 'NF>=15 && $1 ~ /^[0-9]+$/ {print}' "${DATFILE}" | sort -k8,8 -rn | head -10 \
        | awk '{print "period="$3"  copies="$4"  score="$8"  start="$1"  end="$2}'
fi

echo "=========================================================================="
echo "[$(date)] NOR CONFIRMATION PIPELINE COMPLETE"
echo "See HyPR01_Chr08_NOR_Confirmation_Report.md for full interpretation."
echo "=========================================================================="
