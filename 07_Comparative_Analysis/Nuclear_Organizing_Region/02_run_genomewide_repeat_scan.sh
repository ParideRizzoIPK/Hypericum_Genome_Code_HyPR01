#!/bin/bash
#
# ==============================================================================
# HyPR01 Genome-Wide Repeat Scan (Chr01-Chr08, both haplotypes)
# ==============================================================================
# Extends the Chr08 NOR confirmation (run_NOR_confirmation_pipeline.sh,
# HyPR01_Chr08_NOR_Confirmation_Report.md) genome-wide: are there other rDNA
# loci / large tandem repeat arrays anywhere else in the assembly, on either
# haplotype?
#
# Two tools, both already validated at whole-chromosome scale in the Chr08
# work, are run on all 16 chromosome sequences (Chr01-Chr08 x Hap1/Hap2):
#   1. barrnap-HGV --kingdom euk   - direct rRNA gene (18S/5.8S/28S/5S) calls
#   2. TRF                         - tandem repeat arrays (period <=2000bp;
#                                     this catches satellite/microsatellite
#                                     arrays such as centromeric repeats, NOT
#                                     the ~10.3kb rDNA monomer itself, which
#                                     barrnap and the Chr08-specific nucmer
#                                     analysis already characterize)
#
# IMPORTANT - DATA HYGIENE: this script never reads, writes, or indexes the
# original genome/GFF3 files in place. It only ever operates on symlinks and
# extracted copies inside genomewide_repeat_scan/. The original
# hap1.masked_ChrUn.fa / hap2.masked_ChrUn.fa and the GFF3 annotations are
# never modified and never receive a companion .fai index file.
#
# Requires the `nor_test` mamba environment and barrnap-HGV, both already set
# up by run_NOR_confirmation_pipeline.sh. This script does not repeat that
# setup; run run_NOR_confirmation_pipeline.sh first if starting fresh.
#
# Output (all under genomewide_repeat_scan/, nothing written elsewhere):
#   chromosomes/Hap{1,2}_Chr0{1-8}.fa            Extracted per-chromosome FASTAs
#   barrnap/Hap{1,2}_Chr0{1-8}_barrnap_euk.gff3  Per-chromosome rRNA gene calls
#   trf/Hap{1,2}_Chr0{1-8}.fa.*.dat              Per-chromosome TRF tandem repeats
#   summary/rRNA_gene_summary.tsv                Per-chromosome rRNA gene tally
#   summary/TRF_top_repeats_summary.tsv          Largest tandem arrays per chromosome
# ==============================================================================

set -euo pipefail
trap 'echo "[FATAL] line ${LINENO} exited with status $?" >&2' ERR

BASE_DIR="/path/to/your/directory"
SCAN_DIR="${BASE_DIR}/genomewide_repeat_scan"
CHROM_DIR="${SCAN_DIR}/chromosomes"
BARRNAP_DIR="${SCAN_DIR}/barrnap"
TRF_DIR="${SCAN_DIR}/trf"
SUMMARY_DIR="${SCAN_DIR}/summary"

ENV_NAME="nor_test"
BARRNAP="${BASE_DIR}/analysis/phyloFlash/barrnap-HGV/bin/barrnap_HGV"

# Original, untouched source genomes
HAP1_GENOME="${BASE_DIR}/hap1.masked_ChrUn.fa"
HAP2_GENOME="${BASE_DIR}/hap2.masked_ChrUn.fa"

mkdir -p "${CHROM_DIR}" "${BARRNAP_DIR}" "${TRF_DIR}" "${SUMMARY_DIR}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

if [ ! -f "${BARRNAP}" ]; then
    echo "[FATAL] barrnap-HGV not found at ${BARRNAP}." >&2
    echo "        Run run_NOR_confirmation_pipeline.sh first (it retrieves and patches it)." >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# STEP 1: SYMLINK + INDEX GENOMES INSIDE THE SCAN DIRECTORY ONLY
# ------------------------------------------------------------------------------
echo "=========================================================================="
echo "[$(date)] STEP 1: SYMLINKING + INDEXING GENOMES (inside ${CHROM_DIR} only)"
echo "=========================================================================="

cd "${CHROM_DIR}"
[ -L hap1.masked_ChrUn.fa ] || ln -s "${HAP1_GENOME}" hap1.masked_ChrUn.fa
[ -L hap2.masked_ChrUn.fa ] || ln -s "${HAP2_GENOME}" hap2.masked_ChrUn.fa
seqkit faidx hap1.masked_ChrUn.fa
seqkit faidx hap2.masked_ChrUn.fa

# ------------------------------------------------------------------------------
# STEP 2: EXTRACT ALL 16 CHROMOSOME FASTAS
# ------------------------------------------------------------------------------
echo "=========================================================================="
echo "[$(date)] STEP 2: EXTRACTING Chr01-Chr08 FOR BOTH HAPLOTYPES"
echo "=========================================================================="

for HAP in 1 2; do
    for CHR in 01 02 03 04 05 06 07 08; do
        OUT="Hap${HAP}_Chr${CHR}.fa"
        if [ ! -s "${OUT}" ]; then
            seqkit faidx "hap${HAP}.masked_ChrUn.fa" "Chr${CHR}" > "${OUT}"
        fi
    done
done
seqkit stats Hap*_Chr*.fa

# ------------------------------------------------------------------------------
# STEP 3: barrnap-HGV --kingdom euk ON ALL 16 CHROMOSOMES
# ------------------------------------------------------------------------------
# KNOWN ISSUE: HMMER's nhmmer auto-detects the target sequence's alphabet by
# sampling an initial window. Several HyPR01 chromosomes are assembled right
# up to the telomere, and the plant telomere repeat in this strand orientation
# ((CCCTAAA)n) uses only A/C/T - zero G - for several kb. If that low-diversity
# run exceeds nhmmer's sampling window, alphabet detection fails outright
# ("Unable to guess alphabet for target sequence database file"). This is a
# tool limitation triggered by genuine, correctly-assembled telomeric sequence,
# not a data problem. Fix: for any chromosome where this occurs, barrnap is run
# on a copy with the low-diversity leading run trimmed off (first position
# where a full A/C/G/T window appears), and the trim offset is added back to
# every reported coordinate afterwards so final coordinates are always in the
# original, untrimmed chromosome's coordinate system. The trim offset used for
# each affected chromosome is recorded alongside its output for full
# transparency. No file under chromosomes/ used elsewhere in this pipeline is
# ever modified - trimming only ever happens on a throwaway temp copy.
# ------------------------------------------------------------------------------
echo "=========================================================================="
echo "[$(date)] STEP 3: barrnap-HGV EUKARYOTIC rRNA SCAN, ALL 16 CHROMOSOMES"
echo "=========================================================================="

find_diversity_offset() {
    # Prints the 0-based offset of the first 5000bp window containing all of
    # A/C/G/T, scanning only the first 50kb (telomeric runs are always << this).
    # Prints 0 if the sequence is already diverse from the start.
    python3 - "$1" << 'PYEOF'
import sys
path = sys.argv[1]
window, scan_limit = 5000, 50000
with open(path) as f:
    f.readline()
    seq, total = [], 0
    for line in f:
        seq.append(line.strip())
        total += len(line.strip())
        if total >= scan_limit:
            break
seq = ''.join(seq).upper()
offset = 0
for start in range(0, min(len(seq), scan_limit) - window):
    if {'A', 'C', 'G', 'T'}.issubset(set(seq[start:start + window])):
        offset = start
        break
print(offset)
PYEOF
}

for HAP in 1 2; do
    for CHR in 01 02 03 04 05 06 07 08; do
        NAME="Hap${HAP}_Chr${CHR}"
        OUT="${BARRNAP_DIR}/${NAME}_barrnap_euk.gff3"
        [ -s "${OUT}" ] && continue

        echo "[$(date)] barrnap: ${NAME}"
        SRC="${CHROM_DIR}/${NAME}.fa"
        LOG="${BARRNAP_DIR}/${NAME}_barrnap_euk.log"

        if perl "${BARRNAP}" --kingdom euk --threads 8 --quiet "${SRC}" > "${OUT}" 2> "${LOG}"; then
            echo "0" > "${BARRNAP_DIR}/${NAME}_trim_offset.txt"
            continue
        fi

        if grep -q "Unable to guess alphabet" "${LOG}"; then
            # +15000bp safety margin beyond the first fully-diverse 5kb window:
            # empirically, nhmmer's own internal sampling needs more headroom
            # than a single window edge-case provides (verified by direct testing).
            OFFSET=$(( $(find_diversity_offset "${SRC}") + 15000 ))
            echo "[$(date)]   -> alphabet-guess failure (telomeric leading run); retrying with ${OFFSET} bp trimmed, coordinates will be shifted back by +${OFFSET}"
            TRIMMED="${CHROM_DIR}/${NAME}_barrnap_trimmed_tmp.fa"
            python3 - "${SRC}" "${OFFSET}" "${TRIMMED}" << 'PYEOF'
import sys
path, offset, outpath = sys.argv[1], int(sys.argv[2]), sys.argv[3]
with open(path) as f:
    header = f.readline().strip()
    seq = ''.join(l.strip() for l in f)
trimmed = seq[offset:]
with open(outpath, 'w') as out:
    out.write(header + "\n")
    for i in range(0, len(trimmed), 60):
        out.write(trimmed[i:i+60] + "\n")
PYEOF
            RAW_OUT="${OUT}.rawtrimmed"
            if ! perl "${BARRNAP}" --kingdom euk --threads 8 --quiet "${TRIMMED}" > "${RAW_OUT}" 2>> "${LOG}"; then
                echo "[FATAL] barrnap still failed on ${NAME} after trimming ${OFFSET}bp - see ${LOG}" >&2
                rm -f "${TRIMMED}" "${RAW_OUT}" "${OUT}"
                exit 1
            fi
            awk -v off="${OFFSET}" 'BEGIN{OFS="\t"} /^#/{print; next} {$4+=off; $5+=off; print}' "${RAW_OUT}" > "${OUT}"
            echo "${OFFSET}" > "${BARRNAP_DIR}/${NAME}_trim_offset.txt"
            rm -f "${TRIMMED}" "${RAW_OUT}"
        else
            echo "[FATAL] barrnap failed on ${NAME} for a reason other than alphabet guessing - see ${LOG}" >&2
            exit 1
        fi
    done
done

# ------------------------------------------------------------------------------
# STEP 4: TRF ON ALL 16 CHROMOSOMES
# ------------------------------------------------------------------------------
echo "=========================================================================="
echo "[$(date)] STEP 4: TRF TANDEM REPEAT SCAN, ALL 16 CHROMOSOMES"
echo "=========================================================================="

cd "${TRF_DIR}"
for HAP in 1 2; do
    for CHR in 01 02 03 04 05 06 07 08; do
        NAME="Hap${HAP}_Chr${CHR}"
        DATFILE="${NAME}.fa.2.7.7.80.10.50.2000.dat"
        if [ ! -s "${DATFILE}" ]; then
            echo "[$(date)] TRF: ${NAME}"
            trf "${CHROM_DIR}/${NAME}.fa" 2 7 7 80 10 50 2000 -d -h -l 5 || true
        fi
    done
done

# ------------------------------------------------------------------------------
# STEP 5: SUMMARY TABLES
# ------------------------------------------------------------------------------
echo "=========================================================================="
echo "[$(date)] STEP 5: BUILDING SUMMARY TABLES"
echo "=========================================================================="

# --- rRNA gene tally per chromosome ---
RRNA_SUMMARY="${SUMMARY_DIR}/rRNA_gene_summary.tsv"
echo -e "haplotype\tchrom\tgene\tn_calls\tn_complete\tn_partial\tfirst_pos\tlast_pos" > "${RRNA_SUMMARY}"
for HAP in 1 2; do
    for CHR in 01 02 03 04 05 06 07 08; do
        NAME="Hap${HAP}_Chr${CHR}"
        GFF="${BARRNAP_DIR}/${NAME}_barrnap_euk.gff3"
        [ -s "${GFF}" ] || continue
        for GENE in 18S_rRNA 5_8S_rRNA 28S_rRNA 5S_rRNA; do
            LINES=$(awk -F'\t' -v g="Name=${GENE}" '$9 ~ g' "${GFF}")
            N=$(echo "${LINES}" | grep -c . || true)
            [ "${N}" -eq 0 ] && continue
            NPART=$(echo "${LINES}" | grep -c "partial" || true)
            NCOMP=$((N - NPART))
            FIRST=$(echo "${LINES}" | awk -F'\t' '{print $4}' | sort -n | head -1)
            LAST=$(echo "${LINES}" | awk -F'\t' '{print $5}' | sort -n | tail -1)
            echo -e "Hap${HAP}\tChr${CHR}\t${GENE}\t${N}\t${NCOMP}\t${NPART}\t${FIRST}\t${LAST}" >> "${RRNA_SUMMARY}"
        done
    done
done
echo "[+] Saved ${RRNA_SUMMARY}"

# --- Largest TRF tandem arrays per chromosome (top 5 by period*copies = approx array footprint) ---
TRF_SUMMARY="${SUMMARY_DIR}/TRF_top_repeats_summary.tsv"
echo -e "haplotype\tchrom\tstart\tend\tperiod\tcopies\tconsensus_size\tpct_match\tscore\tarray_span_bp" > "${TRF_SUMMARY}"
for HAP in 1 2; do
    for CHR in 01 02 03 04 05 06 07 08; do
        NAME="Hap${HAP}_Chr${CHR}"
        DATFILE="${TRF_DIR}/${NAME}.fa.2.7.7.80.10.50.2000.dat"
        [ -s "${DATFILE}" ] || continue
        # Written to a temp file rather than piped into head: head reading a
        # live pipe from a large producer can SIGPIPE the producer once it has
        # enough lines, which trips set -o pipefail. Reading from a completed
        # file has no such race.
        SORTED_TMP=$(mktemp)
        awk -v hap="Hap${HAP}" -v chr="Chr${CHR}" \
            'NF>=15 && $1 ~ /^[0-9]+$/ {span=$2-$1+1; print hap"\t"chr"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8"\t"span}' \
            "${DATFILE}" | sort -k10,10 -rn > "${SORTED_TMP}"
        head -5 "${SORTED_TMP}" >> "${TRF_SUMMARY}"
        rm -f "${SORTED_TMP}"
    done
done
echo "[+] Saved ${TRF_SUMMARY}"

echo "=========================================================================="
echo "[$(date)] GENOME-WIDE REPEAT SCAN COMPLETE"
echo "=========================================================================="
