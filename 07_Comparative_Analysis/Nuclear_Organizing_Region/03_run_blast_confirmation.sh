#!/bin/bash
#
# ==============================================================================
# BLAST Confirmation: Chr07 solitary rDNA unit + Chr01/Chr06 5S loci
# ==============================================================================
# Follow-up #1 from genomewide_repeat_scan/HyPR01_Genomewide_Repeat_Landscape_Report.md
# section 10: bring the newly-discovered dispersed rDNA unit (Chr07) and the
# two major 5S rDNA loci (Chr01, Chr06) to the same standard of evidence as the
# Chr08 NOR (same-species BLASTn identity, not just HMM/TRF pattern-matching).
#
# References used:
#   - 45S (18S-ITS1-5.8S-ITS2-28S): ON685357.1, H. perforatum (already fetched
#     for the Chr08 confirmation; reused here via the existing HpITS_db)
#   - 5S nuclear:  DQ110884.1 (H. perforatum), DQ110883.1 (H. androsaemum),
#                  DQ110885.1 (H. glandulosum) - congeners included because
#                  nuclear 5S rDNA is highly conserved even between species
#   - 5S plastid:  JX661983.1 (H. perforatum chloroplast rrn5) - included
#                  deliberately as a NEGATIVE CONTROL: a genuine nuclear 5S
#                  array should match the nuclear references well and the
#                  plastid-type sequence poorly/not at all, distinguishing a
#                  true nuclear 5S locus from a NUMT-like plastid insertion.
#
# DATA HYGIENE: reads only from genomewide_repeat_scan/chromosomes/ (already
# extracted, symlink-based copies from the genome-wide scan - never the
# original hap1/hap2 genome files). All new output goes under
# genomewide_repeat_scan/blast_confirmation/.
# ==============================================================================

set -euo pipefail
trap 'echo "[FATAL] line ${LINENO} exited with status $?" >&2' ERR

BASE_DIR="/path/to/your/directory"
SCAN_DIR="${BASE_DIR}/genomewide_repeat_scan"
BC_DIR="${SCAN_DIR}/blast_confirmation"
REGIONS_DIR="${BC_DIR}/regions"
REF5S_DIR="${BC_DIR}/reference_5S"
RESULTS_DIR="${BC_DIR}/blast_results"
HPITS_DB="${BASE_DIR}/analysis/reference_rDNA/HpITS_db"

mkdir -p "${REGIONS_DIR}" "${REF5S_DIR}" "${RESULTS_DIR}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nor_test

# ------------------------------------------------------------------------------
# STEP 1: FETCH 5S REFERENCE SEQUENCES + BUILD BLAST DB
# ------------------------------------------------------------------------------
echo "=========================================================================="
echo "[$(date)] STEP 1: FETCHING 5S rDNA REFERENCES"
echo "=========================================================================="

cd "${REF5S_DIR}"
for ACC in DQ110884.1 DQ110883.1 DQ110885.1 JX661983.1; do
    if [ ! -s "${ACC}.fa" ]; then
        curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${ACC}&rettype=fasta&retmode=text" > "${ACC}.fa"
        sleep 0.4
    fi
done

if [ ! -s Hp5S_db.nsq ]; then
    cat DQ110884.1.fa DQ110883.1.fa DQ110885.1.fa JX661983.1.fa > Hp_5S_combined.fa
    makeblastdb -in Hp_5S_combined.fa -dbtype nucl -out Hp5S_db
fi

# ------------------------------------------------------------------------------
# STEP 2: EXTRACT TARGET REGIONS (from genomewide_repeat_scan/chromosomes/ only)
# ------------------------------------------------------------------------------
echo "=========================================================================="
echo "[$(date)] STEP 2: EXTRACTING TARGET REGIONS"
echo "=========================================================================="

CHROM_DIR="${SCAN_DIR}/chromosomes"
cd "${REGIONS_DIR}"

extract () {
    local out="$1" chrom="$2" range="$3" src="$4"
    [ -s "${out}" ] || seqkit subseq --chr "${chrom}" -r "${range}" "${src}" > "${out}"
}

# Chr07 solitary rDNA unit (+/- ~2kb margin around the barrnap-called triplet)
extract Hap1_Chr07_rDNAunit.fa Chr07 18300000:18309000 "${CHROM_DIR}/Hap1_Chr07.fa"
extract Hap2_Chr07_rDNAunit.fa Chr07 18786000:18795000 "${CHROM_DIR}/Hap2_Chr07.fa"

# Chr06 (Hap1-only) dispersed 18S+5.8S pair, no 28S (+/- ~2kb margin)
extract Hap1_Chr06_dispersedpair.fa Chr06 16146000:16153000 "${CHROM_DIR}/Hap1_Chr06.fa"

# Chr01 5S: main pericentromeric cluster complex + distal secondary cluster
extract Hap1_Chr01_5S_main.fa   Chr01 20150000:22180000 "${CHROM_DIR}/Hap1_Chr01.fa"
extract Hap1_Chr01_5S_distal.fa Chr01 48610000:48665000 "${CHROM_DIR}/Hap1_Chr01.fa"
extract Hap2_Chr01_5S_main.fa   Chr01 21350000:23200000 "${CHROM_DIR}/Hap2_Chr01.fa"
extract Hap2_Chr01_5S_distal.fa Chr01 50190000:50280000 "${CHROM_DIR}/Hap2_Chr01.fa"

# Chr06 5S locus
extract Hap1_Chr06_5S.fa Chr06 25280000:26200000 "${CHROM_DIR}/Hap1_Chr06.fa"
extract Hap2_Chr06_5S.fa Chr06 26365000:26656000 "${CHROM_DIR}/Hap2_Chr06.fa"

seqkit stats *.fa

# ------------------------------------------------------------------------------
# STEP 3: BLASTn - Chr07 solitary unit vs 45S reference
# ------------------------------------------------------------------------------
echo "=========================================================================="
echo "[$(date)] STEP 3: BLASTn - Chr07 solitary rDNA unit vs 45S (ON685357.1)"
echo "=========================================================================="

cd "${BC_DIR}"
for Q in regions/Hap1_Chr07_rDNAunit.fa regions/Hap2_Chr07_rDNAunit.fa regions/Hap1_Chr06_dispersedpair.fa; do
    NAME=$(basename "${Q}" .fa)
    OUT="blast_results/${NAME}_vs_45S.tsv"
    [ -s "${OUT}" ] && continue
    blastn -query "${Q}" -db "${HPITS_DB}" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -max_target_seqs 5 | sort -k4,4 -rn > "${OUT}"
done

# ------------------------------------------------------------------------------
# STEP 4: BLASTn - Chr01/Chr06 5S loci vs nuclear+plastid 5S reference set
# ------------------------------------------------------------------------------
echo "=========================================================================="
echo "[$(date)] STEP 4: BLASTn - Chr01/Chr06 5S loci vs 5S references"
echo "=========================================================================="

for Q in regions/Hap1_Chr01_5S_main.fa regions/Hap1_Chr01_5S_distal.fa \
         regions/Hap2_Chr01_5S_main.fa regions/Hap2_Chr01_5S_distal.fa \
         regions/Hap1_Chr06_5S.fa regions/Hap2_Chr06_5S.fa; do
    NAME=$(basename "${Q}" .fa)
    OUT="blast_results/${NAME}_vs_5S.tsv"
    [ -s "${OUT}" ] && continue
    blastn -query "${Q}" -db "${REF5S_DIR}/Hp5S_db" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -max_target_seqs 4 -word_size 11 > "${OUT}"
done

# ------------------------------------------------------------------------------
# STEP 5: SUMMARY
# ------------------------------------------------------------------------------
echo "=========================================================================="
echo "[$(date)] STEP 5: SUMMARY"
echo "=========================================================================="

SUMMARY="${BC_DIR}/blast_confirmation_summary.tsv"
echo -e "query\tsubject\tn_hits\tmean_pident\tmin_pident\tmax_pident" > "${SUMMARY}"

for f in blast_results/*_vs_45S.tsv; do
    NAME=$(basename "${f}" _vs_45S.tsv)
    awk -v q="${NAME}" -F'\t' '{cnt[$2]++; sum[$2]+=$3; if($3>max[$2]||!(($2) in max))max[$2]=$3; if(!(($2) in min)||$3<min[$2])min[$2]=$3}
        END{for(s in cnt) printf "%s\t%s\t%d\t%.1f\t%.1f\t%.1f\n", q, s, cnt[s], sum[s]/cnt[s], min[s], max[s]}' "${f}" >> "${SUMMARY}"
done

for f in blast_results/*_vs_5S.tsv; do
    NAME=$(basename "${f}" _vs_5S.tsv)
    awk -v q="${NAME}" -F'\t' '{cnt[$2]++; sum[$2]+=$3; if($3>max[$2]||!(($2) in max))max[$2]=$3; if(!(($2) in min)||$3<min[$2])min[$2]=$3}
        END{for(s in cnt) printf "%s\t%s\t%d\t%.1f\t%.1f\t%.1f\n", q, s, cnt[s], sum[s]/cnt[s], min[s], max[s]}' "${f}" >> "${SUMMARY}"
done

echo "[+] Saved ${SUMMARY}"
column -t "${SUMMARY}"

echo "=========================================================================="
echo "[$(date)] BLAST CONFIRMATION COMPLETE"
echo "=========================================================================="
