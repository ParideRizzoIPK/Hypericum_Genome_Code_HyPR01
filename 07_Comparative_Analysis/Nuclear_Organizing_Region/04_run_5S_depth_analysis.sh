#!/bin/bash
#
# ==============================================================================
# 5S rDNA Locus Copy-Number Estimation from Existing mosdepth Depth
# ==============================================================================
# Follow-up from genomewide_repeat_scan/HyPR01_Genomewide_Repeat_Landscape_Report.md
# section 10/12: pull the existing 10kb-window mosdepth depth (from the
# competitive HiFi coverage analysis, 20260725_HyPR01_Genome_Coverage/) for the
# Chr01 and Chr06 5S rDNA loci coordinates identified in this project's
# genome-wide repeat scan, and combine it with the assembled tandem-copy counts
# (from barrnap-HGV, already in summary/rRNA_gene_summary.tsv) to estimate
# absolute haploid copy number - the same method used for the Chr08 NOR in the
# companion Chr08-specific report.
#
# DATA HYGIENE: the source mosdepth bed.gz is only ever read (gzip -dc to
# stdout); it is never modified, and no file is ever written back into
# 20260725_HyPR01_Genome_Coverage/. All output goes to
# genomewide_repeat_scan/depth_analysis/.
# ==============================================================================

set -euo pipefail
trap 'echo "[FATAL] line ${LINENO} exited with status $?" >&2' ERR

BASE_DIR="/path/to/your/directory"
SCAN_DIR="${BASE_DIR}/genomewide_repeat_scan"
OUT_DIR="${SCAN_DIR}/depth_analysis"

COV_DIR="/path/to/input/directory"
MOSDEPTH_BED_GZ="${COV_DIR}/HyPR01_diploid_competitive_mosdepth.regions.bed.gz"
BASELINE=23.24   # haploid nuclear median depth, from HyPR01_Coverage_Pipeline_Report.md

mkdir -p "${OUT_DIR}"

BED="${OUT_DIR}/mosdepth_10kb_windows.bed"
if [ ! -s "${BED}" ]; then
    echo "[$(date)] Decompressing mosdepth regions bed (read-only from source, never written back)..."
    gzip -dc "${MOSDEPTH_BED_GZ}" > "${BED}"
fi

mean_depth () {
    local chrom="$1" start="$2" end="$3"
    awk -v c="${chrom}" -v s="${start}" -v e="${end}" -F'\t' '
        $1==c && $2<e && $3>s {sum+=$4; n++}
        END{ if(n>0) printf "%.2f\t%d", sum/n, n; else print "NA\t0" }' "${BED}"
}

# Locus definitions: name, chrom, start, end, assembled_copy_count
# (coordinates and copy counts from summary/rRNA_gene_summary.tsv / the
# genome-wide report's cluster breakdown, section 4)
LOCI_FILE=$(mktemp)
cat > "${LOCI_FILE}" << 'LOCI'
Chr01_Hap1_main	Chr01_Hap1	21871837	22070205	508
Chr01_Hap1_satellite	Chr01_Hap1	20163558	22164296	209
Chr01_Hap1_distal	Chr01_Hap1	48623438	48654430	89
Chr01_Hap2_main	Chr01_Hap2	21365094	21447987	196
Chr01_Hap2_satellite	Chr01_Hap2	22952535	23184493	307
Chr01_Hap2_distal	Chr01_Hap2	50200593	50270756	151
Chr06_Hap1_A	Chr06_Hap1	25290595	25380187	255
Chr06_Hap1_B	Chr06_Hap1	26008713	26189527	476
Chr06_Hap2_A	Chr06_Hap2	26375032	26453918	223
Chr06_Hap2_B	Chr06_Hap2	26479879	26518233	110
Chr06_Hap2_C	Chr06_Hap2	26542599	26645888	294
LOCI

echo "=========================================================================="
echo "[$(date)] Computing mean depth per locus (baseline = ${BASELINE}x haploid nuclear)"
echo "=========================================================================="

SUMMARY="${OUT_DIR}/5S_locus_depth_summary.tsv"
echo -e "locus\tchrom\tstart\tend\tn_windows\tmean_depth\tfold_vs_baseline\tassembled_copies\timplied_true_copies" > "${SUMMARY}"

while IFS=$'\t' read -r name chrom start end copies; do
    result=$(mean_depth "${chrom}" "${start}" "${end}")
    depth=$(echo "${result}" | cut -f1)
    nwin=$(echo "${result}" | cut -f2)
    if [ "${depth}" = "NA" ]; then
        echo -e "${name}\t${chrom}\t${start}\t${end}\t0\tNA\tNA\t${copies}\tNA" >> "${SUMMARY}"
        continue
    fi
    fold=$(python3 -c "print(f'{${depth}/${BASELINE}:.2f}')")
    implied=$(python3 -c "print(round(${copies}*${fold}))")
    echo -e "${name}\t${chrom}\t${start}\t${end}\t${nwin}\t${depth}\t${fold}\t${copies}\t${implied}" >> "${SUMMARY}"
done < "${LOCI_FILE}"

rm -f "${LOCI_FILE}"

echo "[+] Saved ${SUMMARY}"
column -t "${SUMMARY}"

echo "=========================================================================="
echo "[$(date)] 5S DEPTH ANALYSIS COMPLETE"
echo "=========================================================================="
