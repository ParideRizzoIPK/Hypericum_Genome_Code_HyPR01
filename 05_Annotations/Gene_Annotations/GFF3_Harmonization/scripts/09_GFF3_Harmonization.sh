#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8                    
#SBATCH --mem=40G                            
#SBATCH --time=1-00:00:00                     
#SBATCH --job-name=gff_harmonize
#SBATCH --output=logs/gff_harmonize_%j.out
#SBATCH --error=logs/gff_harmonize_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  # Defer strict -u variable enforcement until profile initialization clears

# ==============================================================================
#                       ENVIRONMENT SETUP & PATHS
# ==============================================================================
echo "[$(date)] Setting up directories and validating entrypoints..."
source ~/.bashrc

set -u  # Activate strict variable enforcement safely after environment profile is loaded

BASE_OUT="/path/to/output/directory"
mkdir -p "${BASE_OUT}/Hap1" "${BASE_OUT}/Hap2"

HAPLOTYPES=("Hap1" "Hap2")

ANNEVO_GFF_H1="/path/to/input/directory/annevo_raw_hap1.gff3"
ANNEVO_GFF_H2="/path/to/input/directory/annevo_raw_hap2.gff3"

BRAKER_GFF_H1="/path/to/input/directory/braker.gff3"
BRAKER_GFF_H2="/path/to/input/directory/braker.gff3"

# Verify physical file footprints before initializing Python sub-shell
for file in "$ANNEVO_GFF_H1" "$ANNEVO_GFF_H2" "$BRAKER_GFF_H1" "$BRAKER_GFF_H2"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Input GFF3 asset missing: $file" >&2
        exit 1
    fi
done

# ==============================================================================
#             DYNAMIC WRITING OF HARMONIZATION ENGINE (PYTHON)
# ==============================================================================
echo "[$(date)] Writing inline Python harmonization engine with requested log updates..."
cat << 'EOF' > "${BASE_OUT}/harmonize_engine.py"
import argparse
import re
import sys
from collections import defaultdict

ID_RE = re.compile(r"(?:^|;)\s*ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)\s*Parent=([^;]+)")

def _attr(pattern, s):
    m = pattern.search(s)
    return m.group(1).strip() if m else None

def parse_gff(path):
    gene_order = []
    gene_info = {}
    gene_lines = defaultdict(list)
    gene_cds = defaultdict(list)
    tx_to_gene = {}
    current_gene = None

    with open(path) as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            ftype = cols[2]
            attrs = cols[8]
            fid = _attr(ID_RE, attrs)
            parent = _attr(PARENT_RE, attrs)

            if ftype == "gene":
                current_gene = fid
                gene_order.append(fid)
                gene_info[fid] = {
                    "seqid": cols[0], "source": cols[1],
                    "start": int(cols[3]), "end": int(cols[4]), "strand": cols[6],
                    "mrna_lines_count": 0
                }
                gene_lines[fid].append(raw.rstrip("\n"))
                continue

            if current_gene is not None:
                gene_lines[current_gene].append(raw.rstrip("\n"))

            if ftype == "mRNA":
                if fid is not None and parent is not None:
                    tx_to_gene[fid] = parent
                if current_gene in gene_info:
                    gene_info[current_gene]["mrna_lines_count"] += 1
            elif ftype == "CDS":
                gid = tx_to_gene.get(parent, current_gene)
                if gid is not None:
                    gene_cds[gid].append((int(cols[3]), int(cols[4])))

    return {"gene_order": gene_order, "gene_info": gene_info, "gene_lines": gene_lines, "gene_cds": gene_cds}

def merge_intervals(ivs):
    if not ivs: return []
    ivs = sorted(ivs)
    out = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]

def total_len(ivs):
    return sum(e - s + 1 for s, e in ivs)

def intersect_intervals(a, b):
    i = j = 0
    out = []
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0])
        e = min(a[i][1], b[j][1])
        if s <= e:
            out.append((s, e))
        if a[i][1] < b[j][1]: i += 1
        else: j += 1
    return out

def inject_attrs(raw_line, extra_pairs):
    cols = raw_line.split("\t")
    attrs = cols[8].rstrip().rstrip(";")
    extra = ";".join(f"{k}={v}" for k, v in extra_pairs)
    cols[8] = f"{attrs};{extra}" if attrs else extra
    return "\t".join(cols)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--braker", required=True)
    parser.add_argument("--annevo", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--hap", required=True)
    args = parser.parse_args()

    b_data = parse_gff(args.braker)
    a_data = parse_gff(args.annevo)

    b_m_ivs = {gid: merge_intervals(b_data["gene_cds"][gid]) for gid in b_data["gene_order"]}
    a_m_ivs = {gid: merge_intervals(a_data["gene_cds"][gid]) for gid in a_data["gene_order"]}

    BIN_SIZE = 100000
    a_bins = defaultdict(list)
    for ag in a_data["gene_order"]:
        info = a_data["gene_info"][ag]
        for s, e in a_m_ivs[ag]:
            for b in range(s // BIN_SIZE, e // BIN_SIZE + 1):
                a_bins[(info["seqid"], info["strand"], b)].append(ag)

    best_match = {}
    max_overlap_bp = defaultdict(int)

    for bg in b_data["gene_order"]:
        b_info = b_data["gene_info"][bg]
        b_intervals = b_m_ivs[bg]
        if not b_intervals: continue
        
        seen_partners = set()
        for s, e in b_intervals:
            for b in range(s // BIN_SIZE, e // BIN_SIZE + 1):
                for ag in a_bins[(b_info["seqid"], b_info["strand"], b)]:
                    if ag in seen_partners: continue
                    seen_partners.add(ag)
                    intersection = intersect_intervals(b_intervals, a_m_ivs[ag])
                    ov_bp = total_len(intersection)
                    if ov_bp > max_overlap_bp[bg]:
                        max_overlap_bp[bg] = ov_bp
                        best_match[bg] = ag

    absorbed_annevo = set()
    braker_cases = {}

    for bg in b_data["gene_order"]:
        blen = total_len(b_m_ivs[bg])
        ag = best_match.get(bg)

        if ag is None or blen == 0:
            braker_cases[bg] = ("Case 2", "High_confidence", None)
            continue

        alen = total_len(a_m_ivs[ag])
        overlap = max_overlap_bp[bg]

        if overlap == blen and alen == blen:
            braker_cases[bg] = ("Case 1", "Highest_confidence", ag)
            absorbed_annevo.add(ag)
        elif overlap == blen and alen > blen:
            ratio = blen / alen
            tier = "High_confidence" if ratio >= 0.70 else "Medium_confidence"
            case_lbl = "Case 4.1" if ratio >= 0.70 else "Case 4.2"
            braker_cases[bg] = (case_lbl, tier, ag)
            absorbed_annevo.add(ag)
        elif overlap > 0:
            braker_cases[bg] = ("Case 5", "High_confidence", ag)
            absorbed_annevo.add(ag)
        else:
            braker_cases[bg] = ("Case 2", "High_confidence", None)

    # Calculate final transcript metrics
    total_braker_mrnas = sum(b_data["gene_info"][bg]["mrna_lines_count"] for bg in b_data["gene_order"])
    total_annevo_mrnas_injected = sum(1 for bg in b_data["gene_order"] if braker_cases[bg][0] in ["Case 4.1", "Case 4.2", "Case 5"])
    case3_count = len(a_data["gene_order"]) - len(absorbed_annevo)
    total_final_mrnas = total_braker_mrnas + total_annevo_mrnas_injected + case3_count

    # Write detailed TSV census
    with open(args.report, "w") as rfh:
        rfh.write("gene_id\tsource\tdecision\tcase_class\tconfidence_tier\tpartner_id\n")
        for bg in b_data["gene_order"]:
            case_lbl, tier, ag = braker_cases[bg]
            rfh.write(f"{bg}\tbraker\tKEPT\t{case_lbl}\t{tier}\t{ag if ag else 'None'}\n")
        for ag in a_data["gene_order"]:
            if ag in absorbed_annevo:
                parent_bg = "None"
                for bg, (_, _, partner_ag) in braker_cases.items():
                    if partner_ag == ag:
                        parent_bg = bg
                        break
                rfh.write(f"{ag}\tannevo\tABSORBED\t-\t-\t{parent_bg}\n")
            else:
                rfh.write(f"{ag}\tannevo\tKEPT\tCase 3\tLow-Predicted\tNone\n")

    # Generate summary counts for log
    counts = defaultdict(int)
    for rec in braker_cases.values():
        counts[rec[0]] += 1
    counts["Case 3"] = case3_count

    # ---- REPORT 2: WRITE DETAILED SUMMARY RUN LOG ----
    with open(args.log, "w") as lfh:
        lfh.write("============================================================\n")
        lfh.write("GFF3 HARMONIZATION DESCRIPTIVE SUMMARY LOG\n")
        lfh.write("============================================================\n")
        lfh.write(f"Haplotype Focus Group      : {args.hap}\n")
        lfh.write("------------------------------------------------------------\n")
        lfh.write("HARMONIZATION CASE DEFINITIONS\n")
        lfh.write("  Reference GFF3 Baseline: BRAKER3\n")
        lfh.write("  Query GFF3 Dataset     : ANNEVO\n\n")
        lfh.write("  Case 1: ANNEVO predicted a CDS that covers 100% of the length of the\n")
        lfh.write("          BRAKER3 CDS and does not exceed it.\n")
        lfh.write("          Confidence class=Highest_confidence\n")
        lfh.write("          Decision: keep only the BRAKER3 CDS (absorbed)\n\n")
        lfh.write("  Case 2: ANNEVO did not predict any CDS, but BRAKER3 did.\n")
        lfh.write("          Confidence class= High_confidence\n")
        lfh.write("          Decision: Keep the BRAKER3 CDS\n\n")
        lfh.write("  Case 3: ANNEVO predicted a CDS, but BRAKER3 made no prediction.\n")
        lfh.write("          Confidence class=Low-Predicted\n")
        lfh.write("          Decision: Keep the ANNEVO CDS\n\n")
        lfh.write("  Case 4.1: ANNEVO predicted a CDS that covers 100% of the length of the\n")
        lfh.write("            BRAKER3 CDS, but it exceeds its length. BRAKER3 covers >= 70%.\n")
        lfh.write("            Confidence class=High_confidence\n")
        lfh.write("            Decision: Keep both gene models\n\n")
        lfh.write("  Case 4.2: ANNEVO predicted a CDS that covers 100% of the length of the\n")
        lfh.write("            BRAKER3 CDS, but it exceeds its length. BRAKER3 covers < 70%.\n")
        lfh.write("            Confidence class=Medium_confidence\n")
        lfh.write("            Decision: Keep both gene models\n\n")
        lfh.write("  Case 5: ANNEVO predicted a CDS that partially covers the length of the\n")
        lfh.write("          BRAKER3 CDS because it is shorter.\n")
        lfh.write("          Confidence class=High_confidence\n")
        lfh.write("          Decision: Keep both gene models\n")
        lfh.write("------------------------------------------------------------\n")
        lfh.write("QUANTITATIVE METRICS & CURATION TALLIES\n\n")
        
        b_total = len(b_data["gene_order"])
        lfh.write(f"Total Input BRAKER3 Genes  : {b_total}\n")
        for c in ["Case 1", "Case 2", "Case 4.1", "Case 4.2", "Case 5"]:
            cnt = counts[c]
            pct = (cnt / b_total * 100) if b_total > 0 else 0.0
            lfh.write(f"  -> {c:<10}: {cnt:>6}   ({pct:.2f}%)\n")
            
        a_total = len(a_data["gene_order"])
        lfh.write(f"\nTotal Input ANNEVO Genes   : {a_total}\n")
        c3_cnt = counts["Case 3"]
        c3_pct = (c3_cnt / a_total * 100) if a_total > 0 else 0.0
        lfh.write(f"  -> Case 3 (Kept Novel)   : {c3_cnt:>6}   ({c3_pct:.2f}%)\n")
        
        abs_cnt = len(absorbed_annevo)
        abs_pct = (abs_cnt / a_total * 100) if a_total > 0 else 0.0
        lfh.write(f"  -> Absorbed & Discarded  : {abs_cnt:>6}   ({abs_pct:.2f}%)\n")
        lfh.write("------------------------------------------------------------\n")
        lfh.write("CONSOLIDATED LOCI SUMMARY\n")
        lfh.write(f"  (i)  Number of Different Genes (Loci)   : {b_total + c3_cnt}\n")
        lfh.write(f"  (ii) Number of Gene Models (Transcripts): {total_final_mrnas}\n")
        lfh.write("============================================================\n")

    # Emit output GFF lines
    final_blocks = []
    for bg in b_data["gene_order"]:
        case_lbl, tier, ag = braker_cases[bg]
        b_info = b_data["gene_info"][bg]
        gene_extra = [("confidence_tier", tier), ("harmonization_case", case_lbl)]
        if ag: gene_extra.append(("annevo_partner", ag))

        lines = []
        for i, raw in enumerate(b_data["gene_lines"][bg]):
            extra = gene_extra if i == 0 else [("confidence_tier", tier)]
            lines.append(inject_attrs(raw, extra))

        if ag and case_lbl in ["Case 4.1", "Case 4.2", "Case 5"]:
            for raw in a_data["gene_lines"][ag]:
                cols = raw.split("\t")
                if cols[2] == "gene": continue
                attrs = cols[8]
                fid = _attr(ID_RE, attrs)
                parent = _attr(PARENT_RE, attrs)

                if cols[2] == "mRNA":
                    new_tx_id = f"{bg}_annevo_{fid}"
                    raw = raw.replace(f"ID={fid}", f"ID={new_tx_id}")
                    raw = re.sub(r"Parent=[^;]+", f"Parent={bg}", raw)
                else:
                    if parent:
                        new_parent_id = f"{bg}_annevo_{parent}"
                        raw = re.sub(r"Parent=[^;]+", f"Parent={new_parent_id}", raw)
                lines.append(inject_attrs(raw, [("confidence_tier", tier), ("isoform_source", "ANNEVO")]))
        final_blocks.append((b_info["seqid"], b_info["start"], b_info["end"], bg, lines))

    for ag in a_data["gene_order"]:
        if ag in absorbed_annevo: continue
        a_info = a_data["gene_info"][ag]
        lines = []
        for i, raw in enumerate(a_data["gene_lines"][ag]):
            extra = [("confidence_tier", "Low-Predicted"), ("harmonization_case", "Case 3")] if i == 0 else [("confidence_tier", "Low-Predicted")]
            lines.append(inject_attrs(raw, extra))
        final_blocks.append((a_info["seqid"], a_info["start"], a_info["end"], ag, lines))

    final_blocks.sort(key=lambda x: (x[0], x[1], x[2]))

    with open(args.out, "w") as out_fh:
        out_fh.write("##gff-version 3\n")
        for _, _, _, _, lines in final_blocks:
            for ln in lines:
                out_fh.write(ln + "\n")
            out_fh.write("###\n")

if __name__ == "__main__":
    main()
EOF
chmod +x "${BASE_OUT}/harmonize_engine.py"

# ==============================================================================
#                          RUN HARMONIZATION LOOP
# ==============================================================================
for HAP in "${HAPLOTYPES[@]}"; do
    echo "------------------------------------------------------------"
    echo "[$(date)] Launching Harmonization Engine for ${HAP}..."
    echo "------------------------------------------------------------"
    
    if [ "$HAP" == "Hap1" ]; then
        IN_ANNEVO="$ANNEVO_GFF_H1"
        IN_BRAKER="$BRAKER_GFF_H1"
    else
        IN_ANNEVO="$ANNEVO_GFF_H2"
        IN_BRAKER="$BRAKER_GFF_H2"
    fi
    
    OUT_GFF="${BASE_OUT}/${HAP}/harmonized_consensus_${HAP}.gff3"
    OUT_TSV="${BASE_OUT}/${HAP}/gff_harmonization_report.tsv"
    OUT_LOG="${BASE_OUT}/${HAP}/gff_harmonization_run.log"
    
    python3 "${BASE_OUT}/harmonize_engine.py" \
        --braker="$IN_BRAKER" \
        --annevo="$IN_ANNEVO" \
        --out="$OUT_GFF" \
        --report="$OUT_TSV" \
        --log="$OUT_LOG" \
        --hap="$HAP"
        
    echo "[$(date)] ${HAP} curation completely compiled."
done

rm -f "${BASE_OUT}/harmonize_engine.py"
echo "[$(date)] Comprehensive pipeline execution finished cleanly."