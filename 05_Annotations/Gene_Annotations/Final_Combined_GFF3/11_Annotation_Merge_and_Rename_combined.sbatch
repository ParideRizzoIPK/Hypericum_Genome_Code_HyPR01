#!/bin/bash
#SBATCH --job-name=HyPR01_Full_Annotate_PurgeChrUn_EMBL
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/Annotate_Reindex_%j.out
#SBATCH --error=logs/Annotate_Reindex_%j.err

# ──────────────────────────────────────────────────────────────────────────────
# HyPR01 Workflow Execution Order:
#   1. Collapse redundant transcript models (structural dedup on full genome).
#   2. Remove transcripts with premature internal stops OR missing start 'M'.
#   3. Merge eggNOG-mapper + InterProScan into GFF3 column 9.
#   4. Rename/reindex locus tags (01G..08G, UnG), position-ordered and strand-aware.
#   5. Recompute each CDS phase to match exact exon coordinates.
#   6. Post-hoc verification of zero internal stops and complete 'M' starts.
#   7. Cross-haplotype casing canonicalization (retaining full genome models).
#   8. Full-genome feature-count summary (including ChrUn).
#   9. POST-STEP 8 PURGING & SUBMISSION PREP:
#      a. Generate ChrUn-purged chromosome FASTAs (hap1/hap2.masked_chromosomes.fa).
#      b. Purge all ChrUn scaffolds, genes, transcripts, and child features.
#      c. Export detailed audit log TSVs tracking every purged feature.
#      d. Consolidate custom attributes into pipe-delimited Note= qualifiers.
#      e. Compute side-by-side comparative feature-count statistics (Full vs Purged vs Retained).
#  10. EMBL conversion: run EMBLmyGFF3 on the submission files to produce .embl.gz flatfiles.
# ──────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ══════════════════════════════════════════════════════════════════════
#  HEADER METADATA & CONFIGURATION
# ══════════════════════════════════════════════════════════════════════
ORGANISM="Hypericum perforatum"
ASSEMBLY="HyPR01"
SPECIES_CODE="Hp"
ASSEMBLY_ACCESSION=""
PROVIDER="Leibniz Institute of Plant Genetics and Crop Plant Research (IPK Gatersleben)"
AUTHORS="Paride Rizzo, Benjamin Chavez"
DOI=""
PIPELINE_VERSION="1.0.0"
PROJECT_ID="PRJEB123953"
TAX_DIV="PLN"
SUBMISSION_EMAIL="rizzo@ipk-gatersleben.de"
STUDY_TITLE="A chromosome-level, haplotype-resolved genome of a diploid sexual Hypericum perforatum L. producing hypericin"
EMBLMYGFF3_BIN="/path/to/your/directory/micromamba/envs/emblmygff3/bin/EMBLmyGFF3"
# ══════════════════════════════════════════════════════════════════════

[[ "${ASSEMBLY}"     =~ ^[A-Za-z0-9]+$ ]] || { echo "[ERROR] ASSEMBLY must be alphanumeric only: '${ASSEMBLY}'" >&2; exit 1; }
[[ "${SPECIES_CODE}" =~ ^[A-Za-z0-9]+$ ]] || { echo "[ERROR] SPECIES_CODE must be alphanumeric only: '${SPECIES_CODE}'" >&2; exit 1; }

ENABLE_STRUCTURAL_DEDUP="true"
ENABLE_INTERNAL_STOP_FILTER="true"

GENOME_FA_HAP1="/path/to/input/directory/hap1.masked_ChrUn.fa"
GENOME_FA_HAP2="/path/to/input/directory/hap2.masked_ChrUn.fa"

FAI_HAP1="${GENOME_FA_HAP1}.fai"
FAI_HAP2="${GENOME_FA_HAP2}.fai"

# 1. Load System Modules & Environments
. /etc/profile.d/modules.sh
module load python/3.12.0

BASE="/path/to/your/directory"
EGG_BASE="${BASE}/HyPR01_Genome_EGGNOG"
IPR_BASE="${BASE}/HyPR01_Genome_InterproScan"
FINAL_DIR="${BASE}/HyPR01_Genome_Final_Annotations"

mkdir -p "${FINAL_DIR}/HyPR01_Hap1" "${FINAL_DIR}/HyPR01_Hap2"

# ──────────────────────────────────────────────────────────────────────────────
# 2. Write Helper Scripts (Steps 1, 2, 3-5, 7)
# ──────────────────────────────────────────────────────────────────────────────

# 2a. Step 1: structural transcript dedup
FILTER_PY="collapse_structural_duplicates.py"
cat << 'PYEOF1' > "$FILTER_PY"
#!/usr/bin/env python3
import sys, argparse
from collections import defaultdict, OrderedDict

def parse_attrs(s):
    d = OrderedDict()
    for item in s.strip().split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            d[k.strip()] = v.strip()
    return d

TX_TYPES = ('mRNA', 'transcript', 'tRNA', 'rRNA', 'ncRNA')
UTR_TYPES = ('five_prime_UTR', 'three_prime_UTR')

def merge_blocks(intervals):
    if not intervals: return ()
    ordered = sorted(intervals)
    out = [list(ordered[0])]
    for s, e in ordered[1:]:
        if s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return tuple((s, e) for s, e in out)

def group_by(items, keyfn):
    buckets = defaultdict(list)
    for it in items: buckets[keyfn(it)].append(it)
    return buckets

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('gff_in')
    ap.add_argument('--report', required=True)
    ap.add_argument('--log', default='')
    ap.add_argument('--write', default='')
    ap.add_argument('--renumber', action='store_true')
    a = ap.parse_args()

    lines, gene_order = [], []
    gene_tx = defaultdict(list)
    tx_gene, tx_locus, tx_order = {}, {}, {}
    tx_cds, tx_utr, tx_exon = defaultdict(list), defaultdict(list), defaultdict(list)

    with open(a.gff_in) as f:
        for i, line in enumerate(f):
            if line.startswith('#') or not line.strip():
                lines.append((line, None, None)); continue
            p = line.rstrip('\n').split('\t')
            if len(p) != 9:
                lines.append((line, None, None)); continue
            attrs = parse_attrs(p[8])
            ftype = p[2]
            lines.append((line, ftype, attrs))
            if ftype == 'gene':
                gid = attrs.get('ID')
                gene_order.append(gid)
                gene_tx.setdefault(gid, [])
            elif ftype in TX_TYPES:
                tid = attrs.get('ID'); gid = attrs.get('Parent')
                tx_gene[tid] = gid
                gene_tx[gid].append(tid)
                tx_locus[tid] = (p[0], p[6])
                tx_order[tid] = i
            elif ftype == 'CDS':
                tx_cds[attrs.get('Parent')].append((int(p[3]), int(p[4])))
            elif ftype in UTR_TYPES:
                tx_utr[attrs.get('Parent')].append((int(p[3]), int(p[4])))
            elif ftype == 'exon':
                tx_exon[attrs.get('Parent')].append((int(p[3]), int(p[4])))

    def full_key(t):
        exon_chain = merge_blocks(tx_exon[t]) or merge_blocks(tx_cds[t] + tx_utr[t])
        return (tx_locus[t], exon_chain, tuple(sorted(tx_cds[t])))

    def cds_key(t): return (tx_locus[t], tuple(sorted(tx_cds[t])))
    def has_utr(t): return bool(tx_utr[t])

    remove, remove_reason = {}, {}
    for gid in gene_tx:
        txs = gene_tx[gid]
        if len(txs) < 2: continue
        structured = [t for t in txs if tx_exon[t] or tx_cds[t] or tx_utr[t]]
        for group in group_by(structured, full_key).values():
            if len(group) < 2: continue
            group.sort(key=lambda t: tx_order[t])
            for t in group[1:]:
                remove[t] = group[0]
                remove_reason[t] = 'exact_full_exon_twin'

        for group in group_by(txs, cds_key).values():
            group = [t for t in group if tx_cds[t] and t not in remove]
            if len(group) < 2: continue
            withutr = [t for t in group if has_utr(t)]
            without = [t for t in group if not has_utr(t)]
            if not withutr or not without: continue
            keep = min(withutr, key=lambda t: tx_order[t])
            for t in without:
                remove[t] = keep
                remove_reason[t] = 'utrless_vs_utr_bearing'

    for t in list(remove):
        seen = {t}
        dest = remove[t]
        while dest in remove and dest not in seen:
            seen.add(dest)
            dest = remove[dest]
        remove[t] = dest

    for gid, txs in gene_tx.items():
        survivors = [t for t in txs if t not in remove]
        if not survivors and txs:
            for t in txs:
                remove.pop(t, None)
                remove_reason.pop(t, None)

    def chain_str(blocks):
        return ";".join(f"{s}-{e}" for s, e in blocks) if blocks else "."

    with open(a.report, 'w') as r:
        r.write("gene\tremoved_transcript\tkept_transcript\treason\tn_CDS\tremoved_CDS_chain\tremoved_exon_chain\tkept_exon_chain\tremoved_has_UTR\tkept_has_UTR\n")
        for t in sorted(remove, key=lambda x: (tx_gene[x], x)):
            k = remove[t]
            cds = tuple(sorted(tx_cds[t]))
            r.write(f"{tx_gene[t]}\t{t}\t{k}\t{remove_reason[t]}\t{len(cds)}\t{chain_str(cds)}\t{chain_str(full_key(t)[1])}\t{chain_str(full_key(k)[1])}\t{'yes' if has_utr(t) else 'no'}\t{'yes' if has_utr(k) else 'no'}\n")

    if not a.write: return

    rename = {}
    if a.renumber:
        for gid, txs in gene_tx.items():
            survivors = [t for t in txs if t not in remove]
            survivors.sort(key=lambda t: tx_order[t])
            for idx, t in enumerate(survivors, 1):
                rename[t] = f"{gid}.t{idx}"

    with open(a.write, 'w') as out:
        for raw, ftype, attrs in lines:
            if ftype is None: out.write(raw); continue
            if ftype in TX_TYPES:
                tid = attrs.get('ID')
                if tid in remove: continue
                if a.renumber and tid in rename:
                    p = raw.rstrip('\n').split('\t')
                    attrs['ID'] = rename[tid]
                    if 'Name' in attrs: attrs['Name'] = rename[tid]
                    p[8] = ";".join(f"{k}={v}" for k, v in attrs.items())
                    out.write("\t".join(p) + "\n"); continue
            elif ftype not in ('gene',):
                parent = attrs.get('Parent')
                if parent in remove: continue
                if a.renumber and parent in rename:
                    p = raw.rstrip('\n').split('\t')
                    old_child = attrs.get('ID', '')
                    attrs['Parent'] = rename[parent]
                    if old_child.startswith(parent + '.'):
                        attrs['ID'] = rename[parent] + old_child[len(parent):]
                    p[8] = ";".join(f"{k}={v}" for k, v in attrs.items())
                    out.write("\t".join(p) + "\n"); continue
            out.write(raw)

if __name__ == '__main__':
    main()
PYEOF1

# 2b. Step 2: internal-stop-codon & canonical start 'M' filter
STOPFILTER_PY="filter_internal_stop_codons.py"
cat << 'PYEOF1B' > "$STOPFILTER_PY"
#!/usr/bin/env python3
import sys, argparse, time
from collections import defaultdict, OrderedDict

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T', 'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
_CODON_CACHE = {}

def resolve_codon(codon):
    if codon in _CODON_CACHE: return _CODON_CACHE[codon]
    if codon in CODON_TABLE:
        _CODON_CACHE[codon] = CODON_TABLE[codon]
        return CODON_TABLE[codon]
    try:
        opts0, opts1, opts2 = IUPAC[codon[0]], IUPAC[codon[1]], IUPAC[codon[2]]
    except (KeyError, IndexError):
        _CODON_CACHE[codon] = 'X'
        return 'X'
    aas = {CODON_TABLE[b0 + b1 + b2] for b0 in opts0 for b1 in opts1 for b2 in opts2}
    result = aas.pop() if len(aas) == 1 else 'X'
    _CODON_CACHE[codon] = result
    return result

def revcomp(s): return s.translate(COMP)[::-1]

def parse_attrs(s):
    d = OrderedDict()
    for item in s.strip().split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            d[k.strip()] = v.strip()
    return d

TX_TYPES = ('mRNA', 'transcript')

def load_genome(path):
    seqs = {}
    name, chunks = None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None: seqs[name] = ''.join(chunks)
                name = line[1:].split()[0].strip()
                chunks = []
            else: chunks.append(line.strip())
        if name is not None: seqs[name] = ''.join(chunks)
    return seqs

def translate(seq):
    n_codons = len(seq) // 3
    aas = [resolve_codon(seq[i * 3:i * 3 + 3].upper()) for i in range(n_codons)]
    n_internal, out = 0, []
    for i, aa in enumerate(aas):
        is_last = (i == len(aas) - 1)
        if aa == '*':
            if is_last: continue
            out.append('.')
            n_internal += 1
        else: out.append(aa)
    return ''.join(out), n_internal

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('gff_in'); ap.add_argument('genome_fa')
    ap.add_argument('--report', required=True); ap.add_argument('--log', required=True)
    ap.add_argument('--proteins', required=True); ap.add_argument('--stops-fa', required=True)
    ap.add_argument('--primary-fa', default=''); ap.add_argument('--write', default='')
    a = ap.parse_args()

    run_date = time.strftime('%Y-%m-%d %H:%M:%S %Z')
    genome = load_genome(a.genome_fa)
    lines, gene_order = [], []
    gene_tx = defaultdict(list)
    tx_gene, tx_seqid, tx_strand = {}, {}, {}
    tx_cds = defaultdict(list)
    tx_canonical = set()

    with open(a.gff_in) as f:
        for i, line in enumerate(f):
            if line.startswith('#') or not line.strip():
                lines.append((line, None, None)); continue
            p = line.rstrip('\n').split('\t')
            if len(p) != 9:
                lines.append((line, None, None)); continue
            attrs = parse_attrs(p[8])
            ftype = p[2]
            lines.append((line, ftype, attrs))
            if ftype == 'gene':
                gid = attrs.get('ID')
                gene_order.append(gid)
                gene_tx.setdefault(gid, [])
            elif ftype in TX_TYPES:
                tid = attrs.get('ID'); gid = attrs.get('Parent')
                tx_gene[tid] = gid
                gene_tx[gid].append(tid)
                tx_seqid[tid] = p[0]
                tx_strand[tid] = p[6]
                if attrs.get('canonical', '').lower() == 'true': tx_canonical.add(tid)
            elif ftype == 'CDS':
                tid = attrs.get('Parent')
                tx_cds[tid].append((int(p[3]), int(p[4]), p[7]))

    proteins, rejected_tx = {}, {}
    n_no_seq = 0

    for tid, segs in tx_cds.items():
        seqid = tx_seqid.get(tid)
        strand = tx_strand.get(tid)
        if seqid not in genome or not segs:
            n_no_seq += 1
            continue
        segs_sorted = sorted(segs, key=lambda x: x[0])
        chrom_seq = genome[seqid]
        spliced = ''.join(chrom_seq[s - 1:e] for s, e, _ in segs_sorted)
        if strand == '-':
            spliced = revcomp(spliced)
            anchor_phase_str = segs_sorted[-1][2]
        else:
            anchor_phase_str = segs_sorted[0][2]
        try: anchor_phase = int(anchor_phase_str)
        except ValueError: anchor_phase = 0
        
        protein, n_internal = translate(spliced[anchor_phase:])
        proteins[tid] = protein
        
        starts_with_m = protein.startswith('M') if protein else False
        has_internal_stop = (n_internal > 0)
        
        if has_internal_stop or not starts_with_m:
            if has_internal_stop and not starts_with_m: reason = "both"
            elif has_internal_stop: reason = "internal_stop"
            else: reason = "missing_start_M"
            rejected_tx[tid] = (reason, n_internal, starts_with_m, protein)

    with open(a.report, 'w') as r:
        r.write("gene\tremoved_transcript\treason\tn_CDS\tcds_len_nt\tn_internal_stops\tstarts_with_M\n")
        for tid in sorted(rejected_tx, key=lambda x: (tx_gene[x], x)):
            reason, n_internal, starts_with_m, _ = rejected_tx[tid]
            cds_len = sum(e - s + 1 for s, e, _ in tx_cds[tid])
            r.write(f"{tx_gene[tid]}\t{tid}\t{reason}\t{len(tx_cds[tid])}\t{cds_len}\t{n_internal}\t{'yes' if starts_with_m else 'no'}\n")

    genes_dropped = set()
    for gid, txs in gene_tx.items():
        survivors = [t for t in txs if t not in rejected_tx]
        if txs and not survivors:
            genes_dropped.add(gid)

    with open(a.log, 'w') as lg:
        lg.write(f"Scanned: {len(proteins)+n_no_seq} | Translated: {len(proteins)} | Rejected: {len(rejected_tx)} | Dropped genes: {len(genes_dropped)}\n")

    with open(a.proteins, 'w') as pf:
        for tid in sorted(proteins):
            pf.write(f">{tid}\n")
            seq = proteins[tid]
            for i in range(0, len(seq), 70): pf.write(seq[i:i + 70] + "\n")

    with open(a.stops_fa, 'w') as sf:
        for tid in sorted(rejected_tx):
            _, _, _, protein = rejected_tx[tid]
            sf.write(f">{tid}\n")
            for i in range(0, len(protein), 70): sf.write(protein[i:i + 70] + "\n")

    if a.primary_fa:
        primary = sorted(t for t in tx_canonical if t in proteins and t not in rejected_tx)
        with open(a.primary_fa, 'w') as pf:
            for tid in primary:
                pf.write(f">{tid} gene={tx_gene.get(tid, '.')}\n")
                seq = proteins[tid]
                for i in range(0, len(seq), 70): pf.write(seq[i:i + 70] + "\n")

    if not a.write: return

    remove_tx = set(rejected_tx)
    remove_gene = genes_dropped
    with open(a.write, 'w') as out:
        for raw, ftype, attrs in lines:
            if ftype is None: out.write(raw); continue
            if ftype == 'gene':
                if attrs.get('ID') in remove_gene: continue
                out.write(raw); continue
            if ftype in TX_TYPES:
                if attrs.get('ID') in remove_tx: continue
                out.write(raw); continue
            parent = attrs.get('Parent')
            if parent in remove_tx: continue
            out.write(raw)

if __name__ == '__main__':
    main()
PYEOF1B

# 2c. Steps 3-5: functional merge, rename, phase-fix (maintains UnG tags for ChrUn)
ANNOT_PY="annotate_and_reindex_combined.py"
cat << 'PYEOF2' > "$ANNOT_PY"
#!/usr/bin/env python3
import sys, os, re, csv, time, argparse
from collections import defaultdict

GO_ROOTS = {"GO:0003674", "GO:0008150", "GO:0005575"}
_ESCAPE_LIST = [('%', '%25'), ('\t', '%09'), ('\n', '%0A'), ('\r', '%0D'), (';', '%3B'), ('=', '%3D'), ('&', '%26'), (',', '%2C')]
_ESCAPE_TEXT = [('%', '%25'), ('\t', '%09'), ('\n', '%0A'), ('\r', '%0D'), (';', '%3B'), ('=', '%3D'), ('&', '%26')]
_SHORTNAME = re.compile(r'^[A-Za-z][A-Za-z0-9.]*(_[A-Za-z0-9.]+)+$')
_ACCESSIONISH = re.compile(r'^(PTHR|PF|SM|SSF|G3DSA|cd|PS|TIGR)[0-9]', re.IGNORECASE)
_LOCUS_TAG = re.compile(r'\(?(?:TAIR\s+)?\b(?:Os\d{2}[Gg]\d{7}|At\dg\d{5}|LOC\d{5,})\b(?:\.\d+)?\)?', re.IGNORECASE)
_BEST_MATCH_PREFIX = re.compile(r'^BEST\s+[A-Z][a-z]+\s+[a-z]+\s+protein\s+match\s+is\s*', re.IGNORECASE)
_DANGLING_TRAILING_WORD = re.compile(r'(?i)\b(of|the|and|or|with|in|to|for|by|from|on|at|as|is)\s*$')
_DANGLING_TAIR_REF = re.compile(r'\(TAIR\s*$', re.IGNORECASE)
_DUP_RELATED_SUFFIX = re.compile(r'(-related){2,}\s*$', re.IGNORECASE)
_PUTATIVE_RELATED_COMBO = re.compile(r'putative.*-related|-related.*putative', re.IGNORECASE)
_GENERIC_ONLY = {'protein', 'expressed protein', 'domain-containing protein', 'family protein', 'protein family'}
_BELONGS_TO_FAMILY = re.compile(r'^belongs to the .+ (family|superfamily)$', re.IGNORECASE)
_ISOFORM_ONLY = re.compile(r'^isoform x\d+$', re.IGNORECASE)
_ENCODED_BY_ONLY = re.compile(r'^encoded by$', re.IGNORECASE)
_GENERIC_SUPERFAMILY_DESCRIPTIONS = {'p-loop containing nucleoside triphosphate hydrolases', 's-adenosyl-l-methionine-dependent methyltransferases', 'transferase family', 'alpha/beta hydrolase family'}
_COMMON_SHORT_WORDS = {'a', 'an', 'the', 'and', 'or', 'nor', 'but', 'if', 'so', 'of', 'in', 'on', 'at', 'to', 'for', 'with', 'from', 'by', 'as', 'is', 'are', 'was', 'were', 'be', 'been', 'being', 'that', 'this', 'these', 'those', 'it', 'its', 'into', 'onto', 'over', 'also', 'very', 'some', 'most', 'many', 'each', 'both', 'when', 'where', 'what', 'which', 'not', 'than', 'then', 'no', 'yes', 'can', 'may', 'will', 'shall', 'do', 'does', 'did', 'has', 'have', 'had', 'you', 'your', 'like', 'such', 'via', 'per', 'non', 'type', 'part', 'more', 'less', 'and/or'}

def gff3_escape_list(val):
    if not val: return ""
    for ch, esc in _ESCAPE_LIST: val = val.replace(ch, esc)
    return val

def gff3_escape_text(val):
    if not val: return ""
    for ch, esc in _ESCAPE_TEXT: val = val.replace(ch, esc)
    return val

def clean_tag(val):
    if not val or val.strip() in ("", "-", "None"): return ""
    return val.strip()

def is_shortname(desc):
    if not desc: return True
    return ' ' not in desc and bool(_SHORTNAME.match(desc) or _ACCESSIONISH.match(desc))

def is_low_value(desc):
    if not desc: return True
    d = desc.strip().lower()
    if 'unknown function' in d or d.startswith('duf') or 'uncharacteri' in d or 'hypothetical' in d or d in _GENERIC_ONLY or d in _GENERIC_SUPERFAMILY_DESCRIPTIONS or _BELONGS_TO_FAMILY.match(d) or _ISOFORM_ONLY.match(d) or _ENCODED_BY_ONLY.match(d) or _PUTATIVE_RELATED_COMBO.search(d):
        return True
    return False

def clean_description(desc):
    if not desc: return desc
    desc = _DUP_RELATED_SUFFIX.sub('-related', desc.strip())
    m = _BEST_MATCH_PREFIX.match(desc)
    if m:
        desc = desc[m.end():].strip()
        if not desc or is_low_value(desc) or _DANGLING_TRAILING_WORD.search(desc): return ""
    if _LOCUS_TAG.search(desc):
        desc = _LOCUS_TAG.sub('', desc)
        desc = re.sub(r'\s{2,}', ' ', desc)
        desc = re.sub(r'\s+([,;])', r'\1', desc).strip(' ,;-')
        if not desc or desc.strip().lower() in _GENERIC_ONLY: return ""
    if _DANGLING_TAIR_REF.search(desc):
        desc = _DANGLING_TAIR_REF.sub('', desc).strip(' ,;-')
        if not desc or desc.strip().lower() in _GENERIC_ONLY: return ""
    return desc

def concise(desc):
    if not desc: return desc
    parts = [p.strip().rstrip(';').rstrip('.').strip() for p in re.split(r'(?<=[.;])\s+', desc.strip()) if p.strip()]
    if not parts: return desc
    first = parts[0]
    if len(parts) > 1 and is_low_value(first):
        rest = [p for p in parts[1:] if p and not is_low_value(p)]
        if rest: return ", ".join([first] + rest)
    if 15 <= len(first) < len(desc): return first
    return desc

def smart_case(desc):
    letters = [c for c in desc if c.isalpha()]
    if not (letters and all(c.isupper() for c in letters) and len(desc) > 3): return desc
    def repl(m):
        w = m.group(0)
        return w if len(w) <= 4 and w.lower() not in _COMMON_SHORT_WORDS else w.capitalize()
    return re.sub(r'[A-Za-z]+', repl, desc)

def parse_eggnog(path):
    data = {}
    if not os.path.exists(path): return data
    with open(path) as f:
        for row in csv.reader(f, delimiter='\t'):
            if not row or row[0].startswith('#'): continue
            try: data[row[0]] = {'desc': clean_tag(row[7]), 'gene_name': clean_tag(row[8]), 'go': clean_tag(row[9]), 'ec': clean_tag(row[10])}
            except IndexError: continue
    return data

def parse_interpro(path):
    data = defaultdict(lambda: {'dbxrefs': set(), 'go': set(), 'descs': []})
    if not os.path.exists(path): return data
    seen = defaultdict(set)
    with open(path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 12: continue
            q = parts[0]
            db, acc, desc, ipr = clean_tag(parts[3]), clean_tag(parts[4]), clean_tag(parts[5]), clean_tag(parts[11])
            if db and acc: data[q]['dbxrefs'].add(f"{db}:{acc}")
            if ipr and ipr != '-': data[q]['dbxrefs'].add(f"InterPro:{ipr}")
            if desc and not is_shortname(desc):
                desc = smart_case(clean_description(desc))
                if desc and desc.lower() not in seen[q]:
                    seen[q].add(desc.lower())
                    data[q]['descs'].append(desc)
            if len(parts) >= 14 and clean_tag(parts[13]):
                for go in parts[13].split('|'):
                    go = clean_tag(go)
                    if go: data[q]['go'].add(go)
    return data

def choose_product(egg, ipr_descs):
    egg_desc = smart_case(clean_description(egg['desc'])) if egg and egg['desc'] else ''
    if egg_desc and not is_low_value(egg_desc): return egg_desc
    good_ipr = [d for d in ipr_descs if not is_low_value(d)]
    if good_ipr: return max(good_ipr, key=len)
    return "hypothetical protein"

def parse_attributes(attr_str):
    d = {}
    for item in attr_str.strip().split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            d[k.strip()] = v.strip()
    return d

def format_attributes(attrs): return ';'.join(f"{k}={v}" for k, v in attrs.items())

def get_chrom_code(chrom_str):
    if "un" in chrom_str.lower(): return "Un"
    m = re.search(r'\d+', chrom_str)
    return f"{int(m.group()):02d}" if m else "Un"

_OWNED = ['Note', 'eggnog_note', 'interpro_note', 'product', 'gene_name', 'EC_number', 'Ontology_term_evidence', 'pfam', 'panther']

def pick_canonical(transcripts, transcript_children):
    def rank(item):
        idx, tr = item
        children = transcript_children.get(tr['attrs'].get('ID', ''), [])
        cds_bp = sum(c['end'] - c['start'] + 1 for c in children if c['type'] == 'CDS')
        exon_bp = sum(c['end'] - c['start'] + 1 for c in children if c['type'] == 'exon')
        utr_sides = len({c['type'] for c in children if c['type'].endswith('_UTR')})
        return (-cds_bp, -exon_bp, -utr_sides, idx)
    if not transcripts: return None
    return min(enumerate(transcripts), key=rank)[0]

def annotate(attrs, tid, egg_map, ipr_map):
    for k in _OWNED: attrs.pop(k, None)
    egg, ipr = egg_map.get(tid), ipr_map.get(tid)
    if not (egg or ipr): return False

    ipr_descs = ipr['descs'] if ipr else []
    product_full = choose_product(egg, ipr_descs)
    prod = gff3_escape_text(concise(product_full))
    attrs['product'] = prod

    egg_desc_clean = smart_case(clean_description(egg['desc'])) if egg and egg['desc'] else ''
    if egg_desc_clean: attrs['eggnog_note'] = gff3_escape_text(egg_desc_clean)

    ipr_ordered = [d for d in ipr_descs if d and not is_shortname(d)]
    if ipr_ordered: attrs['interpro_note'] = ",".join(gff3_escape_list(d) for d in ipr_ordered)

    dbx = set(x for x in attrs.get('Dbxref', '').split(',') if x)
    if ipr: dbx.update(ipr['dbxrefs'])
    if dbx:
        attrs['Dbxref'] = ",".join(sorted(dbx))
        pfam_ids = sorted(x.split(':', 1)[1] for x in dbx if x.startswith('Pfam:'))
        if pfam_ids: attrs['pfam'] = ",".join(pfam_ids)
        panther_ids = sorted(x.split(':', 1)[1] for x in dbx if x.startswith('PANTHER:'))
        if panther_ids: attrs['panther'] = ",".join(panther_ids)

    go = set(x.strip() for x in attrs.get('Ontology_term', '').split(',') if x.strip())
    if egg and egg['go']: go.update(x.strip() for x in egg['go'].split(','))
    if ipr: go.update(ipr['go'])
    go = {g for g in go if g and g != '-' and g not in GO_ROOTS}
    if go:
        attrs['Ontology_term'] = ",".join(sorted(go))
        attrs['Ontology_term_evidence'] = 'IEA'

    if egg and egg['gene_name']: attrs['gene_name'] = gff3_escape_text(egg['gene_name'])
    if egg and egg['ec']: attrs['EC_number'] = gff3_escape_text(egg['ec'])
    return True

def build_header(meta_headers, chrom_max_end, fai_lengths, meta):
    lines = ["##gff-version 3\n"]
    for seqid in sorted(chrom_max_end):
        length = fai_lengths[seqid] if fai_lengths and seqid in fai_lengths else chrom_max_end[seqid]
        lines.append(f"##sequence-region {seqid} 1 {length}\n")
    lines.append(f"#!genome-build {meta['assembly']}\n")
    if meta['assembly_accession']: lines.append(f"#!genome-build-accession {meta['assembly_accession']}\n")
    lines.append(f"#!genotype {meta['assembly']}\n")
    lines.append(f"#!organism {meta['organism']}\n")
    lines.append(f"#!haplotype {meta['haplotype']}\n")
    lines.append(f"#!study-accession {meta['project_id']}\n")
    lines.append(f"#!annotation-provider {meta['provider']}\n")
    lines.append(f"#!annotation-date {meta['date']}\n")
    lines.append(f"#!annotation-pipeline {meta['pipeline']}\n")
    lines.append(f"# Functional evidence: {meta['evidence']}\n")
    if meta['authors']: lines.append(f"# Contributors: {meta['authors']}\n")
    if meta['doi']: lines.append(f"# Citation/DOI: {meta['doi']}\n")
    return lines

def load_fai(path):
    lengths = {}
    if not path or not os.path.exists(path): return lengths
    with open(path) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 2 and p[1].isdigit(): lengths[p[0]] = int(p[1])
    return lengths

def compute_phase_corrections(cds_by_mrna):
    corrections = []
    for mrna, records in cds_by_mrna.items():
        strand = records[0]['strand']
        ordered = sorted(records, key=lambda r: r['start'], reverse=(strand == '-'))
        if any(r['phase'] is None for r in ordered): continue
        anchor = ordered[0]['phase']
        adjusted_cum = anchor
        for i, rec in enumerate(ordered):
            length = rec['end'] - rec['start'] + 1
            if i > 0:
                expected = (3 - (adjusted_cum % 3)) % 3
                if rec['phase'] != expected:
                    corrections.append({'body_idx': rec['body_idx'], 'new_phase': expected})
            adjusted_cum += length
    return corrections

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('gff_in'); ap.add_argument('eggnog'); ap.add_argument('interpro'); ap.add_argument('gff_out'); ap.add_argument('mapping_tsv'); ap.add_argument('audit_log'); ap.add_argument('prefix')
    ap.add_argument('--fai', default=''); ap.add_argument('--organism', default=''); ap.add_argument('--assembly', default='HyPR01'); ap.add_argument('--assembly-accession', default=''); ap.add_argument('--haplotype', default=''); ap.add_argument('--provider', default=''); ap.add_argument('--authors', default=''); ap.add_argument('--doi', default=''); ap.add_argument('--pipeline-version', default='1.0.0'); ap.add_argument('--project-id', default='')
    a = ap.parse_args()

    run_date = time.strftime('%Y-%m-%d %H:%M:%S %Z')
    meta = {'organism': a.organism, 'assembly': a.assembly, 'assembly_accession': a.assembly_accession, 'haplotype': a.haplotype or a.prefix, 'provider': a.provider, 'authors': a.authors, 'doi': a.doi, 'date': run_date, 'pipeline': f"combined-annotate-reindex-phasefix v{a.pipeline_version}", 'evidence': "eggNOG-mapper + InterProScan", 'project_id': a.project_id}
    egg_map, ipr_map, fai_lengths = parse_eggnog(a.eggnog), parse_interpro(a.interpro), load_fai(a.fai)

    chrom_genes, gene_transcripts, transcript_children, chrom_max_end = defaultdict(list), defaultdict(list), defaultdict(list), defaultdict(int)
    meta_headers, body_lines, cds_by_mrna = [], [], defaultdict(list)

    with open(a.gff_in) as f:
        for line in f:
            if line.startswith('#'): continue
            if not line.strip(): continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 9: continue
            chrom, source, ftype, start, end, score, strand, phase, attr_str = parts
            end_i = int(end)
            if end_i > chrom_max_end[chrom]: chrom_max_end[chrom] = end_i
            fd = {'chrom': chrom, 'source': source, 'type': ftype, 'start': int(start), 'end': end_i, 'score': score, 'strand': strand, 'phase': phase, 'attrs': parse_attributes(attr_str)}
            if ftype == 'gene': chrom_genes[chrom].append(fd); gene_transcripts.setdefault(fd['attrs'].get('ID'), [])
            elif ftype in ('mRNA', 'transcript', 'tRNA', 'rRNA'):
                if fd['attrs'].get('Parent'): gene_transcripts[fd['attrs']['Parent']].append(fd)
                transcript_children.setdefault(fd['attrs'].get('ID'), [])
            else:
                if fd['attrs'].get('Parent'): transcript_children[fd['attrs']['Parent']].append(fd)

    with open(a.mapping_tsv, 'w') as tsv_out:
        tsv_out.write("Feature_Type\tOld_Identifier\tNew_Identifier\tChrom_Column1\tStart_Coord\tEnd_Coord\n")
        header_lines = build_header(meta_headers, chrom_max_end, fai_lengths, meta)
        for chrom in sorted(chrom_genes.keys()):
            chrom_code = get_chrom_code(chrom)
            genes = sorted(chrom_genes[chrom], key=lambda x: x['start'])
            for g_idx, gene in enumerate(genes):
                old_gene_id = gene['attrs'].get('ID', 'UNKNOWN')
                new_gene_id = f"{a.prefix}_{chrom_code}G{(g_idx + 1) * 10:06d}"
                gene['attrs']['ID'] = new_gene_id
                if 'Name' in gene['attrs']: gene['attrs']['Name'] = new_gene_id

                transcripts_sorted = sorted(gene_transcripts.get(old_gene_id, []), key=lambda x: x['start'])
                for tr in transcripts_sorted: annotate(tr['attrs'], tr['attrs'].get('ID', ''), egg_map, ipr_map)

                canon_idx = pick_canonical(transcripts_sorted, transcript_children)
                if canon_idx is not None: gene['attrs']['canonical_transcript'] = f"{new_gene_id}.t{canon_idx + 1}"
                rep = (transcripts_sorted[canon_idx]['attrs'].get('product') if canon_idx is not None else None)
                gene['attrs']['description'] = rep if rep else "hypothetical protein"

                body_lines.append("\t".join([gene['chrom'], gene['source'], gene['type'], str(gene['start']), str(gene['end']), gene['score'], gene['strand'], gene['phase'], format_attributes(gene['attrs'])]) + "\n")
                tsv_out.write(f"gene\t{old_gene_id}\t{new_gene_id}\t{gene['chrom']}\t{gene['start']}\t{gene['end']}\n")

                for t_idx, tr in enumerate(transcripts_sorted):
                    old_tx = tr['attrs'].get('ID', 'UNKNOWN')
                    new_tx = f"{new_gene_id}.t{t_idx + 1}"
                    tr['attrs']['ID'] = new_tx
                    tr['attrs']['Parent'] = new_gene_id
                    if 'Name' in tr['attrs']: tr['attrs']['Name'] = new_tx
                    if t_idx == canon_idx: tr['attrs']['canonical'] = 'true'
                    body_lines.append("\t".join([tr['chrom'], tr['source'], tr['type'], str(tr['start']), str(tr['end']), tr['score'], tr['strand'], tr['phase'], format_attributes(tr['attrs'])]) + "\n")
                    tsv_out.write(f"transcript\t{old_tx}\t{new_tx}\t{tr['chrom']}\t{tr['start']}\t{tr['end']}\n")

                    children = transcript_children.get(old_tx, [])
                    by_type = defaultdict(list)
                    for c in children: by_type[c['type']].append(c)
                    for c_type, tc in by_type.items():
                        rev = (tr['strand'] == '-')
                        for idx, c in enumerate(sorted(tc, key=lambda x: x['start'], reverse=rev)):
                            c['attrs']['ID'] = f"{new_tx}.{c_type}{idx + 1}"
                            c['attrs']['Parent'] = new_tx
                    for c in sorted(children, key=lambda x: x['start']):
                        body_lines.append("\t".join([c['chrom'], c['source'], c['type'], str(c['start']), str(c['end']), c['score'], c['strand'], c['phase'], format_attributes(c['attrs'])]) + "\n")
                        if c['type'] == 'CDS':
                            try: phase_int = int(c['phase'])
                            except ValueError: phase_int = None
                            cds_by_mrna[new_tx].append({'start': c['start'], 'end': c['end'], 'strand': c['strand'], 'phase': phase_int, 'body_idx': len(body_lines) - 1})

    corrections = compute_phase_corrections(cds_by_mrna)
    for c in corrections:
        parts = body_lines[c['body_idx']].rstrip('\n').split('\t')
        parts[7] = str(c['new_phase'])
        body_lines[c['body_idx']] = '\t'.join(parts) + '\n'

    with open(a.gff_out, 'w') as gff_out:
        gff_out.writelines(header_lines)
        gff_out.writelines(body_lines)

if __name__ == '__main__':
    main()
PYEOF2

# 2d. Step 7: cross-haplotype casing canonicalization
CANON_PY="canonicalize_casing.py"
cat << 'PYEOF3' > "$CANON_PY"
#!/usr/bin/env python3
import argparse, re
from collections import Counter, defaultdict

CURATED_OVERRIDES = {
    'cytochrome p450': 'Cytochrome P450', 'f-box protein': 'F-box protein', 'zinc finger protein': 'Zinc finger protein',
    'transcription factor': 'Transcription factor', 'ring finger protein': 'RING finger protein', 'wd40 repeat protein': 'WD40 repeat protein',
    'wd40 repeat-containing protein': 'WD40 repeat-containing protein', 'leucine-rich repeat protein': 'Leucine-rich repeat protein',
    'leucine-rich repeat receptor-like protein kinase': 'Leucine-rich repeat receptor-like protein kinase', 'nac domain-containing protein': 'NAC domain-containing protein',
    'myb domain-containing protein': 'MYB domain-containing protein', 'myb-related protein': 'MYB-related protein', 'myb family transcription factor': 'MYB family transcription factor',
    'bzip transcription factor': 'bZIP transcription factor', 'ap2/erf domain-containing protein': 'AP2/ERF domain-containing protein', 'wrky transcription factor': 'WRKY transcription factor',
    'gdsl esterase/lipase': 'GDSL esterase/lipase', 'f-box/lrr-repeat protein': 'F-box/LRR-repeat protein', 'atp-binding cassette transporter': 'ATP-binding cassette transporter',
    'dna-binding protein': 'DNA-binding protein', 'rna-binding protein': 'RNA-binding protein', 'dna repair protein': 'DNA repair protein',
    'atp synthase subunit': 'ATP synthase subunit', 'gtp-binding protein': 'GTP-binding protein', 'nadh dehydrogenase': 'NADH dehydrogenase',
    'nad(p)h dehydrogenase': 'NAD(P)H dehydrogenase', 'nb-arc domain-containing protein': 'NB-ARC domain-containing protein', 'tir-nbs-lrr protein': 'TIR-NBS-LRR protein',
}

_DESC_RE = re.compile(r'description=([^;]*)')
_PROD_RE = re.compile(r'product=([^;]*)')
_EGGNOG_NOTE_RE = re.compile(r'eggnog_note=([^;]*)')
_INTERPRO_NOTE_RE = re.compile(r'interpro_note=([^;]*)')

def load_lines(path):
    with open(path) as f: return f.readlines()

def collect_vocabulary(all_lines_by_file, counter):
    for lines in all_lines_by_file:
        for line in lines:
            if line.startswith('#') or not line.strip(): continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9: continue
            ftype, attrs = p[2], p[8]
            if ftype == 'gene':
                m = _DESC_RE.search(attrs)
                if m and m.group(1): counter[m.group(1).strip().lower()][m.group(1).strip()] += 1
            elif ftype == 'mRNA':
                m = _PROD_RE.search(attrs)
                if m and m.group(1): counter[m.group(1).strip().lower()][m.group(1).strip()] += 1
                m = _EGGNOG_NOTE_RE.search(attrs)
                if m and m.group(1): counter[m.group(1).strip().lower()][m.group(1).strip()] += 1
                m = _INTERPRO_NOTE_RE.search(attrs)
                if m and m.group(1):
                    for item in m.group(1).split(','):
                        if item.strip(): counter[item.strip().lower()][item.strip()] += 1

def build_canonical_map(counter):
    canon = {}
    for key, variants in counter.items():
        if len(variants) <= 1: continue
        if key in CURATED_OVERRIDES: canon[key] = CURATED_OVERRIDES[key]
        else:
            max_count = max(variants.values())
            tied = sorted(v for v, c in variants.items() if c == max_count)
            canon[key] = tied[0]
    return canon

def apply_canonical_map(lines, canon):
    out = []
    for line in lines:
        if line.startswith('#') or not line.strip(): out.append(line); continue
        p = line.rstrip('\n').split('\t')
        if len(p) < 9: out.append(line); continue
        ftype, attrs = p[2], p[8]
        changed = False
        if ftype == 'gene':
            m = _DESC_RE.search(attrs)
            if m and m.group(1):
                new = canon.get(m.group(1).strip().lower())
                if new and new != m.group(1).strip(): attrs = attrs[:m.start(1)] + new + attrs[m.end(1):]; changed = True
        elif ftype == 'mRNA':
            m = _PROD_RE.search(attrs)
            if m and m.group(1):
                new = canon.get(m.group(1).strip().lower())
                if new and new != m.group(1).strip(): attrs = attrs[:m.start(1)] + new + attrs[m.end(1):]; changed = True
            m = _EGGNOG_NOTE_RE.search(attrs)
            if m and m.group(1):
                new = canon.get(m.group(1).strip().lower())
                if new and new != m.group(1).strip(): attrs = attrs[:m.start(1)] + new + attrs[m.end(1):]; changed = True
            m = _INTERPRO_NOTE_RE.search(attrs)
            if m and m.group(1):
                new_items = [canon.get(i.strip().lower(), i.strip()) for i in m.group(1).split(',') if i.strip()]
                new_note = ",".join(new_items)
                if new_note != m.group(1): attrs = attrs[:m.start(1)] + new_note + attrs[m.end(1):]; changed = True
        if changed: p[8] = attrs; out.append('\t'.join(p) + '\n')
        else: out.append(line)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('gff1'); ap.add_argument('gff2'); ap.add_argument('--out1', required=True); ap.add_argument('--out2', required=True); ap.add_argument('--log', default='')
    a = ap.parse_args()
    lines1, lines2 = load_lines(a.gff1), load_lines(a.gff2)
    counter = defaultdict(Counter)
    collect_vocabulary([lines1, lines2], counter)
    canon = build_canonical_map(counter)
    out1 = apply_canonical_map(lines1, canon)
    out2 = apply_canonical_map(lines2, canon)
    with open(a.out1, 'w') as f: f.writelines(out1)
    with open(a.out2, 'w') as f: f.writelines(out2)

if __name__ == '__main__':
    main()
PYEOF3

# ──────────────────────────────────────────────────────────────────────────────
# 3. Execution: Steps 1-5 Per Haplotype (Full Genome with ChrUn)
# ──────────────────────────────────────────────────────────────────────────────
echo "[$(date)] Launching Steps 1-5 on full genomes (including ChrUn)..."

for HAP in Hap1 Hap2; do
    if [[ "${HAP}" == "Hap1" ]]; then 
        HAP_CODE="H1"
        FAI="${FAI_HAP1}"
        GENOME_FA="${GENOME_FA_HAP1}"
    else 
        HAP_CODE="H2"
        FAI="${FAI_HAP2}"
        GENOME_FA="${GENOME_FA_HAP2}"
    fi
    PREFIX="${SPECIES_CODE}_${ASSEMBLY}_${HAP_CODE}"
    echo "[$(date)] ===== Processing ${HAP} (prefix ${PREFIX}) ====="

    GFF_IN="${EGG_BASE}/HyPR01_${HAP}/HyPR01_${HAP}_final_annotations.gff3"
    EGG_TSV="${EGG_BASE}/HyPR01_${HAP}/HyPR01_${HAP}_emapper.emapper.annotations"
    IPR_TSV="${IPR_BASE}/HyPR01_${HAP}/HyPR01_${HAP}_interproscan.tsv"
    OUT_DIR="${FINAL_DIR}/HyPR01_${HAP}"
    mkdir -p "${OUT_DIR}"

    if [[ ! -f "${GFF_IN}" ]]; then
        echo "[ERROR] Missing input GFF3: ${GFF_IN}. Skipping." >&2
        continue
    fi

    ANNOT_INPUT="${GFF_IN}"
    DEDUP_TMP=""
    STOPFILTER_TMP=""

    if [[ "${ENABLE_STRUCTURAL_DEDUP}" == "true" ]]; then
        echo "[$(date)] ${HAP}: collapsing redundant transcript models..."
        DEDUP_TMP="${OUT_DIR}/.HyPR01_${HAP}_structural_deduplicated.tmp.gff3"
        python3 "${FILTER_PY}" "${GFF_IN}" \
            --report "${OUT_DIR}/HyPR01_${HAP}_structural_dedup_report.tsv" \
            --log    "${OUT_DIR}/HyPR01_${HAP}_structural_dedup.log" \
            --write  "${DEDUP_TMP}"
        ANNOT_INPUT="${DEDUP_TMP}"
    fi

    if [[ "${ENABLE_INTERNAL_STOP_FILTER}" == "true" ]]; then
        echo "[$(date)] ${HAP}: filtering transcripts with internal stop codons or missing start 'M'..."
        STOPFILTER_TMP="${OUT_DIR}/.HyPR01_${HAP}_stopfiltered.tmp.gff3"
        python3 "${STOPFILTER_PY}" "${ANNOT_INPUT}" "${GENOME_FA}" \
            --report "${OUT_DIR}/HyPR01_${HAP}_internal_stop_removal_report.tsv" \
            --log    "${OUT_DIR}/HyPR01_${HAP}_internal_stop_removal.log" \
            --proteins "${OUT_DIR}/HyPR01_${HAP}_proteins.fa" \
            --stops-fa "${OUT_DIR}/HyPR01_${HAP}_proteins_internal_stops.fa" \
            --write  "${STOPFILTER_TMP}"
        if [[ -n "${DEDUP_TMP}" ]]; then rm -f "${DEDUP_TMP}"; fi
        ANNOT_INPUT="${STOPFILTER_TMP}"
    fi

    echo "[$(date)] ${HAP}: annotating + reindexing + correcting CDS phase..."
    python3 "${ANNOT_PY}" \
        "${ANNOT_INPUT}" "${EGG_TSV}" "${IPR_TSV}" \
        "${OUT_DIR}/HyPR01_${HAP}_final_annotations_functional.gff3" \
        "${OUT_DIR}/HyPR01_${HAP}_gene_mapping.tsv" \
        "${OUT_DIR}/HyPR01_${HAP}_annotation_audit.log" \
        "${PREFIX}" \
        --organism "${ORGANISM}" \
        --assembly "${ASSEMBLY}" \
        --assembly-accession "${ASSEMBLY_ACCESSION}" \
        --haplotype "${HAP}" \
        --provider "${PROVIDER}" \
        --authors "${AUTHORS}" \
        --doi "${DOI}" \
        --pipeline-version "${PIPELINE_VERSION}" \
        --project-id "${PROJECT_ID}" \
        --fai "${FAI}"

    if [[ -n "${STOPFILTER_TMP}" ]]; then rm -f "${STOPFILTER_TMP}"; fi
    echo "[$(date)] ${HAP}: initial full-genome processing complete."
done

# ──────────────────────────────────────────────────────────────────────────────
# 4. Step 7: Cross-Haplotype Casing Canonicalization (Full Genome)
# ──────────────────────────────────────────────────────────────────────────────
HAP1_FINAL="${FINAL_DIR}/HyPR01_Hap1/HyPR01_Hap1_final_annotations_functional.gff3"
HAP2_FINAL="${FINAL_DIR}/HyPR01_Hap2/HyPR01_Hap2_final_annotations_functional.gff3"

if [[ -f "${HAP1_FINAL}" && -f "${HAP2_FINAL}" ]]; then
    echo "[$(date)] Canonicalizing product and description casing across both haplotypes..."
    HAP1_CANON_TMP="${FINAL_DIR}/HyPR01_Hap1/.HyPR01_Hap1_canon.tmp.gff3"
    HAP2_CANON_TMP="${FINAL_DIR}/HyPR01_Hap2/.HyPR01_Hap2_canon.tmp.gff3"
    python3 "${CANON_PY}" "${HAP1_FINAL}" "${HAP2_FINAL}" \
        --out1 "${HAP1_CANON_TMP}" \
        --out2 "${HAP2_CANON_TMP}" \
        --log  "${FINAL_DIR}/HyPR01_casing_canonicalization.log"
    mv "${HAP1_CANON_TMP}" "${HAP1_FINAL}"
    mv "${HAP2_CANON_TMP}" "${HAP2_FINAL}"
    echo "[$(date)] Casing canonicalization done."
fi

# ──────────────────────────────────────────────────────────────────────────────
# 5. Step 6: Post-Hoc Verification & Representative Protein FASTA Generation
# ──────────────────────────────────────────────────────────────────────────────
echo "[$(date)] Verifying transcripts and generating full & representative protein files..."
for HAP in Hap1 Hap2; do
    if [[ "${HAP}" == "Hap1" ]]; then GENOME_FA="${GENOME_FA_HAP1}"; else GENOME_FA="${GENOME_FA_HAP2}"; fi
    HAP_FINAL="${FINAL_DIR}/HyPR01_${HAP}/HyPR01_${HAP}_final_annotations_functional.gff3"
    VERIFY_REPORT="${FINAL_DIR}/HyPR01_${HAP}/HyPR01_${HAP}_FINAL_stop_codon_verification_report.tsv"
    
    python3 "${STOPFILTER_PY}" "${HAP_FINAL}" "${GENOME_FA}" \
        --report "${VERIFY_REPORT}" \
        --log    "${FINAL_DIR}/HyPR01_${HAP}/HyPR01_${HAP}_FINAL_stop_codon_verification.log" \
        --proteins "${FINAL_DIR}/HyPR01_${HAP}/HyPR01_${HAP}_FINAL_proteins.fa" \
        --stops-fa "${FINAL_DIR}/HyPR01_${HAP}/HyPR01_${HAP}_FINAL_proteins_internal_stops.fa" \
        --primary-fa "${FINAL_DIR}/HyPR01_${HAP}/HyPR01_${HAP}_FINAL_proteins_representative.fa"
        
    N_REMAINING=$(( $(wc -l < "${VERIFY_REPORT}") - 1 ))
    if [[ "${N_REMAINING}" -gt 0 ]]; then
        echo "[ERROR] ${HAP}: ${N_REMAINING} transcript(s) fail validation. See ${VERIFY_REPORT}." >&2
        exit 1
    fi
    echo "[$(date)] ${HAP}: verified -- all transcripts start with 'M' and contain zero internal stops."
done

# ──────────────────────────────────────────────────────────────────────────────
# 6. Step 8: Compute Feature-Count Summary (Full Genome with ChrUn)
# ──────────────────────────────────────────────────────────────────────────────
SUMMARY_FILE="${FINAL_DIR}/HyPR01_feature_counts_summary.txt"
count_features() {
    awk -F'\t' '
        $3=="gene"{g++} $3=="mRNA"{m++}
        $3=="CDS"{c++; n = split($9, parts, ";"); for (i=1; i<=n; i++) { if (parts[i] ~ /^Parent=/) { parent = substr(parts[i], 8); if (!(parent in seen)) { seen[parent] = 1; cdstx++ }; break } } }
        $3=="exon"{e++}
        END{printf "%d\t%d\t%d\t%d\t%d\n", g+0, m+0, c+0, cdstx+0, e+0}
    ' "$1"
}

count_described() {
    awk -F'\t' '
        $3=="gene" { gt++; i = index($9, "description="); if (i > 0) { rest = substr($9, i + length("description=")); j = index(rest, ";"); val = (j > 0) ? substr(rest, 1, j-1) : rest; if (tolower(val) != "hypothetical protein") gd++ } }
        $3=="mRNA" { mt++; i = index($9, "product="); if (i > 0) { rest = substr($9, i + length("product=")); j = index(rest, ";"); val = (j > 0) ? substr(rest, 1, j-1) : rest; if (tolower(val) != "hypothetical protein") md++ } }
        END{printf "%d\t%d\t%d\t%d\n", gd+0, gt+0, md+0, mt+0}
    ' "$1"
}

count_utr_lengths() {
    python3 - "$@" << 'PYEOF_UTR'
import sys
from collections import defaultdict
utr5, utr3 = defaultdict(int), defaultdict(int)
for path in sys.argv[1:]:
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            p = line.rstrip('\n').split('\t')
            if len(p) != 9 or p[2] not in ('five_prime_UTR', 'three_prime_UTR'): continue
            parent = next((kv.split('=', 1)[1] for kv in p[8].split(';') if kv.startswith('Parent=')), None)
            if parent is None: continue
            (utr5 if p[2] == 'five_prime_UTR' else utr3)[parent] += (int(p[4]) - int(p[3]) + 1)
def stats(d):
    vals = sorted(v for v in d.values() if v > 0)
    if not vals: return (0, 0.0, 0.0)
    n = len(vals)
    return (n, sum(vals)/n, vals[n//2] if n%2 else (vals[n//2-1]+vals[n//2])/2)
n5, mean5, median5 = stats(utr5); n3, mean3, median3 = stats(utr3)
print(f"{n5}\t{mean5:.1f}\t{median5:.1f}\t{n3}\t{mean3:.1f}\t{median3:.1f}")
PYEOF_UTR
}

{
    echo "======================================================================"
    echo "FULL GENOME FEATURE COUNTS (INCLUDING CHRUN)"
    echo "======================================================================"
    echo "Run date : $(date)"
    echo ""
    TOTAL_GENES=0; TOTAL_MRNA=0; TOTAL_CDSSEG=0; TOTAL_CDSTX=0; TOTAL_EXON=0
    TOTAL_GENES_DESC=0; TOTAL_MRNA_DESC=0
    EXISTING_GFFS=()
    for HAP in Hap1 Hap2; do
        HAP_GFF="${FINAL_DIR}/HyPR01_${HAP}/HyPR01_${HAP}_final_annotations_functional.gff3"
        if [[ -f "${HAP_GFF}" ]]; then
            EXISTING_GFFS+=("${HAP_GFF}")
            read -r G M CSEG CTX E <<< "$(count_features "${HAP_GFF}")"
            read -r GD GT MD MT <<< "$(count_described "${HAP_GFF}")"
            read -r N5 MEAN5 MED5 N3 MEAN3 MED3 <<< "$(count_utr_lengths "${HAP_GFF}")"
            echo "${HAP}:"
            echo "    Genes : ${G}"
            echo "    mRNA  : ${M}"
            echo "    Complete CDS (unique coding transcripts) : ${CTX}"
            echo "    CDS segments (coding exons)             : ${CSEG}"
            echo "    Exon                                    : ${E}"
            echo "    Genes with description                  : ${GD} / ${GT} ($(awk -v a="$GD" -v b="$GT" 'BEGIN{printf "%.1f", (b>0)?100*a/b:0}')%)"
            echo "    mRNA  with product                      : ${MD} / ${MT} ($(awk -v a="$MD" -v b="$MT" 'BEGIN{printf "%.1f", (b>0)?100*a/b:0}')%)"
            echo "    mRNA with 5' UTR : ${N5} / ${M} ($(awk -v a="$N5" -v b="$M" 'BEGIN{printf "%.1f", (b>0)?100*a/b:0}')%)   mean=${MEAN5} bp   median=${MED5} bp"
            echo "    mRNA with 3' UTR : ${N3} / ${M} ($(awk -v a="$N3" -v b="$M" 'BEGIN{printf "%.1f", (b>0)?100*a/b:0}')%)   mean=${MEAN3} bp   median=${MED3} bp"
            echo ""
            TOTAL_GENES=$((TOTAL_GENES + G)); TOTAL_MRNA=$((TOTAL_MRNA + M))
            TOTAL_CDSSEG=$((TOTAL_CDSSEG + CSEG)); TOTAL_CDSTX=$((TOTAL_CDSTX + CTX)); TOTAL_EXON=$((TOTAL_EXON + E))
            TOTAL_GENES_DESC=$((TOTAL_GENES_DESC + GD)); TOTAL_MRNA_DESC=$((TOTAL_MRNA_DESC + MD))
        fi
    done
    read -r TOTAL_N5 TOTAL_MEAN5 TOTAL_MED5 TOTAL_N3 TOTAL_MEAN3 TOTAL_MED3 <<< "$(count_utr_lengths "${EXISTING_GFFS[@]}")"
    echo "Combined (Hap1 + Hap2):"
    echo "    Genes : ${TOTAL_GENES}"
    echo "    mRNA  : ${TOTAL_MRNA}"
    echo "    Complete CDS (unique coding transcripts) : ${TOTAL_CDSTX}"
    echo "    CDS segments (coding exons)             : ${TOTAL_CDSSEG}"
    echo "    Exon                                    : ${TOTAL_EXON}"
    echo "    Genes with description                  : ${TOTAL_GENES_DESC} / ${TOTAL_GENES} ($(awk -v a="$TOTAL_GENES_DESC" -v b="$TOTAL_GENES" 'BEGIN{printf "%.1f", (b>0)?100*a/b:0}')%)"
    echo "    mRNA  with product                      : ${TOTAL_MRNA_DESC} / ${TOTAL_MRNA} ($(awk -v a="$TOTAL_MRNA_DESC" -v b="$TOTAL_MRNA" 'BEGIN{printf "%.1f", (b>0)?100*a/b:0}')%)"
    echo "    mRNA with 5' UTR : ${TOTAL_N5} / ${TOTAL_MRNA} ($(awk -v a="$TOTAL_N5" -v b="$TOTAL_MRNA" 'BEGIN{printf "%.1f", (b>0)?100*a/b:0}')%)   mean=${TOTAL_MEAN5} bp   median=${TOTAL_MED5} bp"
    echo "    mRNA with 3' UTR : ${TOTAL_N3} / ${TOTAL_MRNA} ($(awk -v a="$TOTAL_N3" -v b="$TOTAL_MRNA" 'BEGIN{printf "%.1f", (b>0)?100*a/b:0}')%)   mean=${TOTAL_MEAN3} bp   median=${TOTAL_MED3} bp"
    echo "======================================================================"
} | tee "${SUMMARY_FILE}"


# ──────────────────────────────────────────────────────────────────────────────
# 7. Step 9: Post-Step 8 ChrUn Purging, Submission GFF3 Prep & Detailed Audit
# ──────────────────────────────────────────────────────────────────────────────
echo "[$(date)] Step 9: Creating chromosome-only FASTAs, purging ChrUn, creating audit logs & generating submission GFF3s..."

PURGE_AUDIT_PY="purge_chrun_and_prep_submission.py"
cat << 'PYEOF4' > "$PURGE_AUDIT_PY"
#!/usr/bin/env python3
import sys, os, re
from collections import defaultdict, OrderedDict

def parse_fasta(fa_path):
    seqs = OrderedDict()
    curr_name = None
    curr_seq = []
    with open(fa_path) as f:
        for line in f:
            if line.startswith('>'):
                if curr_name: seqs[curr_name] = "".join(curr_seq)
                curr_name = line[1:].split()[0].strip()
                curr_seq = []
            else: curr_seq.append(line.strip())
        if curr_name: seqs[curr_name] = "".join(curr_seq)
    return seqs

def write_fasta(seqs, out_fa):
    with open(out_fa, 'w') as out:
        for hdr, seq in seqs.items():
            out.write(f">{hdr}\n")
            for i in range(0, len(seq), 70):
                out.write(seq[i:i+70] + "\n")

def parse_attrs(attr_str):
    d = OrderedDict()
    for item in attr_str.strip().split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            d[k.strip()] = v.strip()
    return d

def format_attrs(d):
    return ";".join(f"{k}={v}" for k, v in d.items())

def build_submission_note(attr_map):
    notes = []
    if attr_map.get("confidence_tier"): notes.append(f"confidence_tier: {attr_map['confidence_tier']}")
    if attr_map.get("description"): notes.append(f"description: {attr_map['description']}")
    if attr_map.get("eggnog_note"): notes.append(f"EggNOG: {attr_map['eggnog_note']}")
    if attr_map.get("interpro_note"): notes.append(f"InterPro: {attr_map['interpro_note']}")
    if attr_map.get("panther"): notes.append(f"PANTHER: {attr_map['panther']}")
    if attr_map.get("pfam"): notes.append(f"Pfam: {attr_map['pfam']}")
    if attr_map.get("canonical_transcript"): notes.append(f"canonical_transcript: {attr_map['canonical_transcript']}")
    if attr_map.get("Ontology_term_evidence"): notes.append(f"GO_evidence: {attr_map['Ontology_term_evidence']}")
    return " | ".join(notes)

def main():
    if len(sys.argv) < 6:
        sys.exit("Usage: purge_chrun_and_prep_submission.py <in_gff> <in_fa> <out_sub_gff> <out_chr_fa> <out_audit_tsv> <hap_label>")
    
    in_gff, in_fa, out_sub_gff, out_chr_fa, out_audit_tsv, hap_label = sys.argv[1:7]
    
    # 1. Filter FASTA to chromosomes only (discard ChrUn)
    raw_seqs = parse_fasta(in_fa)
    chr_seqs = OrderedDict((k, v) for k, v in raw_seqs.items() if "chrun" not in k.lower())
    write_fasta(chr_seqs, out_chr_fa)
    
    # 2. Filter GFF3, log purged ChrUn features, build submission records
    with open(in_gff) as f, open(out_sub_gff, 'w') as gff_out, open(out_audit_tsv, 'w') as audit_out:
        audit_out.write("Haplotype\tFeature_Type\tID\tParent\tChromosome\tStart\tEnd\tStrand\tProduct_or_Description\tFull_Note_Details\n")
        
        for line in f:
            if line.startswith('#'):
                # Strip sequence regions of ChrUn and PASA update comments
                if line.startswith('##sequence-region') and 'ChrUn' in line: continue
                if hap_label == 'Hap1' and 'PASA_UPDATE: g18420.t1' in line: continue
                if hap_label == 'Hap2' and 'PASA_UPDATE: g24783.t1' in line: continue
                gff_out.write(line)
                continue
            
            p = line.rstrip('\n').split('\t')
            if len(p) != 9: continue
            
            chrom, source, ftype, start, end, score, strand, phase, attr_str = p
            attrs = parse_attrs(attr_str)
            fid = attrs.get('ID', '.')
            fparent = attrs.get('Parent', '.')
            
            # Check if ChrUn feature
            if 'chrun' in chrom.lower() or '_UnG' in fid or '_UnG' in fparent:
                desc = attrs.get('product') or attrs.get('description') or '.'
                note_payload = build_submission_note(attrs) or '.'
                audit_out.write(f"{hap_label}\t{ftype}\t{fid}\t{fparent}\t{chrom}\t{start}\t{end}\t{strand}\t{desc}\t{note_payload}\n")
                continue
            
            # Submission attribute consolidation
            note_str = build_submission_note(attrs)
            new_attrs = OrderedDict()
            for k, v in attrs.items():
                if k in {'confidence_tier', 'eggnog_note', 'interpro_note', 'panther', 'pfam',
                        'canonical_transcript', 'canonical', 'description', 'Ontology_term_evidence', 'Note'}:
                    continue
                new_attrs[k] = v
            
            if note_str: new_attrs['Note'] = note_str
            
            p[8] = format_attrs(new_attrs)
            gff_out.write("\t".join(p) + "\n")

if __name__ == '__main__':
    main()
PYEOF4

# Execute Step 9 for both haplotypes
for HAP in Hap1 Hap2; do
    if [[ "${HAP}" == "Hap1" ]]; then 
        RAW_FA="${GENOME_FA_HAP1}"
    else 
        RAW_FA="${GENOME_FA_HAP2}"
    fi
    
    HAP_DIR="${FINAL_DIR}/HyPR01_${HAP}"
    INPUT_FULL_GFF="${HAP_DIR}/HyPR01_${HAP}_final_annotations_functional.gff3"
    SUBMISSION_GFF="${HAP_DIR}/HyPR01_${HAP}_final_annotations_functional_submission.gff3"
    CHR_ONLY_FA="${BASE}/hap$([[ "${HAP}" == "Hap1" ]] && echo "1" || echo "2").masked_chromosomes.fa"
    PURGED_AUDIT_TSV="${HAP_DIR}/HyPR01_${HAP}_ChrUn_purged_features_audit.tsv"
    
    python3 "${PURGE_AUDIT_PY}" \
        "${INPUT_FULL_GFF}" \
        "${RAW_FA}" \
        "${SUBMISSION_GFF}" \
        "${CHR_ONLY_FA}" \
        "${PURGED_AUDIT_TSV}" \
        "${HAP}"
        
    echo "[$(date)] ${HAP}: ChrUn purging complete."
    echo "  - Submission GFF3 : ${SUBMISSION_GFF}"
    echo "  - Chr-only FASTA  : ${CHR_ONLY_FA}"
    echo "  - Purged Audit Log: ${PURGED_AUDIT_TSV}"
done

# ──────────────────────────────────────────────────────────────────────────────
# Step 9b: Side-by-Side Comparative Feature Count (Full vs Purged vs Submission)
# ──────────────────────────────────────────────────────────────────────────────
COMPARE_SUMMARY="${FINAL_DIR}/HyPR01_submission_vs_full_feature_counts.txt"

python3 - << 'PY_COMPARE'
import sys
from collections import defaultdict

def parse_gff(path):
    genes = set()
    transcripts = set()
    cds_unique_tx = set()
    cds_segments = 0
    exons = 0
    genes_described = 0
    transcripts_described = 0
    utr5_by_tx = defaultdict(int)
    utr3_by_tx = defaultdict(int)
    
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            p = line.rstrip('\n').split('\t')
            if len(p) != 9: continue
            ftype = p[2]
            start, end = int(p[3]), int(p[4])
            attr_pairs = [kv.split('=', 1) for kv in p[8].split(';') if '=' in kv]
            attrs = {k.strip(): v.strip() for k, v in attr_pairs}
            
            fid = attrs.get('ID')
            fparent = attrs.get('Parent')
            
            if ftype == 'gene':
                genes.add(fid)
                desc = attrs.get('description', '').lower()
                if desc and desc != 'hypothetical protein': genes_described += 1
            elif ftype in ('mRNA', 'transcript'):
                transcripts.add(fid)
                prod = attrs.get('product', '').lower()
                if prod and prod != 'hypothetical protein': transcripts_described += 1
            elif ftype == 'CDS':
                cds_segments += 1
                if fparent: cds_unique_tx.add(fparent)
            elif ftype == 'exon':
                exons += 1
            elif ftype == 'five_prime_UTR':
                if fparent: utr5_by_tx[fparent] += (end - start + 1)
            elif ftype == 'three_prime_UTR':
                if fparent: utr3_by_tx[fparent] += (end - start + 1)
                
    return {
        'genes': len(genes),
        'tx': len(transcripts),
        'cds_tx': len(cds_unique_tx),
        'cds_seg': cds_segments,
        'exons': exons,
        'genes_desc': genes_described,
        'tx_desc': transcripts_described,
        'n5': len(utr5_by_tx),
        'n3': len(utr3_by_tx)
    }

base = "/path/to/your/directory"
haps = ['HyPR01_Hap1', 'HyPR01_Hap2']

report_lines = []
report_lines.append("="*85)
report_lines.append(f"{'FEATURE':<30} | {'FULL GENOME':<16} | {'PURGED (ChrUn)':<16} | {'SUBMISSION (Chr01-08)':<18}")
report_lines.append("="*85)

for h in haps:
    full_gff = f"{base}/{h}/{h}_final_annotations_functional.gff3"
    sub_gff = f"{base}/{h}/{h}_final_annotations_functional_submission.gff3"
    
    full_c = parse_gff(full_gff)
    sub_c = parse_gff(sub_gff)
    
    report_lines.append(f"\n--- {h} ---")
    keys = [
        ('Genes', 'genes'),
        ('Transcripts (mRNA)', 'tx'),
        ('Complete Coding Models (CDS)', 'cds_tx'),
        ('CDS Segments (Coding Exons)', 'cds_seg'),
        ('Exons', 'exons'),
        ('Genes with Functional Desc.', 'genes_desc'),
        ('Transcripts with Product', 'tx_desc'),
        ('mRNA with 5\' UTR', 'n5'),
        ('mRNA with 3\' UTR', 'n3')
    ]
    
    for label, k in keys:
        f_val = full_c[k]
        s_val = sub_c[k]
        diff = f_val - s_val
        report_lines.append(f"{label:<30} | {f_val:<16} | {diff:<16} | {s_val:<18}")

report_lines.append("\n" + "="*85)

summary_text = "\n".join(report_lines)
print(summary_text)

with open(f"{base}/HyPR01_submission_vs_full_feature_counts.txt", 'w') as out:
    out.write(summary_text + "\n")
PY_COMPARE

# ──────────────────────────────────────────────────────────────────────────────
# 8. Step 10: Run EMBLmyGFF3 Conversion for ENA Flatfile Generation
# ──────────────────────────────────────────────────────────────────────────────
if [[ -x "${EMBLMYGFF3_BIN}" ]]; then
    echo "[$(date)] Step 10: Executing EMBLmyGFF3 flatfile generation using chromosome FASTAs..."
    for HAP in Hap1 Hap2; do
        if [[ "${HAP}" == "Hap1" ]]; then 
            HAP_CODE="H1"
            GENOME_FA="${BASE}/hap1.masked_chromosomes.fa"
            ENTRY_DESC="Hypericum perforatum HyPR01 genome assembly, haplotype 1"
        else 
            HAP_CODE="H2"
            GENOME_FA="${BASE}/hap2.masked_chromosomes.fa"
            ENTRY_DESC="Hypericum perforatum HyPR01 genome assembly, haplotype 2"
        fi
        
        HAP_DIR="${FINAL_DIR}/HyPR01_${HAP}"
        SUBMISSION_GFF="${HAP_DIR}/HyPR01_${HAP}_final_annotations_functional_submission.gff3"
        DUMMY_TAG="Hp${ASSEMBLY}${HAP_CODE}"
        
        cat << 'MAP_EOF' > "${HAP_DIR}/mapping_attributes.json"
{
  "ID": "ID",
  "Parent": "Parent",
  "Name": "gene",
  "gene_name": "gene",
  "product": "product",
  "EC_number": "EC_number",
  "Note": "note",
  "Dbxref": "db_xref",
  "Ontology_term": "db_xref"
}
MAP_EOF

        (
            cd "${HAP_DIR}"
            "${EMBLMYGFF3_BIN}" "${SUBMISSION_GFF}" "${GENOME_FA}" \
                -o "HyPR01_${HAP}_${PROJECT_ID}.embl.gz" \
                -s "${ORGANISM}" \
                -p "${PROJECT_ID}" \
                -x "${TAX_DIV}" \
                -i "${DUMMY_TAG}" \
                -m "genomic DNA" \
                -t linear \
                -r 1 \
                -d WGS \
                --de "${ENTRY_DESC}" \
                --author "Rizzo P." "Chavez B." \
                --rg "${PROVIDER}" \
                --rt "${STUDY_TITLE}" \
                --email "${SUBMISSION_EMAIL}" \
                --use_attribute_value_as_locus_tag ID \
                --translate \
                --uncompressed_log \
                -z \
                2>&1 | tee "HyPR01_${HAP}_emblmygff3.log"
        )
        echo "[$(date)] ${HAP}: EMBL submission flatfile generated: ${HAP_DIR}/HyPR01_${HAP}_${PROJECT_ID}.embl.gz"
    done
else
    echo "[ERROR] EMBLmyGFF3 executable not found or not executable at ${EMBLMYGFF3_BIN}." >&2
    exit 1
fi

rm -f "${FILTER_PY}" "${STOPFILTER_PY}" "${ANNOT_PY}" "${CANON_PY}" "${PURGE_AUDIT_PY}"
echo "[*] End-to-end pseudomolecule pipeline completed successfully at $(date)"