#!/usr/bin/env python3
"""
Build publication summary tables (CSV + Markdown) for the HyPR01
*Hypericum perforatum* organelle genomes:

  1. Chloroplast -- full quadripartite plastome summary (size, LSC/IR/SSC,
     unique gene counts, GC overall + per region).
  2. Mitochondrion -- the applicable rows (plant mitogenomes are NOT
     quadripartite, so LSC/IR/SSC do not apply); draft, 3 unresolved contigs,
     with a per-contig size/GC breakdown.

Everything is computed from the sequence + annotation files (nothing is
hand-typed), so the tables are reproducible and auditable.

Counting rules (documented, applied consistently):
  - "unique" = distinct genes, IR/segmental duplicate copies collapsed.
  - rRNA "unique" is counted by distinct /product species, which correctly
    (a) merges the chloroplast partial "-fragment" IR copies into their
    species and (b) merges the two mitochondrial 5S copies (rrn5, rrn5S).
  - protein-coding "unique" = distinct CDS /gene names, excluding partial
    "-fragment" annotations (e.g. the trans-spliced rps12 3' segment, so
    rps12 is counted once).
  - tRNA "unique" = distinct tRNA /gene names (anticodon-specific).

Chloroplast IR boundaries are taken from the GeSeq/OGDraw repeat_region
annotation and were independently confirmed by a blastn self-alignment
(16,558 bp, 100% identity, minus strand) during development.

Run inside `pycirclize_env` (Biopython available).
"""

import csv
import re
import warnings
from pathlib import Path
from Bio import SeqIO

warnings.filterwarnings("ignore")

ORG = Path("/path/to/your/directory")

SPECIES = "Hypericum perforatum (HyPR01)"
ACCESSION = "This study"

# --- chloroplast inputs / outputs ------------------------------------------
CP_DIR = ORG / "HyPR01_Chloroplast_data"
CP_FASTA = CP_DIR / "HyPR01_Chloroplast_Genome_rotated.fasta"
CP_GB = CP_DIR / ("HyPR01_Chloroplast_final_job-results-2026720154524/"
                  "20260720_HyPR01_Chloroplast_Geseq_annotations_"
                  "HyPR01_Chloroplast_Genome_rotated-circular-length%3D139646_GenBank.gb")
CP_CSV = CP_DIR / "HyPR01_Chloroplast_Genome_Summary.csv"
CP_MD = CP_DIR / "HyPR01_Chloroplast_Genome_Summary.md"

# --- mitochondrion inputs / outputs ----------------------------------------
MT_DIR = ORG / "HyPR01_Mitochondria"
MT_FLAT_GB = MT_DIR / "OGDraw_local_maps/fixed_gbk/combined_flattened.gb"
MT_CONTIG_GB = {
    "Contig 1": MT_DIR / "OGDraw_local_maps/fixed_gbk/contig1.gb",
    "Contig 2": MT_DIR / "OGDraw_local_maps/fixed_gbk/contig2.gb",
    "Contig 3": MT_DIR / "OGDraw_local_maps/fixed_gbk/contig3.gb",
}
MT_CSV = MT_DIR / "HyPR01_Mitochondria_Genome_Summary.csv"
MT_MD = MT_DIR / "HyPR01_Mitochondria_Genome_Summary.md"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def read_fasta(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def read_gb_origin(path):
    seq, in_origin = [], False
    for line in open(path):
        if line.startswith("ORIGIN"):
            in_origin = True
            continue
        if in_origin:
            if line.startswith("//"):
                break
            seq.append(re.sub(r"[^A-Za-z]", "", line))
    return "".join(seq).upper()


def gc(s):
    return (s.count("G") + s.count("C")) / len(s) * 100 if s else 0.0


def locus_length(path):
    for line in open(path):
        if line.startswith("LOCUS"):
            m = re.search(r"\s(\d+)\s+bp", line)
            if m:
                return int(m.group(1))
    raise ValueError(f"no LOCUS length in {path}")


def gene_counts(gb_path):
    """Return (protein_unique, trna_unique, rrna_species, details) with the
    documented counting rules applied."""
    rec = next(SeqIO.parse(str(gb_path), "genbank"))

    def qual(f, key):
        return f.qualifiers.get(key, [None])[0]

    is_frag = lambda n: bool(n) and re.search(r"frag|pseudo|partial", n, re.I)

    cds = {qual(f, "gene") for f in rec.features if f.type == "CDS"}
    cds = {g for g in cds if g and not is_frag(g)}

    trna = {qual(f, "gene") for f in rec.features if f.type == "tRNA"}
    trna = {g for g in trna if g and not is_frag(g)}

    # rRNA counted by distinct product species (merges fragments & copies)
    rrna_species = {qual(f, "product") for f in rec.features if f.type == "rRNA"}
    rrna_species = {p for p in rrna_species if p}

    details = {"cds": sorted(cds), "trna": sorted(trna),
               "rrna_species": sorted(rrna_species)}
    return len(cds), len(trna), len(rrna_species), details


def write_table(csv_path, md_path, title, header, rows, notes):
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            w.writerow(r)
    with open(md_path, "w") as fh:
        fh.write(f"# {title}\n\n")
        fh.write("| " + " | ".join(header) + " |\n")
        fh.write("|" + "|".join(["---"] * len(header)) + "|\n")
        for r in rows:
            fh.write("| " + " | ".join(str(x) for x in r) + " |\n")
        if notes:
            fh.write("\n**Notes**\n\n")
            for n in notes:
                fh.write(f"- {n}\n")
    print(f"wrote {csv_path.name} and {md_path.name}")


# ---------------------------------------------------------------------------
# chloroplast
# ---------------------------------------------------------------------------
def build_chloroplast():
    seq = read_fasta(CP_FASTA)
    L = len(seq)

    # IR boundaries from the GeSeq/OGDraw repeat_region annotation
    txt = open(CP_GB).read()
    irs = [(int(a), int(b)) for a, b in
           re.findall(r"repeat_region\s+(\d+)\.\.(\d+)", txt)]
    irs.sort()
    (irb_s, irb_e), (ira_s, ira_e) = irs[0], irs[1]
    ir_len = irb_e - irb_s + 1
    assert (ira_e - ira_s + 1) == ir_len, "IR copies differ in length"

    irb = seq[irb_s - 1:irb_e]
    ira = seq[ira_s - 1:ira_e]
    ssc = seq[irb_e:ira_s - 1]              # between IRB end and IRA start
    lsc = seq[ira_e:] + seq[:irb_s - 1]     # wraps the origin
    assert len(lsc) + len(irb) + len(ssc) + len(ira) == L

    prot, trna, rrna, det = gene_counts(CP_GB)
    total_genes = prot + trna + rrna

    header = ["Metric", SPECIES]
    rows = [
        ["Accession number", ACCESSION],
        ["Total chloroplast genome size (bp)", f"{L:,}"],
        ["Large single copy (LSC) (bp)", f"{len(lsc):,}"],
        ["Inverted repeat (IR) (bp)", f"{ir_len:,}"],
        ["Small single copy (SSC) (bp)", f"{len(ssc):,}"],
        ["Total number of genes (unique)", total_genes],
        ["Protein-coding gene (unique)", prot],
        ["tRNA (unique)", trna],
        ["rRNA (unique)", rrna],
        ["GC content (%)", f"{gc(seq):.2f}%"],
        ["LSC (%)", f"{gc(lsc):.2f}%"],
        ["IR (%)", f"{gc(ira):.2f}%"],
        ["SSC (%)", f"{gc(ssc):.2f}%"],
    ]
    notes = [
        "Assembled/annotated in this study; sequence rotated to the standard "
        "LSC-IRB-SSC-IRA orientation.",
        f"IR boundaries: IRB {irb_s:,}-{irb_e:,}, IRA {ira_s:,}-{ira_e:,} "
        f"(each {ir_len:,} bp), from the GeSeq/OGDraw annotation and "
        "independently confirmed by blastn self-alignment (100% identity).",
        "LSC spans the origin (IRA end -> genome end + start -> IRB start).",
        "Gene counts are unique genes from the GeSeq final annotation with IR "
        "duplicates collapsed; the trans-spliced rps12 is counted once "
        "(rps12-fragment excluded); rRNA counted as distinct species so partial "
        "IR-copy fragments (rrn16, rrn23) merge into their species.",
        f"rRNA species: {', '.join(det['rrna_species'])}.",
    ]
    write_table(CP_CSV, CP_MD, "HyPR01 chloroplast genome summary",
                header, rows, notes)


# ---------------------------------------------------------------------------
# mitochondrion
# ---------------------------------------------------------------------------
def build_mitochondrion():
    seq = read_gb_origin(MT_FLAT_GB)
    L = len(seq)
    contig_len = {n: locus_length(p) for n, p in MT_CONTIG_GB.items()}
    total = sum(contig_len.values())
    assert total == L, f"contig sum {total} != flattened length {L}"

    # per-contig slices (concatenated in order contig1+2+3)
    bounds, acc = {}, 0
    for n in ("Contig 1", "Contig 2", "Contig 3"):
        bounds[n] = (acc, acc + contig_len[n])
        acc += contig_len[n]

    prot, trna, rrna, det = gene_counts(MT_FLAT_GB)
    total_genes = prot + trna + rrna

    header = ["Metric", SPECIES]
    rows = [
        ["Accession number", ACCESSION],
        ["Assembly status", "Draft (3 unresolved contigs)"],
        ["Total mitochondrial genome size (bp)", f"{L:,}"],
        ["Number of contigs", len(contig_len)],
        ["Contig 1 size (bp)", f"{contig_len['Contig 1']:,}"],
        ["Contig 2 size (bp)", f"{contig_len['Contig 2']:,}"],
        ["Contig 3 size (bp)", f"{contig_len['Contig 3']:,}"],
        ["Total number of genes (unique)", total_genes],
        ["Protein-coding gene (unique)", prot],
        ["tRNA (unique)", trna],
        ["rRNA (unique)", rrna],
        ["GC content (%)", f"{gc(seq):.2f}%"],
        ["Contig 1 GC (%)", f"{gc(seq[slice(*bounds['Contig 1'])]):.2f}%"],
        ["Contig 2 GC (%)", f"{gc(seq[slice(*bounds['Contig 2'])]):.2f}%"],
        ["Contig 3 GC (%)", f"{gc(seq[slice(*bounds['Contig 3'])]):.2f}%"],
    ]
    notes = [
        "Assembled/annotated (PMGA) in this study.",
        "Plant mitochondrial genomes are NOT organized into the quadripartite "
        "LSC/IR/SSC structure of plastomes, so those rows are omitted; a "
        "per-contig size/GC breakdown is given instead.",
        "Draft assembly: the three contigs are not resolved into a single "
        "circular molecule and their mutual connectivity/order is unconfirmed "
        "(see HyPR01_Mitochondria_Assembly_Methods_and_Audit.md, sections 8, 13). "
        "'Total size' is the sum of the three contigs.",
        "Gene counts are unique genes; rRNA counted as distinct species -- two "
        "5S rRNA copies (rrn5, rrn5S) are present and collapse to one species "
        f"({', '.join(det['rrna_species'])}).",
    ]
    write_table(MT_CSV, MT_MD, "HyPR01 mitochondrial genome summary",
                header, rows, notes)


if __name__ == "__main__":
    build_chloroplast()
    build_mitochondrion()
