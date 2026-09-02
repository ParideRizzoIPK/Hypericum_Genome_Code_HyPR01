#!/usr/bin/env python3
"""
Circular gene maps for the 3 HyPR01 mitochondrial contigs, from PMGA's
per-contig GenBank annotation files, using pyCirclize.

Environment: mamba env `pycirclize_env` (python=3.11, pycirclize, biopython).

Important note on circularity: none of the 3 contigs was assembled as a
confirmed closed circle by Flye/TIPPo/HiMT (see
HyPR01_Mitochondria_Assembly_Methods_and_Audit.md, sections 2-13) -- their
GenBank LOCUS lines already say "circular" as a PMGA template default, not
a verified structural claim. These maps are drawn in circular style for
visual consistency with standard organelle genome figures, per instruction,
but that is a cosmetic choice, not a claim of confirmed circularity --
state this explicitly in any figure legend/caption these are used in.

Outputs:
  - One combined figure: all 3 contigs as separate sectors in one Circos
    plot (mito_contigs_combined.png)
  - Three separate per-contig circular maps (mito_contig1/2/3_map.png)
"""

from pathlib import Path
from pycirclize import Circos
from pycirclize.parser import Genbank

BASE = Path("/path/to/your/directory")
GBK_DIR = BASE / "PMGA_Annotation_of_Mitochondria" / "01.Annotation"
OUT_DIR = BASE / "pyCirclize_maps"
OUT_DIR.mkdir(exist_ok=True)

CONTIGS = ["contig1", "contig2", "contig3"]
GBK_PATHS = {c: GBK_DIR / f"20250917150619818593.{c}.gb" for c in CONTIGS}

# Gene-category color scheme, matching the OGDraw legend convention already
# used for the HyPR01 chloroplast/mitochondrion figures in this project.
# Classified by gene-name prefix (case-insensitive), same groupings as
# PMGA's own report (see HyPR01_Mitochondria_Assembly_Methods_and_Audit.md
# section 4.1).
CATEGORY_COLORS = {
    "complex_i":        "#FFD700",  # nad*      - NADH dehydrogenase - yellow
    "complex_ii":       "#3CB371",  # sdh*      - succinate dehydrogenase - green
    "complex_iii":      "#9ACD32",  # cob       - cytochrome bc1 - yellow-green
    "complex_iv":       "#FFB6C1",  # cox*      - cytochrome c oxidase - pink
    "atp_synthase":     "#808000",  # atp*      - ATP synthase - olive
    "ribosomal_ssu":    "#DEB887",  # rps*      - ribosomal protein SSU - tan
    "ribosomal_lsu":    "#8B7355",  # rpl*      - ribosomal protein LSU - brown
    "maturase":         "#FF8C00",  # matR      - maturase - orange
    "other":            "#9370DB",  # ccm*, mttB, etc. - other genes - purple
    "trna":             "#191970",  # tRNA feature type - navy
    "rrna":             "#DC143C",  # rRNA feature type - red
    "unclassified_cds": "#B0B0B0",  # fallback - grey
}

PREFIX_CATEGORY = {
    "nad": "complex_i",
    "sdh": "complex_ii",
    "cob": "complex_iii",
    "cox": "complex_iv",
    "atp": "atp_synthase",
    "rps": "ribosomal_ssu",
    "rpl": "ribosomal_lsu",
    "matr": "maturase",
}


def classify_gene(gene_name: str, feature_type: str) -> str:
    if feature_type == "tRNA":
        return "trna"
    if feature_type == "rRNA":
        return "rrna"
    if not gene_name:
        return "unclassified_cds"
    name_lc = gene_name.lower()
    for prefix, category in PREFIX_CATEGORY.items():
        if name_lc.startswith(prefix):
            return category
    return "other"


def get_gene_name(feat) -> str:
    for key in ("gene", "product", "locus_tag"):
        if key in feat.qualifiers:
            return feat.qualifiers[key][0]
    return ""


def build_sector_tracks(sector, features):
    """Add forward/reverse feature tracks + tick track to one Circos sector."""
    f_track = sector.add_track((86, 95))
    r_track = sector.add_track((76, 85))
    tick_track = sector.add_track((95, 97))
    tick_track.axis(fc="none", ec="black", lw=1)
    tick_track.xticks_by_interval(
        interval=max(1, sector.size // 10),
        label_formatter=lambda v: f"{int(v/1000)} kb",
        label_size=6,
        show_endlabel=False,  # avoids the 0 kb / final tick colliding at the sector boundary
    )

    plotted_labels = []
    for feat in features:
        if feat.type not in ("CDS", "tRNA", "rRNA"):
            continue
        gene_name = get_gene_name(feat)
        category = classify_gene(gene_name, feat.type)
        color = CATEGORY_COLORS[category]
        strand = feat.location.strand
        track = f_track if strand == 1 else r_track
        track.genomic_features(feat, plotstyle="arrow", fc=color, lw=0.3, ec="black")
        if gene_name:
            mid = (int(feat.location.start) + int(feat.location.end)) / 2
            plotted_labels.append((mid, gene_name))

    # Gene name labels placed just outside the tick track via Sector.text(),
    # which natively supports r > 100 (no invalid-range track hack needed).
    for pos, name in plotted_labels:
        sector.text(name, x=pos, r=100, size=5, orientation="vertical")


def make_legend_handles():
    import matplotlib.patches as mpatches
    labels = {
        "complex_i": "Complex I (NADH dehydrogenase)",
        "complex_ii": "Complex II (succinate dehydrogenase)",
        "complex_iii": "Complex III (cytochrome bc1)",
        "complex_iv": "Complex IV (cytochrome c oxidase)",
        "atp_synthase": "ATP synthase",
        "ribosomal_ssu": "Ribosomal protein (SSU)",
        "ribosomal_lsu": "Ribosomal protein (LSU)",
        "maturase": "Maturase",
        "other": "Other genes",
        "trna": "tRNA",
        "rrna": "rRNA",
    }
    return [mpatches.Patch(color=CATEGORY_COLORS[k], label=v) for k, v in labels.items()]


def plot_single_contig(contig_name: str, gbk_path: Path, out_path: Path):
    gbk = Genbank(str(gbk_path))
    seqid2size = gbk.get_seqid2size()
    seqid2features = gbk.get_seqid2features(feature_type=None)

    circos = Circos(seqid2size, space=0)
    for sector in circos.sectors:
        build_sector_tracks(sector, seqid2features[sector.name])

    fig = circos.plotfig()
    fig.legend(handles=make_legend_handles(), loc="center", fontsize=6, frameon=False)
    fig.suptitle(
        f"HyPR01 mitochondrion — {contig_name} ({seqid2size[list(seqid2size)[0]]:,} bp)\n"
        f"circular-style rendering; assembly circularity not independently confirmed",
        fontsize=8,
    )
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"  wrote {out_path}")


def plot_combined(gbk_paths: dict, out_path: Path):
    seqid2size = {}
    seqid2features = {}
    for contig_name, path in gbk_paths.items():
        gbk = Genbank(str(path))
        s2s = gbk.get_seqid2size()
        s2f = gbk.get_seqid2features(feature_type=None)
        # Genbank parser uses the LOCUS name (e.g. "contig1") as seqid already
        for seqid in s2s:
            seqid2size[seqid] = s2s[seqid]
            seqid2features[seqid] = s2f[seqid]

    circos = Circos(seqid2size, space=3)
    for sector in circos.sectors:
        build_sector_tracks(sector, seqid2features[sector.name])
        sector.text(sector.name, size=8, r=115)

    fig = circos.plotfig()
    fig.legend(handles=make_legend_handles(), loc="center", fontsize=6, frameon=False)
    fig.suptitle(
        "HyPR01 mitochondrion — 3 assembled contigs (combined view)\n"
        "circular-style rendering; assembly circularity not independently confirmed; "
        "contigs shown separately, not joined",
        fontsize=8,
    )
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"  wrote {out_path}")


def main():
    print("Individual per-contig maps:")
    for contig in CONTIGS:
        out_path = OUT_DIR / f"mito_{contig}_map.png"
        plot_single_contig(contig, GBK_PATHS[contig], out_path)

    print("\nCombined map:")
    plot_combined(GBK_PATHS, OUT_DIR / "mito_contigs_combined.png")


if __name__ == "__main__":
    main()
