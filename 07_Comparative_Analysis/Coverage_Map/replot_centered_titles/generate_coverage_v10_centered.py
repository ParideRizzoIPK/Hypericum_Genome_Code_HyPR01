#!/usr/bin/env python3
"""
==============================================================================
HyPR01 Coverage Visualization Engine (v10 - Centered Titles Version)
==============================================================================
Processes competitive mapping BED and stats files into publication figures
and QC tables with centered figure titles and customized output paths.
"""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

# Global Styling
PLT_PARAMS = {
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 10.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "figure.titlesize": 12,
    "figure.titleweight": "bold",
    "axes.edgecolor": "#333333",
    "axes.linewidth": 0.8,
    "axes.grid": True,
    "grid.color": "#e0e0e0",
    "grid.linestyle": "-",
    "grid.linewidth": 0.5,
    "savefig.facecolor": "white",
}
plt.rcParams.update(PLT_PARAMS)

COLORS = {
    "Hap1": "#1f77b4",  # Steel Blue
    "Hap2": "#e65c00",  # Deep Orange
    "CP":   "#27ae60",  # Forest Green
    "MT":   "#c0392b",  # Crimson Red
    "ref":  "#444444",  # Slate Baseline
    "flag": "#d62728",  # Alert Red
    "gap":  "#7f7f7f",  # Assembly Gap
    "LSC":  "#2c3e50",  # Dark Slate Blue
    "IRa":  "#e74c3c",  # Crimson
    "IRb":  "#e74c3c",  # Crimson
    "SSC":  "#2980b9",  # Ocean Blue
}


def parse_chrom_meta(chrom: str) -> tuple[str, str, str]:
    """Parses full chromosome header into (base_chrom, haplotype, seq_class)."""
    if chrom.endswith("_Hap1"):
        return chrom[:-5], "Hap1", "nuclear"
    elif chrom.endswith("_Hap2"):
        return chrom[:-5], "Hap2", "nuclear"
    elif chrom.startswith("CP_"):
        return chrom, "Organelle", "chloroplast"
    elif chrom.startswith("MT_"):
        return chrom, "Organelle", "mitochondrion"
    else:
        return chrom, "Unknown", "nuclear"


def load_regions(bed_path: str, window: int):
    """Reads mosdepth regions BED and attaches sequence metadata."""
    if not os.path.exists(bed_path):
        sys.exit(f"[!] FATAL: mosdepth BED file not found: {bed_path}")

    df = pd.read_csv(
        bed_path,
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 3],
        names=["chrom", "start", "end", "depth"],
        dtype={0: str, 1: np.int64, 2: np.int64, 3: np.float64},
        compression="gzip",
    )

    lengths = df.groupby("chrom")["end"].max().astype(np.int64).to_dict()

    df["width"] = df["end"] - df["start"]
    df = df[df["width"] >= window].copy()

    df["mid_mb"] = (df["start"] + df["end"]) / 2e6
    df["mid_kb"] = (df["start"] + df["end"]) / 2e3

    parsed = [parse_chrom_meta(c) for c in df["chrom"]]
    df["base_chrom"] = [p[0] for p in parsed]
    df["haplotype"] = [p[1] for p in parsed]
    df["seq_class"] = [p[2] for p in parsed]

    df["in_gap"] = False
    return df, lengths


def parse_geseq_gff3(gff_path: str) -> list[tuple[float, float, str, str, float]]:
    """Parses GeSeq GFF3 to extract OGDRAW quadripartite structural features."""
    raw_regions = []
    if os.path.exists(gff_path):
        print(f"[*] Parsing GeSeq GFF3: {os.path.basename(gff_path)}")
        try:
            with open(gff_path, "r") as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 9 and parts[1] == "OGDRAW":
                        start_kb = int(parts[3]) / 1e3
                        end_kb = int(parts[4]) / 1e3
                        attr = parts[8]
                        size_kb = end_kb - start_kb
                        if "repreg1" in attr or "inverted repeat b" in attr.lower():
                            raw_regions.append((start_kb, end_kb, "IRb", COLORS["IRb"], size_kb))
                        elif "repreg2" in attr or "inverted repeat a" in attr.lower():
                            raw_regions.append((start_kb, end_kb, "IRa", COLORS["IRa"], size_kb))
                        elif "lsc" in attr or "large single copy" in attr.lower():
                            raw_regions.append((start_kb, end_kb, "LSC", COLORS["LSC"], size_kb))
                        elif "ssc" in attr or "small single copy" in attr.lower():
                            raw_regions.append((start_kb, end_kb, "SSC", COLORS["SSC"], size_kb))
        except Exception as e:
            print(f"[!] GFF3 parsing error ({e}). Falling back to defaults.")

    if not raw_regions:
        raw_regions = [
            (0.001, 16.557, "IRb", COLORS["IRb"], 16.556),
            (16.558, 112.225, "LSC", COLORS["LSC"], 95.667),
            (112.226, 128.783, "IRa", COLORS["IRa"], 16.557),
            (128.784, 139.646, "SSC", COLORS["SSC"], 10.862),
        ]

    return sorted(raw_regions, key=lambda r: r[0])


def summarise(df, lengths, haploid_ref_depth: float, edge_trim_bp: int) -> pd.DataFrame:
    """Per-sequence coverage metrics relative to the competitive haploid nuclear median."""
    rows = []
    for chrom, sub in df.groupby("chrom", sort=False):
        d = sub["depth"].to_numpy()
        length = int(lengths.get(chrom, sub["end"].max()))
        is_nuclear = sub["seq_class"].iloc[0] == "nuclear"

        keep = (sub["start"] >= edge_trim_bp) & (sub["end"] <= length - edge_trim_bp)
        if keep.sum() >= max(3, 0.4 * len(sub)):
            dt = sub.loc[keep, "depth"].to_numpy()
            mean_trim, med_trim, n_trim = dt.mean(), float(np.median(dt)), int(keep.sum())
        else:
            mean_trim, med_trim, n_trim = np.nan, np.nan, 0

        rows.append({
            "haplotype": sub["haplotype"].iloc[0],
            "chrom": chrom,
            "base_chrom": sub["base_chrom"].iloc[0],
            "seq_class": sub["seq_class"].iloc[0],
            "length_bp": length,
            "length_mb": length / 1e6,
            "n_windows": d.size,
            "mean_depth": d.mean(),
            "median_depth": float(np.median(d)),
            "n_windows_trimmed": n_trim,
            "mean_depth_trimmed": mean_trim,
            "median_depth_trimmed": med_trim,
            "sd_depth": d.std(ddof=1) if d.size > 1 else np.nan,
            "cv_depth": (d.std(ddof=1) / d.mean()) if d.size > 1 and d.mean() else np.nan,
            "p01_depth": np.percentile(d, 1),
            "p99_depth": np.percentile(d, 99),
            "max_depth": d.max(),
            "ratio_to_haploid_nuclear": d.mean() / haploid_ref_depth if haploid_ref_depth else np.nan,
            "pct_windows_in_0.5_to_2x_band":
                100 * np.mean((d >= 0.5 * haploid_ref_depth) & (d <= 2.0 * haploid_ref_depth)) if is_nuclear else np.nan,
            "pct_windows_lt_0.5x_ref":
                100 * np.mean(d < 0.5 * haploid_ref_depth) if is_nuclear else np.nan,
            "pct_windows_gt_2.0x_ref":
                100 * np.mean(d > 2.0 * haploid_ref_depth) if is_nuclear else np.nan,
            "n_zero_depth_windows": int(np.sum(d == 0)),
        })

    out = pd.DataFrame(rows)
    return out.sort_values(["seq_class", "haplotype", "base_chrom"]).reset_index(drop=True)


def parse_samtools_stats(path: str) -> dict:
    """Parses SN block from samtools stats output."""
    if not path or not os.path.exists(path):
        return {}
    keep = {
        "raw total sequences": "reads_total",
        "reads mapped": "reads_mapped",
        "reads unmapped": "reads_unmapped",
        "bases mapped (cigar)": "bases_mapped_cigar",
        "total length": "bases_total",
        "average length": "read_length_mean",
        "maximum length": "read_length_max",
        "error rate": "error_rate",
    }
    out = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("SN\t"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                key = parts[1].rstrip(":").strip()
                if key in keep:
                    try:
                        out[keep[key]] = float(parts[2].split("#")[0].strip())
                    except ValueError:
                        pass
    return out


def mapping_stats_row(stats, summary_df, haploid_ref_depth) -> dict:
    """Calculates overall mapping rate and implied nuclear genome size."""
    nuc = summary_df[summary_df["seq_class"] == "nuclear"]
    org = summary_df[summary_df["seq_class"] != "nuclear"]
    assembled_nuc = float(nuc["length_bp"].sum())
    assembled_org = float(org["length_bp"].sum())

    row = {
        "assembled_nuclear_total_bp": assembled_nuc,
        "assembled_haploid_equiv_bp": assembled_nuc / 2.0,
        "assembled_organellar_bp": assembled_org,
        "haploid_nuclear_median_depth": haploid_ref_depth,
    }
    row.update({k: stats.get(k, np.nan) for k in
                ["reads_total", "reads_mapped", "reads_unmapped", "bases_total",
                 "bases_mapped_cigar", "read_length_mean", "read_length_max", "error_rate"]})

    if stats.get("reads_total"):
        row["pct_reads_mapped"] = 100.0 * stats.get("reads_mapped", np.nan) / stats["reads_total"]
    else:
        row["pct_reads_mapped"] = np.nan

    bases = stats.get("bases_mapped_cigar", np.nan)
    if np.isfinite(bases) and haploid_ref_depth:
        org_bases = float((org["length_bp"] * org["mean_depth"]).sum())
        implied_haploid = (bases - org_bases) / haploid_ref_depth
        row["organellar_bases_mapped_est"] = org_bases
        row["implied_haploid_genome_bp"] = implied_haploid
        row["assembled_over_implied"] = (assembled_nuc / 2.0) / implied_haploid if implied_haploid > 0 else np.nan
    else:
        row["organellar_bases_mapped_est"] = np.nan
        row["implied_haploid_genome_bp"] = np.nan
        row["assembled_over_implied"] = np.nan
    return row


def detect_blocks(df, haploid_ref_depth: float, min_windows: int = 3) -> pd.DataFrame:
    """Merges consecutive anomalous windows into continuous blocks."""
    rows = []
    for chrom, sub in df[df["seq_class"] == "nuclear"].groupby("chrom", sort=False):
        sub = sub.sort_values("start").reset_index(drop=True)
        d = sub["depth"].to_numpy()
        for label, mask in (("high_gt_2x_haploid", d > 2.0 * haploid_ref_depth),
                            ("zero_depth", d == 0)):
            if not mask.any():
                continue
            idx = np.flatnonzero(mask)
            splits = np.flatnonzero(np.diff(idx) > 1) + 1
            for run in np.split(idx, splits):
                if run.size < min_windows:
                    continue
                seg = sub.iloc[run]
                rows.append({
                    "haplotype": seg["haplotype"].iloc[0],
                    "chrom": chrom,
                    "base_chrom": seg["base_chrom"].iloc[0],
                    "seq_class": sub["seq_class"].iloc[0],
                    "start": int(seg["start"].iloc[0]),
                    "end": int(seg["end"].iloc[-1]),
                    "span_kb": (int(seg["end"].iloc[-1]) - int(seg["start"].iloc[0])) / 1e3,
                    "n_windows": int(run.size),
                    "block_type": label,
                    "mean_depth": float(seg["depth"].mean()),
                    "max_depth": float(seg["depth"].max()),
                    "mean_over_haploid_nuclear": float(seg["depth"].mean()) / haploid_ref_depth if haploid_ref_depth else np.nan,
                })
    out = pd.DataFrame(rows)
    return out if not out.empty else pd.DataFrame(
        columns=["haplotype", "chrom", "base_chrom", "seq_class", "start", "end",
                 "span_kb", "n_windows", "block_type", "mean_depth", "max_depth", "mean_over_haploid_nuclear"]
    )


def plot_nuclear(df, lengths, chroms, hap, color, haploid_ref_depth, smooth, out_stem, formats):
    """Plots nuclear chromosome tracks with centered title and embedded legend."""
    print(f"[*] Plotting Nuclear Coverage with Legend Box for {hap}...")

    lens = {c: lengths[c] for c in chroms}
    max_len_mb = max(lens.values()) / 1e6

    n = len(chroms)
    panel_h_in, gap_in = 0.55, 0.12
    left_in, right_in, top_in, bottom_in = 0.85, 0.35, 0.80, 1.10
    plot_w_in = 8.8

    fig_w = left_in + plot_w_in + right_in
    fig_h = top_in + bottom_in + n * panel_h_in + (n - 1) * gap_in

    fig = plt.figure(figsize=(fig_w, fig_h))

    ymax = float(min(2.6 * haploid_ref_depth, np.percentile(df["depth"], 99.9) * 1.05))
    ymax = max(ymax, 1.8 * haploid_ref_depth)

    for i, chrom in enumerate(chroms):
        sub = df[df["chrom"] == chrom].sort_values("start")
        frac = (lens[chrom] / 1e6) / max_len_mb

        ax = fig.add_axes([
            left_in / fig_w,
            1 - (top_in + (i + 1) * panel_h_in + i * gap_in) / fig_h,
            (plot_w_in * frac) / fig_w,
            panel_h_in / fig_h,
        ])

        x = sub["mid_mb"].to_numpy()
        d = sub["depth"]

        med = d.rolling(smooth, center=True, min_periods=1).median().to_numpy()
        lo = d.rolling(smooth, center=True, min_periods=1).quantile(0.10).rolling(smooth, center=True, min_periods=1).mean().to_numpy()
        hi = d.rolling(smooth, center=True, min_periods=1).quantile(0.90).rolling(smooth, center=True, min_periods=1).mean().to_numpy()

        ax.axhline(haploid_ref_depth, color=COLORS["ref"], lw=0.7, zorder=1)
        for m in (0.5, 2.0):
            ax.axhline(haploid_ref_depth * m, color=COLORS["ref"], lw=0.6, ls=(0, (4, 3)), alpha=0.55, zorder=1)

        ax.fill_between(x, lo, hi, color=color, alpha=0.18, lw=0, rasterized=True, zorder=2)
        ax.plot(x, med, color=color, lw=0.85, zorder=3)

        over = sub.loc[sub["depth"] > ymax, "mid_mb"].to_numpy()
        if over.size:
            ax.plot(over, np.full(over.size, ymax * 0.985), ls="none", marker="v", ms=2.5, color=COLORS["flag"], clip_on=False, zorder=4)

        for _, w in sub[(sub["depth"] == 0)].iterrows():
            ax.axvline(w["mid_mb"], lw=0.6, alpha=0.85, zorder=4, color=COLORS["flag"])

        ax.set_xlim(0, lens[chrom] / 1e6)
        ax.set_ylim(0, ymax * 1.05)
        ax.set_yticks(sorted({0, int(round(haploid_ref_depth)), int(ymax // 10 * 10)}))
        
        base_name = sub["base_chrom"].iloc[0]
        ax.set_ylabel(base_name, fontsize=9.5, fontweight="bold", rotation=0, ha="right", va="center", labelpad=10)
        ax.grid(axis="y", alpha=0.45)
        ax.grid(axis="x", alpha=0.25)

        if i == 0:
            ax.xaxis.set_label_position("top")
            ax.xaxis.tick_top()
            ax.set_xlabel("Genomic position (Mb)", fontsize=10, fontweight="bold", labelpad=8)
        else:
            ax.set_xticklabels([])

    # CENTERED TITLE ADJUSTMENT
    fig.suptitle(
        f"HyPR01 {hap} \u2014 Competitive PacBio HiFi Read Depth Across Pseudomolecules",
        x=0.5, ha="center", y=0.985, fontweight="bold", fontsize=11,
    )

    ax_leg = fig.add_axes([left_in / fig_w, 0.08 / fig_h, plot_w_in / fig_w, 0.75 / fig_h])
    ax_leg.axis("off")

    legend_elements = [
        Line2D([0], [0], color=color, lw=1.2, label=f"Rolling Median ({smooth * 10} kb)"),
        Line2D([0], [0], color=color, alpha=0.35, lw=6, label="10th\u201390th Percentile Ribbon"),
        Line2D([0], [0], color=COLORS["ref"], lw=0.8, ls="-", label=f"Haploid Baseline ({haploid_ref_depth:.1f}\u00d7)"),
        Line2D([0], [0], color=COLORS["ref"], lw=0.6, ls="--", label="0.5\u00d7 / 2.0\u00d7 Expectations"),
        Line2D([0], [0], marker="v", color="w", markerfacecolor=COLORS["flag"], markersize=6, label=f"Off-scale Peak (>{ymax:.0f}\u00d7)"),
        Line2D([0], [0], color=COLORS["flag"], lw=1.0, label="Zero-coverage Window (0\u00d7)"),
    ]

    ax_leg.legend(handles=legend_elements, loc="center", ncol=3, frameon=True, facecolor="#f8f9fa", edgecolor="#cccccc", fontsize=8.0)

    for ext in formats:
        path = f"{out_stem}.{ext}"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        print(f"[+] Saved {path}")
    plt.close(fig)


def plot_unplaced(df, lengths, chroms, hap, color, haploid_ref_depth, out_stem, formats):
    """Plots unplaced scaffolds with a centered title on a dedicated log-scale axis."""
    if not chroms:
        return
    print(f"[*] Plotting Unplaced Sequence Coverage for {hap}...")

    n = len(chroms)
    fig, axes = plt.subplots(n, 1, figsize=(11, 2.9 * n), squeeze=False)
    axes = axes[:, 0]
    FLOOR = 0.5

    for ax, chrom in zip(axes, chroms):
        sub = df[df["chrom"] == chrom].sort_values("start")
        length_mb = lengths[chrom] / 1e6
        d = sub["depth"].to_numpy()

        ax.step(sub["mid_mb"], np.clip(d, FLOOR, None), where="mid", color=color, lw=0.7, zorder=3)
        ax.axhline(haploid_ref_depth, color=COLORS["ref"], lw=0.8, zorder=1)
        for m in (0.5, 2.0):
            ax.axhline(haploid_ref_depth * m, color=COLORS["ref"], lw=0.6, ls=(0, (4, 3)), alpha=0.55, zorder=1)

        ax.set_yscale("log")
        ax.set_ylim(FLOOR, max(d.max() * 1.6, haploid_ref_depth * 4))
        ax.set_xlim(0, length_mb)

        ax.set_title(
            f"{chrom} ({hap})   |   {length_mb:.2f} Mb   |   {len(sub):,} windows   |   "
            f"median {np.median(d):,.1f}\u00d7   |   max {d.max():,.0f}\u00d7",
            loc="left", fontsize=9, fontweight="bold",
        )
        ax.set_ylabel("Depth (\u00d7, log)", fontweight="bold", fontsize=9)
        ax.set_xlabel("Genomic position (Mb)", fontweight="bold", fontsize=9)
        ax.grid(axis="y", which="both", alpha=0.35)

    # CENTERED TITLE ADJUSTMENT
    fig.suptitle(f"HyPR01 {hap} \u2014 Unplaced Sequence Read Depth", x=0.5, ha="center", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    for ext in formats:
        path = f"{out_stem}.{ext}"
        fig.savefig(path, dpi=300)
        print(f"[+] Saved {path}")
    plt.close(fig)


def plot_hap_comparison(dists_df, chrom_order, haploid_ref_depth, out_stem, formats):
    """Renders side-by-side per-chromosome boxplots with centered figure title."""
    print("[*] Generating Competitive Haplotype Comparison Figure...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.4, 1]})

    sub_df = dists_df[dists_df["base_chrom"].isin(chrom_order)].copy()

    sns.boxplot(
        data=sub_df, x="base_chrom", y="depth", hue="haplotype",
        palette={"Hap1": COLORS["Hap1"], "Hap2": COLORS["Hap2"]},
        ax=ax1, showfliers=False, linewidth=0.8, width=0.6, order=chrom_order,
    )

    ax1.axhline(haploid_ref_depth, color=COLORS["ref"], lw=0.8, ls="-", zorder=1, label=f"Haploid Median ({haploid_ref_depth:.1f}\u00d7)")
    ax1.axhline(haploid_ref_depth * 0.5, color=COLORS["ref"], lw=0.6, ls="--", alpha=0.6, zorder=1)
    ax1.axhline(haploid_ref_depth * 2.0, color=COLORS["ref"], lw=0.6, ls="--", alpha=0.6, zorder=1)

    ax1.set_ylim(0, max(np.percentile(sub_df["depth"], 99.0) * 1.15, haploid_ref_depth * 2.5))
    ax1.set_ylabel("10 kb Window Depth (\u00d7)", fontweight="bold")
    ax1.set_xlabel("Chromosome", fontweight="bold")
    ax1.set_title("Per-chromosome Competitive Depth Distributions", loc="left", fontweight="bold")
    ax1.legend(frameon=True, facecolor="white", edgecolor="none", fontsize=8, loc="upper right")
    ax1.tick_params(axis="x", rotation=35)

    bins = np.linspace(0, max(np.percentile(sub_df["depth"], 99.0) * 1.15, haploid_ref_depth * 2.5), 100)
    for hap in ["Hap1", "Hap2"]:
        v = sub_df.loc[sub_df["haplotype"] == hap, "depth"].to_numpy()
        if v.size:
            ax2.hist(v, bins=bins, histtype="step", lw=1.3, density=True, color=COLORS[hap], label=hap)

    ax2.axvline(haploid_ref_depth, color=COLORS["ref"], lw=0.8, ls="-", zorder=1)
    ax2.set_xlabel("Window Depth (\u00d7)", fontweight="bold")
    ax2.set_ylabel("Density", fontweight="bold")
    ax2.set_title("Competitive Depth Density (10 kb)", loc="left", fontweight="bold")
    ax2.legend(frameon=True, facecolor="white", edgecolor="none", fontsize=8, loc="upper right")

    # CENTERED TITLE ADJUSTMENT
    fig.suptitle("HyPR01 Competitive Haplotype Read Partitioning (Hap1 vs Hap2)", x=0.5, ha="center", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    for ext in formats:
        path = f"{out_stem}.{ext}"
        fig.savefig(path, dpi=300)
        print(f"[+] Saved {path}")
    plt.close(fig)


def plot_organelles(df_org, lengths_org, haploid_ref_depth, geseq_gff3, out_stem, formats):
    """Renders organellar coverage at fine 200 bp resolution with GeSeq annotations."""
    chroms = sorted(df_org["chrom"].unique(), key=lambda c: (not c.startswith("CP_"), c))
    print(f"[*] Generating Fine-Resolution Organellar Figure ({len(chroms)} contigs)...")

    raw_plastid_regions = parse_geseq_gff3(geseq_gff3)
    
    n = len(chroms)
    fig, axes = plt.subplots(n, 1, figsize=(10.5, 2.8 * n), squeeze=False)
    axes = axes[:, 0]

    for ax, chrom in zip(axes, chroms):
        sub = df_org[df_org["chrom"] == chrom].sort_values("start")
        is_cp = chrom.startswith("CP_")
        base_color = COLORS["CP"] if is_cp else COLORS["MT"]
        length_kb = lengths_org[chrom] / 1e3

        x = sub["mid_kb"].to_numpy()
        y = sub["depth"].to_numpy()

        ax.fill_between(x, 0, y, step="mid", color=base_color, alpha=0.18, lw=0)
        ax.step(x, y, where="mid", color=base_color, lw=1.2)

        mean_d = sub["depth"].mean()
        ax.axhline(mean_d, color=COLORS["ref"], lw=0.8, ls=(0, (4, 3)))

        label = "Chloroplast" if is_cp else "Mitochondrion"
        clean = chrom.split("_", 1)[1]
        ymax_val = sub["depth"].max() * (1.35 if (is_cp and raw_plastid_regions) else 1.15)

        if is_cp and raw_plastid_regions:
            bar_y = sub["depth"].max() * 1.12
            bar_h = sub["depth"].max() * 0.14
            for s_k, e_k, reg_name, reg_col, orig_len in raw_plastid_regions:
                block_w = e_k - s_k
                if block_w > 0:
                    ax.add_patch(Rectangle((s_k, bar_y), block_w, bar_h, color=reg_col, alpha=0.92, ec="white", lw=0.8))
                    ax.text((s_k + e_k) / 2, bar_y + bar_h / 2,
                            f"{reg_name} ({orig_len:.1f} kb)" if block_w >= 12.0 else reg_name,
                            color="white", fontweight="bold", fontsize=8.0, ha="center", va="center")

        ax.set_ylim(0, ymax_val)
        ax.set_xlim(0, length_kb)
        ax.set_title(
            f"{label} \u2014 {clean}   |   {length_kb:.1f} kb   |   "
            f"mean depth {mean_d:,.0f}\u00d7   |   {mean_d / haploid_ref_depth:,.0f}\u00d7 fold haploid nuclear",
            loc="left", fontsize=9.2, fontweight="bold",
        )
        ax.set_ylabel("Depth (\u00d7)", fontweight="bold", fontsize=9)
        ax.set_xlabel("Position (kb)", fontweight="bold", fontsize=9)

        ax_fold = ax.twinx()
        curr_ylim = ax.get_ylim()
        ax_fold.set_ylim(curr_ylim[0] / haploid_ref_depth, curr_ylim[1] / haploid_ref_depth)
        ax_fold.set_ylabel("Fold Haploid Nuclear", fontweight="bold", fontsize=8.5, color="#555555")
        ax_fold.tick_params(axis="y", labelsize=8)
        ax_fold.grid(False)

    fig.suptitle("HyPR01 Organellar Read Depth (Competitive Mapping, 200 bp windows)", x=0.5, ha="center")
    fig.tight_layout(rect=(0, 0, 1, 0.965))

    for ext in formats:
        path = f"{out_stem}.{ext}"
        fig.savefig(path, dpi=300)
        print(f"[+] Saved {path}")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--bed-10kb", required=True)
    p.add_argument("--bed-200bp", required=True)
    p.add_argument("--geseq-gff3", required=True)
    p.add_argument("--samtools-stats", required=True)
    p.add_argument("--smooth", type=int, default=10)
    p.add_argument("--edge-trim-bp", type=int, default=15000)
    p.add_argument("--formats", nargs="+", default=["png", "svg"])
    args = p.parse_args()

    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    print(f"[*] Loading competitive 10 kb genome-wide profiles: {args.bed_10kb}")
    df_10k, lengths_10k = load_regions(args.bed_10kb, 10000)

    nuc = df_10k[df_10k["seq_class"] == "nuclear"]
    primary_nuc = nuc[nuc["base_chrom"].str.match(r"^Chr0[1-8]$")]
    haploid_ref_depth = float(np.median(primary_nuc["depth"]))
    print(f"[*] Haploid Nuclear Median Depth: {haploid_ref_depth:.2f}x")

    summary_df = summarise(df_10k, lengths_10k, haploid_ref_depth, args.edge_trim_bp)
    tsv_path = os.path.join(out_dir, "HyPR01_coverage_summary.tsv")
    summary_df.to_csv(tsv_path, sep="\t", index=False, float_format="%.4f")
    print(f"[+] Saved {tsv_path}")

    stats_dict = parse_samtools_stats(args.samtools_stats)
    mstats = mapping_stats_row(stats_dict, summary_df, haploid_ref_depth)
    mstats_path = os.path.join(out_dir, "HyPR01_mapping_stats.tsv")
    pd.DataFrame([mstats]).to_csv(mstats_path, sep="\t", index=False, float_format="%.4f")
    print(f"[+] Saved {mstats_path}")

    blocks_df = detect_blocks(df_10k, haploid_ref_depth)
    blocks_path = os.path.join(out_dir, "HyPR01_depth_blocks.tsv")
    blocks_df.to_csv(blocks_path, sep="\t", index=False, float_format="%.4f")
    print(f"[+] Saved {blocks_path}")

    chrom_order = [f"Chr0{i}" for i in range(1, 9)]
    piv = (summary_df[summary_df["base_chrom"].isin(chrom_order)]
           .pivot(index="base_chrom", columns="haplotype", values="length_bp")
           .reindex(chrom_order))
    if "Hap1" in piv.columns and "Hap2" in piv.columns:
        piv["delta_bp"] = piv["Hap2"] - piv["Hap1"]
        piv["delta_mb"] = piv["delta_bp"] / 1e6
        piv["pct_delta"] = 100.0 * piv["delta_bp"] / piv["Hap1"]
    length_tsv = os.path.join(out_dir, "HyPR01_haplotype_length_comparison.tsv")
    piv.to_csv(length_tsv, sep="\t", float_format="%.4f")
    print(f"[+] Saved {length_tsv}")

    for hap in ["Hap1", "Hap2"]:
        sub_hap = df_10k[df_10k["haplotype"] == hap]
        chroms_std = [f"Chr0{i}_{hap}" for i in range(1, 9) if f"Chr0{i}_{hap}" in sub_hap["chrom"].values]
        
        plot_nuclear(
            sub_hap[sub_hap["chrom"].isin(chroms_std)], lengths_10k, chroms_std, hap,
            COLORS[hap], haploid_ref_depth, args.smooth,
            os.path.join(out_dir, f"HyPR01_{hap}_nuclear_coverage"), args.formats,
        )

        unplaced = [c for c in sub_hap["chrom"].unique() if c not in chroms_std and sub_hap.loc[sub_hap["chrom"] == c, "seq_class"].iloc[0] == "nuclear"]
        if unplaced:
            plot_unplaced(
                sub_hap[sub_hap["chrom"].isin(unplaced)], lengths_10k, unplaced, hap,
                COLORS[hap], haploid_ref_depth, os.path.join(out_dir, f"HyPR01_{hap}_unplaced_coverage"), args.formats,
            )

    plot_hap_comparison(
        nuc[nuc["base_chrom"].isin(chrom_order)], chrom_order, haploid_ref_depth,
        os.path.join(out_dir, "HyPR01_haplotype_coverage_comparison"), args.formats,
    )

    if os.path.exists(args.bed_200bp):
        print(f"[*] Loading 200 bp organellar profiles: {args.bed_200bp}")
        df_200, lengths_200 = load_regions(args.bed_200bp, 200)
        df_org = df_200[df_200["seq_class"].isin(["chloroplast", "mitochondrion"])]
        if not df_org.empty:
            plot_organelles(
                df_org, lengths_200, haploid_ref_depth, args.geseq_gff3,
                os.path.join(out_dir, "HyPR01_organellar_coverage"), args.formats,
            )

    print(f"[+] COMPLETE: ALL RE-PLOTTED FIGURES SAVED TO {out_dir}")


if __name__ == "__main__":
    main()
