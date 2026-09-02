#!/bin/bash
# ============================================================================
# HyPR01 - Gene-density vs TE-density co-plots  (LOCAL macOS version)
# ----------------------------------------------------------------------------
# Converts the original per-haplotype RStudio workflow into a single runnable
# script. Extracts coordinates from the gene (harmonized consensus) and TE
# (EDTA) GFF3 files, bins them into 500 kb windows, and renders paired
# gene/TE density tracks for both haplotypes as PNG (300 dpi), SVG, and PDF.
#
# Run:   bash HyPR01_Local_Gene_vs_TE_Density.sh
#   (or) chmod +x HyPR01_Local_Gene_vs_TE_Density.sh && ./HyPR01_Local_Gene_vs_TE_Density.sh
# ============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# 0. PATHS  -- edit here if your files move
# ---------------------------------------------------------------------------
BASE_DIR="/path/to/your/directory"
OUT_DIR="$BASE_DIR/Gene_density_vs_TE_density"
TSV_DIR="$OUT_DIR/intermediate_tsv"
SCRIPTS_DIR="$OUT_DIR/scripts"

# Gene annotations (harmonized consensus) + TE annotations (EDTA)
HAP1_GENE_GFF="$BASE_DIR/harmonized_consensus_Hap1.gff3"
HAP2_GENE_GFF="$BASE_DIR/harmonized_consensus_Hap2.gff3"
HAP1_TE_GFF="$BASE_DIR/Hap1_curated_chr_un_gapclosed.fasta.mod.EDTA.TEanno.gff3"
HAP2_TE_GFF="$BASE_DIR/Hap2_curated_chr_un_gapclosed.fasta.mod.EDTA.TEanno.gff3"
# NOTE: the *.masked_unchr.fa FASTAs are NOT needed for a density plot
#       (counts come from the GFF3 coordinates), so they are not used here.

mkdir -p "$OUT_DIR" "$TSV_DIR" "$SCRIPTS_DIR"

# ---------------------------------------------------------------------------
# 1. CHECK R IS INSTALLED
# ---------------------------------------------------------------------------
if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript not found on PATH." >&2
  echo "Install R for macOS from https://cran.r-project.org/bin/macosx/ then re-run." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# 2. CHECK INPUT FILES EXIST
# ---------------------------------------------------------------------------
for f in "$HAP1_GENE_GFF" "$HAP2_GENE_GFF" "$HAP1_TE_GFF" "$HAP2_TE_GFF"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing input file: $f" >&2
    exit 1
  fi
done

# ---------------------------------------------------------------------------
# 3. EXTRACT COORDINATES FROM GFF3  (chromosome, start, end)
# ---------------------------------------------------------------------------
echo "==> Extracting coordinates..."

extract_genes() {   # $1 = gff3 in, $2 = tsv out
  awk -F'\t' '$3=="gene"{print $1"\t"$4"\t"$5}' "$1" > "$2"
  if [ ! -s "$2" ]; then
    echo "  WARNING: no 'gene' features found in $1" >&2
    echo "  Feature types present (col 3), most common first:" >&2
    awk -F'\t' '!/^#/{print $3}' "$1" | sort | uniq -c | sort -rn | head >&2
  fi
}
extract_tes() {     # $1 = gff3 in, $2 = tsv out
  awk -F'\t' '!/^#/ && NF>=5 {print $1"\t"$4"\t"$5}' "$1" > "$2"
}

extract_genes "$HAP1_GENE_GFF" "$TSV_DIR/hap1_genes.tsv"
extract_genes "$HAP2_GENE_GFF" "$TSV_DIR/hap2_genes.tsv"
extract_tes   "$HAP1_TE_GFF"   "$TSV_DIR/hap1_repeats.tsv"
extract_tes   "$HAP2_TE_GFF"   "$TSV_DIR/hap2_repeats.tsv"

for f in hap1_genes hap1_repeats hap2_genes hap2_repeats; do
  n=$(wc -l < "$TSV_DIR/$f.tsv" | tr -d ' ')
  echo "    $f.tsv: $n rows"
done

# ---------------------------------------------------------------------------
# 4. WRITE THE R PLOTTER
# ---------------------------------------------------------------------------
cat << 'EOF' > "$SCRIPTS_DIR/render_density_local.R"
#!/usr/bin/env Rscript
# HyPR01 local density plotter -- one haplotype per invocation.

suppressWarnings(suppressMessages({
  ## ---- package setup (CRAN binaries on macOS; no compilation needed) ----
  lib <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(lib)) lib <- file.path(path.expand("~"), "Library", "R", "hypr01-lib")
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))

  need <- c("ggplot2", "dplyr", "tidyr", "ggnewscale", "scales")
  for (p in need) if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, lib = lib, repos = "https://cloud.r-project.org")
  # svglite is optional (nicer SVG); base svg() is used if unavailable
  if (!requireNamespace("svglite", quietly = TRUE))
    try(install.packages("svglite", lib = lib, repos = "https://cloud.r-project.org"),
        silent = TRUE)

  library(ggplot2); library(dplyr); library(tidyr); library(ggnewscale)
}))

args        <- commandArgs(trailingOnly = TRUE)
hap_label   <- args[1]
genes_tsv   <- args[2]
repeats_tsv <- args[3]
out_dir     <- args[4]
window_size <- 500000

read_coords <- function(path, what) {
  if (!file.exists(path) || file.info(path)$size == 0)
    stop("No ", what, " coordinates found in: ", path,
         "\n  (check the GFF3 feature types / paths)")
  df <- read.table(path, header = FALSE, sep = "\t",
                   col.names = c("Chromosome", "Start", "End"))
  df$Midpoint <- (df$Start + df$End) / 2
  df
}

genes   <- read_coords(genes_tsv,   "gene")
repeats <- read_coords(repeats_tsv, "TE")

## ---- normalise chromosome labels to lowercase chrNN / chrun ----
clean_labels <- function(x) {
  x <- tolower(x)
  x <- gsub("hpchr", "chr", x)
  x <- gsub("hypr01_hap1_chr", "chr", x)
  x <- gsub("hypr01_hap2_chr", "chr", x)
  x <- gsub("chrun_chr0", "chrun", x)
  x
}
genes$Clean   <- clean_labels(genes$Chromosome)
repeats$Clean <- clean_labels(repeats$Chromosome)

cat("\n[", hap_label, "] gene chromosomes: ",
    paste(sort(unique(genes$Clean)), collapse = ", "), "\n", sep = "")
cat("[", hap_label, "] TE   chromosomes: ",
    paste(sort(unique(repeats$Clean)), collapse = ", "), "\n", sep = "")
og  <- setdiff(unique(genes$Clean), unique(repeats$Clean))
orp <- setdiff(unique(repeats$Clean), unique(genes$Clean))
if (length(og))  cat("  WARNING: present in genes but not TEs: ", paste(og,  collapse = ", "), "\n")
if (length(orp)) cat("  WARNING: present in TEs but not genes: ", paste(orp, collapse = ", "), "\n")

## anything matching this is treated as "unplaced" and placed at the bottom
unplaced_rx <- "un|scaffold|contig|ctg|scf|tig|random|patch"

save_plot <- function(p, out_dir, suffix) {
  # PNG (300 dpi) via ggsave -- uses macOS quartz/ragg, high quality
  ggsave(file.path(out_dir, paste0(suffix, ".png")),
         plot = p, width = 14, height = 9, dpi = 300)
  # SVG (vector) -- svglite if available, else base svg() (cairo on macOS)
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file.path(out_dir, paste0(suffix, ".svg")), width = 14, height = 9)
  } else {
    svg(file.path(out_dir, paste0(suffix, ".svg")), width = 14, height = 9)
  }
  print(p); dev.off()
  # PDF (vector) -- base device, always reliable
  pdf(file.path(out_dir, paste0(suffix, ".pdf")), width = 14, height = 9)
  print(p); dev.off()
}

execute_plot_run <- function(gd, rd, suffix, with_unplaced) {
  if (nrow(gd) == 0 || nrow(rd) == 0) {
    cat("  (skipped ", suffix, ": no data)\n", sep = ""); return(invisible())
  }

  master_grid <- bind_rows(
    gd %>% select(Clean, Midpoint),
    rd %>% select(Clean, Midpoint)
  ) %>%
    mutate(BinStart = floor(Midpoint / window_size) * window_size) %>%
    group_by(Clean) %>%
    summarise(MaxBin = max(BinStart), .groups = "drop") %>%
    rowwise() %>%
    reframe(Clean = Clean, BinStart = seq(0, MaxBin, by = window_size))

  df_genes <- gd %>%
    mutate(BinStart = floor(Midpoint / window_size) * window_size) %>%
    group_by(Clean, BinStart) %>%
    summarise(Count = n(), .groups = "drop") %>%
    right_join(master_grid, by = c("Clean", "BinStart")) %>%
    mutate(Count = replace_na(Count, 0), Type = "Gene")

  df_repeats <- rd %>%
    mutate(BinStart = floor(Midpoint / window_size) * window_size) %>%
    group_by(Clean, BinStart) %>%
    summarise(Count = n(), .groups = "drop") %>%
    right_join(master_grid, by = c("Clean", "BinStart")) %>%
    mutate(Count = replace_na(Count, 0), Type = "Repeat")

  combined <- bind_rows(df_genes, df_repeats)

  ## Ordering: chr01 at the TOP descending to chrNN, unplaced at the very
  ## bottom. The y-axis is numeric and stacks BOTTOM-UP (first level -> y=1 ->
  ## bottom), so the level order from bottom to top must be:
  ##   unplaced, chrNN, ..., chr02, chr01
  ## i.e. unplaced first, then chromosomes highest-number-first. Sort by the
  ## numeric part so chr10+ never sorts before chr2.
  lv       <- unique(combined$Clean)
  std_all  <- lv[!grepl(unplaced_rx, lv)]
  unpl     <- sort(lv[grepl(unplaced_rx, lv)])
  key      <- suppressWarnings(as.numeric(gsub("[^0-9]", "", std_all)))
  std_asc  <- std_all[order(key, std_all)]          # chr01 -> chrNN
  final_levels <- c(unpl, rev(std_asc))
  combined$Chromosome_Display <- factor(combined$Clean, levels = final_levels)

  title_string <- paste("Gene vs. TE Density Across Gap-Filled HyPR01", toupper(hap_label))
  if (with_unplaced) title_string <- paste0(title_string, " (With Unplaced Scaffolds)")

  p <- ggplot() +
    geom_tile(data = subset(combined, Type == "Gene"),
              aes(x = BinStart, y = as.numeric(Chromosome_Display), fill = Count),
              height = 0.35) +
    scale_fill_viridis_c(option = "viridis", name = "Genes per\n500kb",
                         limits = c(0, 100), oob = scales::squish,
                         guide = guide_colorbar(order = 1)) +
    new_scale_fill() +
    geom_tile(data = subset(combined, Type == "Repeat"),
              aes(x = BinStart, y = as.numeric(Chromosome_Display) - 0.40, fill = Count),
              height = 0.35) +
    scale_fill_viridis_c(option = "inferno", name = "Transposable\nElements\nper 500kb",
                         limits = c(0, 800), oob = scales::squish,
                         guide = guide_colorbar(order = 2)) +
    scale_y_continuous(breaks = seq_along(final_levels) - 0.20, labels = final_levels,
                       expand = expansion(add = c(0.5, 0.5))) +
    scale_x_continuous(labels = function(x) paste0(x / 1e6, "Mb"),
                       expand = c(0, 0), position = "top") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(face = "bold", size = 15, color = "black"),
          axis.text.x = element_text(size = 14, color = "black"),
          axis.title = element_blank(),
          axis.ticks.x = element_line(color = "black"),
          axis.ticks.y = element_blank(),
          axis.line.x.top = element_line(color = "black", linewidth = 1.2),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 20,
                                    margin = margin(b = 15)),
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA)) +
    labs(title = title_string)

  save_plot(p, out_dir, suffix)
  cat("  saved: ", suffix, ".{png,svg,pdf}\n", sep = "")
}

base <- paste0("HyPR01_", hap_label, "_Gene_vs_TE_Density_500kb")

# Version A: all sequences incl. unplaced scaffolds
execute_plot_run(genes, repeats, paste0(base, "_with_unplaced"), with_unplaced = TRUE)

# Version B: anchored chromosomes only (unplaced removed)
gd_clean <- genes   %>% filter(!grepl(unplaced_rx, Clean))
rd_clean <- repeats %>% filter(!grepl(unplaced_rx, Clean))
execute_plot_run(gd_clean, rd_clean, paste0(base, "_chr_only"), with_unplaced = FALSE)

cat("[", hap_label, "] complete.\n", sep = "")
EOF

# ---------------------------------------------------------------------------
# 5. RENDER BOTH HAPLOTYPES
# ---------------------------------------------------------------------------
echo "==> Rendering Hap1..."
Rscript "$SCRIPTS_DIR/render_density_local.R" \
  "Hap1" "$TSV_DIR/hap1_genes.tsv" "$TSV_DIR/hap1_repeats.tsv" "$OUT_DIR"

echo "==> Rendering Hap2..."
Rscript "$SCRIPTS_DIR/render_density_local.R" \
  "Hap2" "$TSV_DIR/hap2_genes.tsv" "$TSV_DIR/hap2_repeats.tsv" "$OUT_DIR"

echo ""
echo "==================================================================="
echo " DONE. Figures written to:"
echo "   $OUT_DIR"
echo " Each haplotype produced two versions (with_unplaced + chr_only),"
echo " each as .png (300 dpi), .svg, and .pdf."
echo "==================================================================="
