#!/bin/bash
# ============================================================================
# HyPR01 - Comparative Gene Structure (Multi-species)  -- LOCAL macOS, R port
# ----------------------------------------------------------------------------
# Parses gene/CDS/exon features from 6 species' GFF3s (HyPR01 Hap2 + 5
# reference genomes), computes gene length / CDS length / exon length /
# exon number / intron length distributions, renders comparative figures,
# and outputs a detailed summary table (TSV) with species metrics breakdown.
# ============================================================================

set -uo pipefail

# ---------------------------------------------------------------------------
# 0. PATHS
# ---------------------------------------------------------------------------
FIG_BASE="/path/to/output/directory"
OUT_DIR="$FIG_BASE/Gene_Structure_plotting"
GFF3_DIR="$OUT_DIR/GFF3_files"
INTER_DIR="$OUT_DIR/intermediate"
LOG_DIR="$OUT_DIR/logs"
SCRIPTS_DIR="$OUT_DIR/scripts"

mkdir -p "$OUT_DIR" "$INTER_DIR" "$LOG_DIR" "$SCRIPTS_DIR"

TS=$(date +"%Y%m%d_%H%M%S")
LOGFILE="$LOG_DIR/run_${TS}.log"
exec > >(tee -a "$LOGFILE") 2>&1

echo "==================================================================="
echo " HyPR01 Comparative Gene Structure -- run started $(date)"
echo " Log file: $LOGFILE"
echo "==================================================================="

# ---------------------------------------------------------------------------
# 1. SPECIES CONFIGURATION
# ---------------------------------------------------------------------------
SP_FILES=(
  "$FIG_BASE/harmonized_consensus_Hap2.gff3"
  "$GFF3_DIR/Manihot_esculenta.M.esculenta_v8.63.gff3"
  "$GFF3_DIR/Erythroxylum_novogranatense.gff3"
  "$GFF3_DIR/Populus_trichocarpa.Pop_tri_v4.63.gff3"
  "$GFF3_DIR/Hevea_brasiliensis.gff3"
  "$GFF3_DIR/Arabidopsis_thaliana.TAIR10.63.gff3"
)
SP_NAMES=(
  "Hypericum perforatum"
  "Manihot esculenta"
  "Erythroxylum novogranatense"
  "Populus trichocarpa"
  "Hevea brasiliensis"
  "Arabidopsis thaliana"
)
SP_COLORS=(
  "#660099"
  "#F1C40F"
  "#00BFFF"
  "#FF7F50"
  "#2ECC71"
  "#708090"
)

# ---------------------------------------------------------------------------
# 2. CHECKS
# ---------------------------------------------------------------------------
if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript not found on PATH." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# 3. OPTIONAL AGAT DETECTION
# ---------------------------------------------------------------------------
AGAT_ENV="agat_env"
AGAT_CMD=""
if command -v mamba >/dev/null 2>&1; then
  for cmd in agat_convert_sp_gxf2gxf.pl agat_convert_sp_gxf2gxf; do
    if mamba run -n "$AGAT_ENV" which "$cmd" >/dev/null 2>&1; then
      AGAT_CMD="$cmd"
      break
    fi
  done
fi

# ---------------------------------------------------------------------------
# 4. WRITE THE R PARSER
# ---------------------------------------------------------------------------
cat << 'RSCRIPT_EOF' > "$SCRIPTS_DIR/parse_gff3.R"
#!/usr/bin/env Rscript
suppressWarnings(suppressMessages({
  lib <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(lib)) lib <- file.path(path.expand("~"), "Library", "R", "hypr01-lib")
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))
  if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table", lib = lib, repos = "https://cloud.r-project.org")
  library(data.table)
}))

args    <- commandArgs(trailingOnly = TRUE)
gff     <- args[1]
species <- args[2]
out_rds <- args[3]
out_qc  <- args[4]

get_attr <- function(x, key) {
  pat <- paste0(".*(?:^|;)", key, "=([^;]*).*")
  hit <- grepl(paste0("(?:^|;)", key, "="), x, ignore.case = TRUE, perl = TRUE)
  val <- rep(NA_character_, length(x))
  val[hit] <- sub(pat, "\\1", x[hit], ignore.case = TRUE, perl = TRUE)
  val
}

explode_parents <- function(d) {
  empty <- data.table(start = integer(), end = integer(), length = double(), parent = character())
  if (nrow(d) == 0) return(empty)
  d <- d[!is.na(parent) & parent != ""]
  if (nrow(d) == 0) return(empty)
  plist <- strsplit(d$parent, ",", fixed = TRUE)
  n_rep <- lengths(plist)
  data.table(
    start  = rep(d$start,  n_rep),
    end    = rep(d$end,    n_rep),
    length = rep(d$length, n_rep),
    parent = unlist(plist)
  )
}

message("Reading: ", gff)
dt <- fread(cmd = paste("grep -v '^#'", shQuote(gff)), sep = "\t", header = FALSE,
            quote = "", fill = TRUE, na.strings = c(),
            col.names = c("seqid","source","type","start","end","score","strand","phase","attributes"))
dt <- dt[!is.na(start) & !is.na(end)]
dt[, `:=`(start = as.numeric(start), end = as.numeric(end))]
dt[, length := end - start + 1]

n_total     <- nrow(dt)
dt[, feat_id     := get_attr(attributes, "ID")]
dt[, feat_parent := get_attr(attributes, "Parent")]

gene_rows <- dt[type == "gene"]
cds       <- dt[type == "CDS"]
exon      <- dt[type == "exon"]

tx_rows <- dt[type != "gene" & type != "CDS" & type != "exon" &
              !is.na(feat_id) & feat_id != "" &
              !is.na(feat_parent) & feat_parent != ""]
tx_map <- data.table(tx = character(), gene = character())
if (nrow(tx_rows) > 0) {
  plist  <- strsplit(tx_rows$feat_parent, ",", fixed = TRUE)
  n_rep  <- lengths(plist)
  tx_map <- data.table(tx = rep(tx_rows$feat_id, n_rep), gene = unlist(plist))
  tx_map <- unique(tx_map)
}

cds[,  parent := feat_parent]
exon[, parent := feat_parent]
cds_e  <- explode_parents(cds)
exon_e <- explode_parents(exon)

coding_tx <- unique(cds_e$parent)
coding_tx <- coding_tx[!is.na(coding_tx) & coding_tx != ""]

if (nrow(tx_map) > 0) {
  coding_genes <- unique(tx_map[tx %in% coding_tx, gene])
} else {
  coding_genes <- coding_tx
}
coding_genes <- coding_genes[!is.na(coding_genes) & coding_genes != ""]

n_genes_all <- nrow(gene_rows)

# Gene length = CDS span (earliest CDS start to latest CDS end per gene) --
# excludes UTRs, unlike the raw GFF3 "gene" feature's own start/end.
if (nrow(tx_map) > 0) {
  cds_gene <- merge(cds_e, tx_map, by.x = "parent", by.y = "tx")
} else {
  cds_gene <- copy(cds_e)
  cds_gene[, gene := parent]
}
cds_gene <- cds_gene[gene %in% coding_genes]
gene_span <- cds_gene[, .(g_start = min(start), g_end = max(end)), by = gene]
gene_span[, gene_len := g_end - g_start + 1]
genes <- gene_span$gene_len
n_genes_dropped <- n_genes_all - nrow(gene_span)

cds_e  <- cds_e[parent %in% coding_tx]
exon_e <- exon_e[parent %in% coding_tx]

cds_len_per_tx <- if (nrow(cds_e)) {
  agg <- cds_e[, .(cds_len = sum(length)), by = parent]
  agg[cds_len > 0, cds_len]
} else numeric(0)

seg <- if (nrow(exon_e) > 0) exon_e else cds_e
used_fallback <- nrow(exon_e) == 0 && nrow(cds_e) > 0

exon_len_values   <- if (nrow(seg)) seg$length else numeric(0)
exon_count_per_tx <- if (nrow(seg)) seg[, .N, by = parent]$N else integer(0)

introns <- numeric(0)
if (nrow(seg) > 0) {
  setorder(seg, parent, start)
  seg[, prev_end := shift(end), by = parent]
  ilen <- seg$start - seg$prev_end - 1
  introns <- ilen[!is.na(ilen) & ilen > 0]
}

metrics <- list(
  species = species, gene_length = genes, cds_length = cds_len_per_tx,
  exon_length = exon_len_values, exon_number = exon_count_per_tx, intron_length = introns,
  n_genes = length(genes), n_genes_all_features = n_genes_all, n_genes_dropped_noncoding = n_genes_dropped,
  n_transcripts_cds = length(cds_len_per_tx), n_transcripts_exon = length(exon_count_per_tx),
  used_cds_fallback_for_exon = used_fallback
)
saveRDS(metrics, out_rds)

status <- if (length(genes) > 0) "PASS" else "FAIL: no genes found"
writeLines(status, out_qc)

cat("\n==================================================\n")
cat(" FEATURE ACCEPTANCE METRICS FOR:", species, "\n")
cat("==================================================\n")
cat("  Total Raw GFF3 Rows Processed      :", n_total, "\n")
cat("  Accepted Protein-Coding Genes      :", length(genes), "\n")
cat("  Accepted Coding Transcripts (mRNAs):", length(coding_tx), "\n")
cat("  Accepted Functional Exon Segments  :", nrow(exon_e), "\n")
cat("  Accepted Functional CDS Segments   :", nrow(cds_e), "\n")
cat("  Accepted Functional Intron Gaps    :", length(introns), "\n")
cat("--------------------------------------------------\n")
cat("  QC Verdict Status                  :", status, "\n\n")
RSCRIPT_EOF

# ---------------------------------------------------------------------------
# 5. WRITE THE R PLOTTER & TSV EXPORTER
# ---------------------------------------------------------------------------
cat << 'RSCRIPT_EOF' > "$SCRIPTS_DIR/plot_gene_structure.R"
#!/usr/bin/env Rscript
suppressWarnings(suppressMessages({
  lib <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(lib)) lib <- file.path(path.expand("~"), "Library", "R", "hypr01-lib")
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))
  need <- c("ggplot2", "dplyr", "scales", "patchwork", "svglite")
  for (p in need) if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, lib = lib, repos = "https://cloud.r-project.org")
  library(ggplot2); library(dplyr); library(patchwork); library(svglite)
}))

args    <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
rest    <- args[-1]
specs <- data.frame(
  rds   = rest[seq(1, length(rest), 3)],
  name  = rest[seq(2, length(rest), 3)],
  color = rest[seq(3, length(rest), 3)],
  stringsAsFactors = FALSE
)
species_order <- specs$name
color_map     <- setNames(specs$color, specs$name)

# Helper function to compute safe statistical metrics without throwing errors on empty vectors
safe_stat <- function(vec, fn, digits = 2) {
  vec <- vec[!is.na(vec) & !is.nan(vec)]
  if (length(vec) == 0) return(NA)
  val <- fn(vec)
  if (is.numeric(val)) round(val, digits) else val
}

long_list    <- list()
summary_rows <- list()

for (i in seq_len(nrow(specs))) {
  m  <- readRDS(specs$rds[i])
  nm <- specs$name[i]
  
  # Extract metric vectors
  g_len <- m$gene_length
  c_len <- m$cds_length
  e_len <- m$exon_length
  e_num <- m$exon_number
  i_len <- m$intron_length

  # Construct a detailed summary row for this organism
  summary_rows[[i]] <- data.frame(
    species                     = nm,
    total_gff_raw_rows          = m$n_genes_all_features,
    accepted_coding_genes       = m$n_genes,
    dropped_noncoding_genes     = m$n_genes_dropped_noncoding,
    accepted_coding_transcripts = m$n_transcripts_cds,
    
    # Gene Length Statistics (bp)
    gene_length_mean            = safe_stat(g_len, mean),
    gene_length_median          = safe_stat(g_len, median),
    gene_length_sd              = safe_stat(g_len, sd),
    gene_length_min             = safe_stat(g_len, min, 0),
    gene_length_max             = safe_stat(g_len, max, 0),
    
    # CDS Length Statistics (bp)
    cds_length_mean             = safe_stat(c_len, mean),
    cds_length_median           = safe_stat(c_len, median),
    cds_length_sd               = safe_stat(c_len, sd),
    cds_length_min              = safe_stat(c_len, min, 0),
    cds_length_max              = safe_stat(c_len, max, 0),
    
    # Exon Length Statistics (bp)
    exon_length_mean            = safe_stat(e_len, mean),
    exon_length_median          = safe_stat(e_len, median),
    exon_length_sd              = safe_stat(e_len, sd),
    exon_length_min             = safe_stat(e_len, min, 0),
    exon_length_max             = safe_stat(e_len, max, 0),
    
    # Exons Per Transcript Statistics
    exons_per_tx_mean           = safe_stat(e_num, mean),
    exons_per_tx_median         = safe_stat(e_num, median),
    exons_per_tx_sd             = safe_stat(e_num, sd),
    exons_per_tx_min            = safe_stat(e_num, min, 0),
    exons_per_tx_max            = safe_stat(e_num, max, 0),
    
    # Intron Length Statistics (bp)
    intron_length_mean          = safe_stat(i_len, mean),
    intron_length_median        = safe_stat(i_len, median),
    intron_length_sd            = safe_stat(i_len, sd),
    intron_length_min           = safe_stat(i_len, min, 0),
    intron_length_max           = safe_stat(i_len, max, 0),
    
    stringsAsFactors = FALSE
  )

  # Prepare long format data frame for ggplot2 density plots
  for (metric in c("gene_length", "cds_length", "exon_length", "exon_number", "intron_length")) {
    v <- m[[metric]]
    if (length(v)) long_list[[length(long_list) + 1]] = data.frame(species = nm, metric = metric, value = as.numeric(v))
  }
}

# --- WRITE COMPREHENSIVE TSV SUMMARY TABLE ---
summary_df <- do.call(rbind, summary_rows)
tsv_out_path <- file.path(out_dir, "comparative_gene_structure_summary.tsv")
write.table(summary_df, file = tsv_out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Detailed comparative TSV written to:", tsv_out_path, "\n")

long_df         <- do.call(rbind, long_list)
long_df$species <- factor(long_df$species, levels = species_order)

panel_spec <- list(
  gene_length   = list(title = "Gene length",   xlab = "Gene length (bp)",   lim = 20000),
  cds_length    = list(title = "CDS length",    xlab = "CDS length (bp)",    lim = 6000),
  exon_length   = list(title = "Exon length",   xlab = "Exon length (bp)",   lim = 600),
  exon_number   = list(title = "Exon number",   xlab = "Exons per transcript", lim = 50),
  intron_length = list(title = "Intron length", xlab = "Intron length (bp)", lim = 1000)
)

base_theme <- theme_classic(base_size = 13) +
  theme(
    plot.title   = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title   = element_text(size = 12),
    legend.text  = element_text(face = "italic", size = 21), 
    legend.title = element_text(face = "bold", size = 18),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.3),
    plot.margin = margin(8, 14, 8, 8)
  )

make_panel <- function(metric_key, show_legend = FALSE) {
  spec <- panel_spec[[metric_key]]
  d <- long_df %>% dplyr::filter(metric == metric_key, value > 0)

  ylab_final <- if (metric_key == "exon_number") {
    "Proportion of transcripts"
  } else {
    expression(Density~(bp^-1))        
  }

  if (metric_key == "exon_number") {
    prop <- d %>%
      dplyr::mutate(value = round(value)) %>%
      dplyr::count(species, value, name = "k") %>%
      dplyr::group_by(species) %>%
      dplyr::mutate(proportion = k / sum(k)) %>%
      dplyr::ungroup()
      
    p <- ggplot(prop, aes(x = value, y = proportion, color = species)) +
      geom_line(linewidth = 0.9, show.legend = show_legend) +
      geom_point(size = 1.2, show.legend = show_legend) +
      scale_color_manual(values = color_map, name = "Species") +
      labs(title = spec$title, x = spec$xlab, y = ylab_final) +
      coord_cartesian(xlim = c(1, spec$lim)) +
      base_theme
      
    return(p)
  }

  # geom_density() shares one evaluation grid across all species; far
  # outliers in any one species stretch it past where the data actually
  # live, under-resolving the rest at the default n=512. Cap the input at
  # 5x the plotted range and raise resolution to n=8192 for stable peaks.
  density_cap <- 5 * spec$lim
  d_dens <- d %>% dplyr::filter(value <= density_cap)

  p <- ggplot(d_dens, aes(x = value, color = species)) +
    geom_density(fill = NA, linewidth = 0.9, key_glyph = "path", show.legend = show_legend, n = 8192) +
    scale_color_manual(values = color_map, name = "Species") +
    labs(title = spec$title, x = spec$xlab, y = ylab_final) +
    coord_cartesian(xlim = c(0, spec$lim)) +
    base_theme

  p
}

build_figure_layout <- function(suffix) {
  dummy_data <- data.frame(
    x = c(1, 2), y = c(1, 2),
    species = factor(rep(species_order, each = 2), levels = species_order)
  )
  
  p_inline_legend <- ggplot(dummy_data, aes(x = x, y = y, color = species)) +
    geom_line(alpha = 0) +
    scale_color_manual(values = color_map, name = NULL) +
    labs(title = "Species", x = NULL, y = NULL) +
    base_theme +
    theme(
      axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.48, 0.42),
      legend.key.height = unit(0.9, "cm"),
      legend.background = element_rect(fill = "transparent", color = NA)
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 2.2)))

  # Tags assigned manually so only the 5 data panels get lettered a-e; the
  # legend panel carries no tag.
  panels <- list(
    make_panel("gene_length",   show_legend = FALSE) + labs(tag = "a"),
    make_panel("cds_length",    show_legend = FALSE) + labs(tag = "b"),
    make_panel("exon_length",   show_legend = FALSE) + labs(tag = "c"),
    make_panel("exon_number",   show_legend = FALSE) + labs(tag = "d"),
    make_panel("intron_length", show_legend = FALSE) + labs(tag = "e"),
    p_inline_legend
  )

  fig <- wrap_plots(panels, ncol = 3) &
    theme(plot.tag = element_text(family = "Times New Roman", face = "bold", size = 34))
  
  png_path <- file.path(out_dir, paste0("HyPR01_GeneStructure_Comparative_", suffix, ".png"))
  svg_path <- file.path(out_dir, paste0("HyPR01_GeneStructure_Comparative_", suffix, ".svg"))
  
  ggsave(png_path, fig, width = 15, height = 9, dpi = 300)
  ggsave(svg_path, fig, width = 15, height = 9, device = svglite::svglite, fix_text_size = FALSE)
  
  cat("Saved graphics layout formats (PNG + SVG):", suffix, "\n")
}

# Linear plot only
build_figure_layout(suffix = "linear")
RSCRIPT_EOF

# ---------------------------------------------------------------------------
# 6. HELPERS: AGAT cleaning + per-species processing
# ---------------------------------------------------------------------------
clean_with_agat() {
  local raw="$1" cleaned="$2"
  [ -z "$AGAT_CMD" ] && return 1
  mamba run -n "$AGAT_ENV" "$AGAT_CMD" -g "$raw" -o "$cleaned" >/dev/null 2>&1
  [ -f "$cleaned" ]
}

PLOT_ARGS_FILE="$INTER_DIR/.plot_args_${TS}.txt"
: > "$PLOT_ARGS_FILE"

process_species() {
  local raw="$1" name="$2" color="$3"
  local safe="${name// /_}"
  local rds="$INTER_DIR/${safe}_metrics_v2.rds"
  local qc="$INTER_DIR/${safe}_v2.qc.txt"

  echo ""
  echo "--- $name ---"
  if [ ! -f "$raw" ]; then
    echo "WARNING: missing file for $name: $raw -- skipping this species." >&2
    return 1
  fi
  if [ -f "$rds" ]; then
    echo "-> using cached metrics: $rds"
    echo "${rds}|${name}|${color}" >> "$PLOT_ARGS_FILE"
    return 0
  fi

  echo "-> parsing $(basename "$raw")..."
  Rscript "$SCRIPTS_DIR/parse_gff3.R" "$raw" "$name" "$rds" "$qc"
  if [ ! -f "$qc" ]; then return 1; fi
  status=$(cat "$qc")

  if [[ "$status" == FAIL:* ]]; then
    if [ -n "$AGAT_CMD" ]; then
      cleaned="$INTER_DIR/${safe}_AGAT_cleaned.gff3"
      if clean_with_agat "$raw" "$cleaned"; then
        Rscript "$SCRIPTS_DIR/parse_gff3.R" "$cleaned" "$name" "$rds" "$qc"
      fi
    fi
  fi
  echo "${rds}|${name}|${color}" >> "$PLOT_ARGS_FILE"
}

# ---------------------------------------------------------------------------
# 7. RUN FOR ALL SPECIES
# ---------------------------------------------------------------------------
echo ""
echo "==================================================================="
echo " STEP 1: PARSE + QC"
echo "==================================================================="
for i in "${!SP_FILES[@]}"; do
  process_species "${SP_FILES[$i]}" "${SP_NAMES[$i]}" "${SP_COLORS[$i]}"
done

# ---------------------------------------------------------------------------
# 8. BUILD FIGURES & TSV TABLE
# ---------------------------------------------------------------------------
echo ""
echo "==================================================================="
echo " STEP 2: RENDERING COMPARATIVE FIGURES & EXPORTING SUMMARY TABLE"
echo "==================================================================="
PLOT_ARGS=()
while IFS='|' read -r rds name color; do
  [ -z "$rds" ] && continue
  PLOT_ARGS+=("$rds" "$name" "$color")
done < "$PLOT_ARGS_FILE"

Rscript "$SCRIPTS_DIR/plot_gene_structure.R" "$OUT_DIR" "${PLOT_ARGS[@]}"
echo "==================================================================="
echo " DONE. All outputs written to $OUT_DIR"
echo "==================================================================="