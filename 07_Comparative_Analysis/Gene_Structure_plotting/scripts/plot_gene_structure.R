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
