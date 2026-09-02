#!/bin/bash
#SBATCH --job-name=Hp_Subgenomes_Telomeres
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/telomeres_%j.out
#SBATCH --error=logs/telomeres_%j.err

# ────────────────────────────────────────────────────────
# SCRIPT METADATA & VERSION CONTROL
# Pipeline Version: 4.0.0 (Tetraploid Subgenome Publication Build)
# Last Modified: 2026-07-16
# Description: Detects plant telomeric repeats (TTTAGGG/CCCTAAA) 
#              on concatenated subgenomes and outputs publication-ready plots.
# ────────────────────────────────────────────────────────

set -euo pipefail

# 1. Load System Modules & Environments
. /etc/profile.d/modules.sh
module load R/4.5.1

# 2. Establish Workspaces
INPUT_DIR="/path/to/input/directory"
OUT_DIR="/path/to/output/directory"
SCRIPTS_DIR="${OUT_DIR}/scripts"
LOGS_DIR="${OUT_DIR}/logs"

mkdir -p "${OUT_DIR}" "${SCRIPTS_DIR}" "${LOGS_DIR}"

# 3. Write Embedded Subgenome-Aware R Plotting Engine
cat << 'RSCRIPT_EOF' > "${SCRIPTS_DIR}/render_telomeres_unified.R"
#!/usr/bin/env Rscript
# Explicitly force cairo bitmap device type to run seamlessly on headless cluster nodes
options(bitmapType = "cairo")

suppressWarnings(suppressMessages({
  lib <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(lib)) lib = file.path(path.expand("~"), "Library", "R", "hypr01-lib")
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))
  library(ggplot2)
}))

args      <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript render_telomeres_unified.R <sub_A.fa> <sub_B.fa> <sub_C.fa> <sub_D.fa> <out_dir>")
}

fasta_sa  <- args[1]
fasta_sb  <- args[2]
fasta_sc  <- args[3]
fasta_sd  <- args[4]
out_dir   <- args[5]

WINDOW     <- 50000L   
MIN_COPIES <- 3L       
FWD <- "TTTAGGG"       
REV <- "CCCTAAA"       

clean_labels <- function(x) {
  x <- tolower(x)
  x <- gsub("hpchr|hypr01_hap1_chr|hypr01_hap2_chr|chrun_chr0", "Chr", x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

max_tandem <- function(s, unit) {
  if (nchar(s) == 0) return(0L)
  m <- gregexpr(sprintf("(?:%s){1,}", unit), s, perl = TRUE)[[1]]
  if (m[1] == -1) return(0L)
  as.integer(max(attr(m, "match.length")) %/% nchar(unit))
}
end_copies <- function(s) max(max_tandem(s, FWD), max_tandem(s, REV))

process_fasta <- function(fasta_path, subgenome_label) {
  message("Reading FASTA: ", fasta_path)
  lines   <- readLines(fasta_path, warn = FALSE)
  hdr_idx <- which(startsWith(lines, ">"))
  raw_names   <- sub("\\s.*$", "", sub("^>", "", lines[hdr_idx]))
  clean_names <- clean_labels(raw_names)
  rec_start   <- hdr_idx + 1L
  rec_end     <- c(hdr_idx[-1] - 1L, length(lines))

  records <- list()
  for (k in seq_along(hdr_idx)) {
    nm <- clean_names[k]
    if (grepl("^Chr[0-9]+[A-Da-d]?$", nm, ignore.case = TRUE)) {
      num <- suppressWarnings(as.integer(gsub("\\D", "", nm)))
      if (!is.na(num) && num >= 1 && num <= 8) {
        canon <- sprintf("Chr%02d", num)
        seqv  <- if (rec_end[k] >= rec_start[k]) toupper(paste0(lines[rec_start[k]:rec_end[k]], collapse = "")) else ""
        records[[canon]] <- seqv
      }
    }
  }

  canon_present <- sprintf("Chr%02d", 1:8)
  canon_present <- canon_present[canon_present %in% names(records)]

  res <- data.frame(Chromosome = character(), Length = numeric(),
                    Left_copies = integer(), Right_copies = integer(),
                    stringsAsFactors = FALSE)
  for (nm in canon_present) {
    s <- records[[nm]]
    L <- nchar(s)
    left_win  <- substr(s, 1, min(WINDOW, L))
    right_win <- substr(s, max(1, L - WINDOW + 1), L)
    res <- rbind(res, data.frame(Chromosome = nm, Length = L,
                                 Left_copies  = end_copies(left_win),
                                 Right_copies = end_copies(right_win),
                                 stringsAsFactors = FALSE))
  }
  res$Left_status  <- ifelse(res$Left_copies  >= MIN_COPIES, "Present", "Absent")
  res$Right_status <- ifelse(res$Right_copies >= MIN_COPIES, "Present", "Absent")
  res$Subgenome    <- subgenome_label
  return(res)
}

df_sa <- process_fasta(fasta_sa, "Subgenome A")
df_sb <- process_fasta(fasta_sb, "Subgenome B")
df_sc <- process_fasta(fasta_sc, "Subgenome C")
df_sd <- process_fasta(fasta_sd, "Subgenome D")
res   <- rbind(df_sa, df_sb, df_sc, df_sd)

# ---- Layout Calibration Math ---------------------------------
res <- res[order(res$Subgenome, as.integer(gsub("\\D", "", res$Chromosome))), ]
res$y <- as.numeric(as.factor(res$Chromosome))
res$y <- max(res$y) - res$y + 1 
ry <- 0.26
res$ymin <- res$y - ry
res$ymax <- res$y + ry

maxLen <- max(res$Length)
xlo <- -0.06 * maxLen; xhi <- 1.08 * maxLen; xrange <- xhi - xlo

ylo <- 0.3; yhi <- 9.2; yrange <- yhi - ylo
A   <- 0.45 
rx  <- xrange * ry * A / yrange

arc_pts <- function(cx, cy, side, grp, sub_lbl) {
  th <- seq(-pi/2, pi/2, length.out = 64)
  x  <- if (side == "right") cx + rx * cos(th) else cx - rx * cos(th)
  data.frame(x = x, y = cy + ry * sin(th), grp = grp, Subgenome = sub_lbl)
}

cap_list <- list()
lab  <- data.frame(x = numeric(), y = numeric(), label = character(),
                   hj = numeric(), Subgenome = character(), stringsAsFactors = FALSE)
gi <- 0

for (i in seq_len(nrow(res))) {
  cy <- res$y[i]; L <- res$Length[i]; sub_lbl <- res$Subgenome[i]
  
  if (res$Left_status[i] == "Present") {
    gi <- gi + 1; cap_list[[gi]] <- arc_pts(0, cy, "left", gi, sub_lbl)
    lab <- rbind(lab, data.frame(x = rx * 0.3, y = cy + ry + 0.16, label = paste0("n = ", res$Left_copies[i]), hj = 0, Subgenome = sub_lbl))
  }
  if (res$Right_status[i] == "Present") {
    gi <- gi + 1; cap_list[[gi]] <- arc_pts(L, cy, "right", gi, sub_lbl)
    lab <- rbind(lab, data.frame(x = L - rx * 0.3, y = cy + ry + 0.16, label = paste0("n = ", res$Right_copies[i]), hj = 1, Subgenome = sub_lbl))
  }
}
cap_df <- if (length(cap_list)) do.call(rbind, cap_list) else data.frame()

blue <- "#1F78B4"; blue_dk <- "#0D3B66"; body <- "#DDDDDD"; edge <- "#555555"
FONT_FAMILY <- "serif" 

p <- ggplot() +
  geom_rect(data = res, aes(xmin = 0, xmax = Length, ymin = ymin, ymax = ymax),
            fill = body, color = edge, linewidth = 0.4)

if (nrow(cap_df) > 0) {
  p <- p + geom_polygon(data = cap_df, aes(x = x, y = y, group = grp),
                        fill = blue, color = blue_dk, linewidth = 0.4)
}

p <- p +
  geom_text(data = lab, aes(x = x, y = y, label = label, hjust = hj),
            size = 3.8, fontface = "bold", color = "grey20", family = FONT_FAMILY) +
  facet_grid(Subgenome ~ .) +
  scale_y_continuous(breaks = 1:8, labels = sprintf("Chr%02d", 8:1), expand = c(0, 0), limits = c(ylo, yhi)) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6, " Mb"),
                     limits = c(xlo, xhi), expand = c(0, 0), position = "top") +
  theme_minimal(base_size = 16, base_family = FONT_FAMILY) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3), 
        strip.text   = element_text(face = "bold", size = 18),
        strip.background = element_rect(fill = "grey93", color = NA),
        axis.text.y  = element_text(face = "bold", size = 15, color = "black"),
        axis.text.x  = element_text(size = 13, color = "black"),
        axis.title   = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.spacing = unit(3, "lines"),
        plot.margin = margin(25, 25, 20, 20),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# ---- Save High-Res Formats (Ordered sequentially per user request) ----
message("Saving layout to PDF...")
ggsave(file.path(out_dir, "chinense_Subgenomes_Comparative_Telomeres.pdf"), p, width = 14, height = 18)

message("Saving layout to JPG...")
ggsave(file.path(out_dir, "chinense_Subgenomes_Comparative_Telomeres.jpg"), p, width = 14, height = 18, dpi = 300)

message("Saving layout to PNG...")
ggsave(file.path(out_dir, "chinense_Subgenomes_Comparative_Telomeres.png"), p, width = 14, height = 18, dpi = 300)

message("Saving layout to SVG...")
ggsave(file.path(out_dir, "chinense_Subgenomes_Comparative_Telomeres.svg"), p, width = 14, height = 18, device = svglite::svglite)

write.csv(res, file.path(out_dir, "chinense_Subgenomes_Comparative_Telomere_Summary.csv"), row.names = FALSE)
cat("Graphics saved to:", out_dir, "\n")
RSCRIPT_EOF

chmod +x "${SCRIPTS_DIR}/render_telomeres_unified.R"

# 4. Execute Analysis Pipeline
echo "[$(date)] ==> Initiating Headless Multi-Subgenome Telomere Plotter Pass..."

Rscript "${SCRIPTS_DIR}/render_telomeres_unified.R" \
  "${INPUT_DIR}/subgenome_A.fa" \
  "${INPUT_DIR}/subgenome_B.fa" \
  "${INPUT_DIR}/subgenome_C.fa" \
  "${INPUT_DIR}/subgenome_D.fa" \
  "${OUT_DIR}"

echo "[$(date)] ==> Done! Review output figures inside: ${OUT_DIR}"