#!/bin/bash
# ============================================================================
# HyPR01 - Unified TE Landscape & Gene Collision Circos Plot Workflow
# ============================================================================

set -uo pipefail

# ---------------------------------------------------------
# 0. DIRECTORY & INPUT CONFIGURATION
# ---------------------------------------------------------
# Master Output Directory
OUT_BASE="/path/to/output/directory"

# Create output subdirectories and log directory
mkdir -p "$OUT_BASE/Hap1" "$OUT_BASE/Hap2" "$OUT_BASE/logs"

# Initialize detailed logging
TS=$(date +"%Y%m%d_%H%M%S")
LOGFILE="$OUT_BASE/logs/pipeline_run_${TS}.log"
exec > >(tee -a "$LOGFILE") 2>&1

echo "==================================================================="
echo " HyPR01 EDTA Circos Pipeline -- Run started $(date)"
echo " Detailed log saved to: $LOGFILE"
echo "==================================================================="

# ---------------------------------------------------------
# TOOL BINARY PATHS
# ---------------------------------------------------------
# Resolved via $PATH — activate the conda/mamba env that provides these
# (tested with bedtools in a `bedtools_env` env and samtools in a `samtools-env` env)
# before running, or override by exporting BEDTOOLS_BIN/SAMTOOLS_BIN yourself.
BEDTOOLS_BIN="${BEDTOOLS_BIN:-bedtools}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"

# Haplotype 1 Inputs
H1_FASTA="/path/to/input/directory/hap1.masked_unchr.fa"
H1_EDTA="/path/to/input/directory/Hap1_curated_chr_un_gapclosed.fasta.mod.EDTA.TEanno.gff3"
H1_GENES="/path/to/input/directory/harmonized_consensus_Hap1.gff3"

# Haplotype 2 Inputs
H2_FASTA="/path/to/input/directory/hap2.masked_unchr.fa"
H2_EDTA="/path/to/input/directory/Hap2_curated_chr_un_gapclosed.fasta.mod.EDTA.TEanno.gff3"
H2_GENES="/path/to/input/directory/harmonized_consensus_Hap2.gff3"

# ---------------------------------------------------------
# R SCRIPT GENERATOR (Dynamic for any Haplotype)
# ---------------------------------------------------------
create_r_script() {
  local work_dir="$1"
  local hap_name="$2"
  local r_script="$work_dir/plot_circos.R"

cat << EOF > "$r_script"
#!/usr/bin/env Rscript
if (!requireNamespace("circlize", quietly = TRUE)) {
    install.packages("circlize", repos = "https://cloud.r-project.org")
}
library(circlize)

setwd("$work_dir")

# Read genome and filter unplaced contigs (e.g., chrUn) for the main plot
genome <- read.table("genome.txt", header = FALSE, stringsAsFactors = FALSE)
genome_df_all <- data.frame(chr = genome\$V1, start = 0, end = genome\$V2)
genome_df_main <- genome_df_all[!grepl("Un", genome_df_all\$chr, ignore.case = TRUE), ]

# Read tracks
collisions <- read.table("track_collisions.txt", header = FALSE)
track_total <- read.table("track_total_TE.txt", header = FALSE)
track_gypsy <- read.table("track_gypsy.txt", header = FALSE)
track_copia <- read.table("track_copia.txt", header = FALSE)
track_line <- read.table("track_line.txt", header = FALSE)
track_dna <- read.table("track_dna.txt", header = FALSE)
track_helitron <- read.table("track_helitron.txt", header = FALSE)

generate_circos_plot <- function(target_genome) {
    circos.clear()
    circos.par("track.height" = 0.1, "start.degree" = 90, gap.degree = 2)
    circos.genomicInitialize(target_genome, plotType = c("axis", "labels"))

    circos.genomicTrack(collisions, ylim = c(0, 1), bg.border = NA, track.height = 0.05,
        panel.fun = function(region, value, ...) {
            circos.rect(CELL_META\$xlim[1], 0, CELL_META\$xlim[2], 1, col = "black", border = "black")
            circos.genomicRect(region, value, ytop = 1, ybottom = 0, col = "#74C476", border = NA)
        }
    )

    add_density_track <- function(data, color) {
        circos.genomicTrack(data, track.height = 0.08, bg.border = "grey90",
            panel.fun = function(region, value, ...) {
                circos.genomicLines(region, value, col = color, type = "h", lwd = 1)
            }
        )
    }

    add_density_track(track_total, "darkgrey")       
    add_density_track(track_gypsy, "#E64B35")        
    add_density_track(track_copia, "#4DBBD5")        
    add_density_track(track_line, "#E69F00")         
    add_density_track(track_dna, "#00A087")          
    add_density_track(track_helitron, "#8491B4")     

    legend_labels <- c("Nuclear Chromosomes", "Gene-TE collisions", "Total TE density", 
                       "Gypsy LTR-retro", "Copia LTR-retro", "LINEs", 
                       "DNA transposons", "Helitrons")
    legend_colors <- c("black", "#74C476", "darkgrey", "#E64B35", 
                       "#4DBBD5", "#E69F00", "#00A087", "#8491B4")

    legend("topright", legend = legend_labels, col = legend_colors, pch = 15,          
           cex = 0.7, pt.cex = 1.5, x.intersp = 0.8, y.intersp = 1.0,      
           bty = "n", inset = c(0.01, 0.01))
}

# Generate Outputs
message("Rendering ", "$hap_name", " plots...")

# Main Chromosomes Only
png(paste0("HyPR01_", "$hap_name", "_TE_Landscape_MainChr.png"), width = 10, height = 10, units = "in", res = 600)
generate_circos_plot(genome_df_main)
dev.off() 

svg(paste0("HyPR01_", "$hap_name", "_TE_Landscape_MainChr.svg"), width = 10, height = 10)
generate_circos_plot(genome_df_main)
dev.off()

# All Chromosomes (including unplaced contigs if present)
if (nrow(genome_df_all) > nrow(genome_df_main)) {
    png(paste0("HyPR01_", "$hap_name", "_TE_Landscape_AllChr.png"), width = 10, height = 10, units = "in", res = 600)
    generate_circos_plot(genome_df_all)
    dev.off()
    
    svg(paste0("HyPR01_", "$hap_name", "_TE_Landscape_AllChr.svg"), width = 10, height = 10)
    generate_circos_plot(genome_df_all)
    dev.off()
}
EOF
  chmod +x "$r_script"
}

# ---------------------------------------------------------
# MASTER PROCESSING FUNCTION
# ---------------------------------------------------------
process_haplotype() {
  local hap_name="$1"
  local work_dir="$2"
  local fasta="$3"
  local genes="$4"
  local edta="$5"

  echo "=================================================="
  echo " Processing $hap_name"
  echo "=================================================="
  cd "$work_dir"

  echo "1. Indexing genome and extracting lengths..."
  "$SAMTOOLS_BIN" faidx "$fasta"
  cut -f1,2 "${fasta}.fai" > genome.txt

  echo "2. Creating 100kb windows..."
  "$BEDTOOLS_BIN" makewindows -g genome.txt -w 100000 > windows_100kb.bed

  echo "3. Splitting TEs by superfamily..."
  grep -v "^#" "$edta" > all_TEs.gff3
  grep "Gypsy" all_TEs.gff3 > gypsy.gff3 || true
  grep "Copia" all_TEs.gff3 > copia.gff3 || true
  grep "LINE" all_TEs.gff3 > line.gff3 || true
  grep -i -E "TIR|mutator|hAT|CACTA|PIF|Harbinger" all_TEs.gff3 > dna_te.gff3 || true
  grep -i "Helitron" all_TEs.gff3 > helitron.gff3 || true

  echo "4. Calculating TE densities..."
  "$BEDTOOLS_BIN" coverage -a windows_100kb.bed -b all_TEs.gff3 | awk '{print $1"\t"$2"\t"$3"\t"$7}' > track_total_TE.txt
  "$BEDTOOLS_BIN" coverage -a windows_100kb.bed -b gypsy.gff3 | awk '{print $1"\t"$2"\t"$3"\t"$7}' > track_gypsy.txt
  "$BEDTOOLS_BIN" coverage -a windows_100kb.bed -b copia.gff3 | awk '{print $1"\t"$2"\t"$3"\t"$7}' > track_copia.txt
  "$BEDTOOLS_BIN" coverage -a windows_100kb.bed -b line.gff3 | awk '{print $1"\t"$2"\t"$3"\t"$7}' > track_line.txt
  "$BEDTOOLS_BIN" coverage -a windows_100kb.bed -b dna_te.gff3 | awk '{print $1"\t"$2"\t"$3"\t"$7}' > track_dna.txt
  "$BEDTOOLS_BIN" coverage -a windows_100kb.bed -b helitron.gff3 | awk '{print $1"\t"$2"\t"$3"\t"$7}' > track_helitron.txt

  echo "5. Finding TE-Gene collisions..."
  awk '$3=="gene" || $3=="mRNA"' "$genes" > genes_only.gff3
  "$BEDTOOLS_BIN" window -a genes_only.gff3 -b all_TEs.gff3 -w 1000 | awk '{print $1"\t"$4"\t"$5}' | sort -k1,1 -k2,2n | uniq > track_collisions.txt

  echo "6. Plotting..."
  create_r_script "$work_dir" "$hap_name"
  Rscript "$work_dir/plot_circos.R"
  
  echo "✔ $hap_name processing complete!"
  echo ""
}

# ---------------------------------------------------------
# EXECUTE
# ---------------------------------------------------------
process_haplotype "Hap1" "$OUT_BASE/Hap1" "$H1_FASTA" "$H1_GENES" "$H1_EDTA"
process_haplotype "Hap2" "$OUT_BASE/Hap2" "$H2_FASTA" "$H2_GENES" "$H2_EDTA"

echo "==================================================================="
echo "ALL PIPELINES SUCCESSFULLY COMPLETED."
echo "OUTPUTS LOCATED IN: $OUT_BASE"
echo "==================================================================="