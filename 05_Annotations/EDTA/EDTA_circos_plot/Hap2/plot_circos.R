#!/usr/bin/env Rscript
if (!requireNamespace("circlize", quietly = TRUE)) {
    install.packages("circlize", repos = "https://cloud.r-project.org")
}
library(circlize)

setwd("/path/to/your/directory")

# Read genome and filter unplaced contigs (e.g., chrUn) for the main plot
genome <- read.table("genome.txt", header = FALSE, stringsAsFactors = FALSE)
genome_df_all <- data.frame(chr = genome$V1, start = 0, end = genome$V2)
genome_df_main <- genome_df_all[!grepl("Un", genome_df_all$chr, ignore.case = TRUE), ]

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
            circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, col = "black", border = "black")
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
message("Rendering ", "Hap2", " plots...")

# Main Chromosomes Only
png(paste0("HyPR01_", "Hap2", "_TE_Landscape_MainChr.png"), width = 10, height = 10, units = "in", res = 600)
generate_circos_plot(genome_df_main)
dev.off() 

svg(paste0("HyPR01_", "Hap2", "_TE_Landscape_MainChr.svg"), width = 10, height = 10)
generate_circos_plot(genome_df_main)
dev.off()

# All Chromosomes (including unplaced contigs if present)
if (nrow(genome_df_all) > nrow(genome_df_main)) {
    png(paste0("HyPR01_", "Hap2", "_TE_Landscape_AllChr.png"), width = 10, height = 10, units = "in", res = 600)
    generate_circos_plot(genome_df_all)
    dev.off()
    
    svg(paste0("HyPR01_", "Hap2", "_TE_Landscape_AllChr.svg"), width = 10, height = 10)
    generate_circos_plot(genome_df_all)
    dev.off()
}
