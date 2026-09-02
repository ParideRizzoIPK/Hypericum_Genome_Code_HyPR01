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
