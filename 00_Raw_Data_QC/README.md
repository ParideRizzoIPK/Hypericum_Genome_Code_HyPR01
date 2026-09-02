# 00 — Raw Data QC

QC and pre-processing on the raw sequencing inputs, before any assembly
step: PacBio HiFi read QC, Hi-C read QC, k-mer profiling for genome-size
estimation, and RNA-seq trimming/QC.

- `Genome_Read_QC/` — NanoPlot QC on the raw HiFi reads.
- `HiC_Read_QC/` — FastQC on the raw Hi-C reads.
- `Kmer_Profiling_KMC/` — KMC k-mer counting (k=21) for genome-size/heterozygosity estimation.
- `RNAseq_Evidence_QC_and_Trimming/` — fastp/cutadapt trimming and QC of the RNA-seq evidence used later for annotation.

Full per-script detail: [`../PIPELINE_OVERVIEW.md`](../PIPELINE_OVERVIEW.md#00_raw_data_qc--input-data-qc-trimming-and-k-mer-profiling).
