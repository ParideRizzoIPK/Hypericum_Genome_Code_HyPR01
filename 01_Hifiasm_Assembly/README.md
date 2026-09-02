# 01 — Hifiasm Assembly

The primary haplotype-resolved assembly step: hifiasm co-assembles the
PacBio HiFi reads with the Hi-C reads, producing the phased Hap1/Hap2
contig graphs everything downstream builds on.

- `01_HyPR01RAW_hifiasm_HiC_script.sbatch` — hifiasm assembly integrating Hi-C reads (`--h1`/`--h2`) for haplotype phasing.
- `Initial_Assembly_QUAST/` — QUAST QC on the raw haplotigs, before any scaffolding — the earliest checkpoint in the pipeline.

Full detail: [`../PIPELINE_OVERVIEW.md`](../PIPELINE_OVERVIEW.md#01_hifiasm_assembly--primary-haplotype-resolved-assembly).
