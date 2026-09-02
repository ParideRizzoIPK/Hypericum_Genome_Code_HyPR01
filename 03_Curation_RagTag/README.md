# 03 — Curation (RagTag)

Reference-guided curation of the YAHS scaffolds: pre-curation BUSCO QC
(eudicots_odb10), length filtering (≥100 kb), bidirectional RagTag
curation (minimap2) between haplotypes, and defining final chromosome IDs.

- `scripts/` — the curation pipeline itself, in run order (QC → length filter → RagTag → define chromosomes → rename/verify Hap2 → validation).
- `Validation_Dotplots/` — before/after synteny dot-plots and per-chromosome Merqury QV bar charts.

Full detail: [`../PIPELINE_OVERVIEW.md`](../PIPELINE_OVERVIEW.md#03_curation_ragtag--reference-guided-curation-and-validation).
