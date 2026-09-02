# 07 — Comparative Analysis

Genome-wide QC, coverage validation, and comparison against other genomes
and species — the manuscript's Technical Validation section.

- `Genome_BUSCO/` — completeness checks against multiple lineage databases.
- `Coverage_Map/` — PacBio HiFi coverage mapping (minimap2 + mosdepth) and visualization.
- `Nuclear_Organizing_Region/` — rDNA/NOR confirmation (barrnap-HGV, BLASTn, nucmer, TRF).
- `HiC_Visualization/` — Hi-C contact-map heatmaps.
- `Comparison_vs_Hp_chinense/` — BUSCO and telomere comparison against the published tetraploid *H. p.* subsp. *chinense* genome.
- `QUAST_Comparison_vs_Hp_chinense/` — QUAST comparison against the same reference.
- `Synteny_Circos/` — MCScanX + Circos synteny plots, within HyPR01 and against Malpighiales reference species.
- `Gene_Structure_plotting/` — comparative gene-structure metrics across HyPR01 and 5 reference species.

Full detail: [`../PIPELINE_OVERVIEW.md`](../PIPELINE_OVERVIEW.md#07_comparative_analysis--genome-wide-qc-coverage-synteny-and-cross-species-comparison).
