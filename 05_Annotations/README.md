# 05 — Annotations

Repeat/transposable-element annotation, genome masking, gene prediction,
and functional annotation — the largest stage in the pipeline.

- `EDTA/` — transposable-element discovery (EDTA) and its Circos visualization.
- `Masking/` — soft-masking the genome from the EDTA output, then detecting and hard-masking organellar-derived (NUMT/NUPT) sequence.
- `Gene_density_vs_TE_density/` — combined gene/TE density visualization.
- `Gene_Annotations/` — the gene-model pipeline: RNA-seq evidence (STAR), ab initio + evidence-based prediction (ANNEVO, BRAKER3), GFF3 harmonization with PASA UTR annotation, functional annotation (eggNOG-mapper, InterProScan), and the final combined, submission-ready GFF3 (ChrUn purging, EMBLmyGFF3).

Full detail: [`../PIPELINE_OVERVIEW.md`](../PIPELINE_OVERVIEW.md#05_annotations--repeatte-masking-gene-prediction-functional-annotation).
