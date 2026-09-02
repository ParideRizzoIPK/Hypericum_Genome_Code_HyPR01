# 04 — Gap Closure and Final QC

Closes remaining assembly gaps on the curated, anchored chromosomes using
the raw HiFi reads, then runs the final whole-assembly QC pass.

- `Gap_Closing/` — TGS-GapCloser run.
- `Final_QC/` — final QUAST + BUSCO + Merqury QC on the gap-closed assembly, including the per-chromosome Merqury QV bar chart.

Full detail: [`../PIPELINE_OVERVIEW.md`](../PIPELINE_OVERVIEW.md#04_gap_closure_and_final_qc--long-read-gap-closing-and-whole-assembly-qc).
