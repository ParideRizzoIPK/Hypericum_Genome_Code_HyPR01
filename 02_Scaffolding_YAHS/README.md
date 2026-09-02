# 02 — Scaffolding (YAHS)

Hi-C scaffolding of the hifiasm contigs into chromosome-scale scaffolds
with YAHS, plus the Juicebox/3D-DNA scripts used to manually review and
correct the Hi-C contact maps.

- `scripts/01_YAHS_scaffolding_HyPR01.sh` — YAHS + bwa-mem2 Hi-C scaffolding, telomere-aware (`CCCTAAA` motif).
- `scripts/` (remaining files) — the 3D-DNA/Juicer helper scripts for preparing and reconstructing the assembly around a manual Juicebox review.

Full detail: [`../PIPELINE_OVERVIEW.md`](../PIPELINE_OVERVIEW.md#02_scaffolding_yahs--hi-c-scaffolding-and-manual-curation).
