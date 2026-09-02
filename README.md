# HyPR01 *Hypericum perforatum* Genome — Analysis Scripts

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22258703.svg)](https://doi.org/10.5281/zenodo.22258703)

Scripts accompanying the chromosome-level, haplotype-resolved genome
assembly and annotation of *Hypericum perforatum* L. (genotype HyPR01),
submitted manuscript. This collection is deposited to make every
computational step reproducible: from raw-read QC through assembly,
scaffolding, curation, annotation, organellar genome reconstruction,
comparative analysis, and the figures/tables in the manuscript.

## Start here

| If you want to... | Go to... |
|---|---|
| Understand what each script does and the order stages run in | **`PIPELINE_OVERVIEW.md`** |
| Know what software (and which version) each stage needs | **`SOFTWARE_VERSIONS.md`** |

Each stage folder (`00`–`08`) also has its own short `README.md` with a
quick summary of what's inside.

## Directory structure

Folders `00`–`08` run in pipeline order:

```
00_Raw_Data_QC/
    PacBio HiFi read QC (NanoPlot) and Hi-C read QC (FastQC); KMC k-mer
    counting (k=21) for genome-size/heterozygosity estimation; RNA-seq
    evidence trimming and QC (fastp, cutadapt).

01_Hifiasm_Assembly/
    hifiasm haplotype-resolved assembly integrating Hi-C reads (--h1/--h2);
    QUAST QC on the raw, pre-scaffolding haplotigs.

02_Scaffolding_YAHS/
    YAHS + 3D-DNA/Juicer Hi-C scaffolding, telomere-aware (CCCTAAA motif);
    Hi-C contact maps manually reviewed and corrected in Juicebox.

03_Curation_RagTag/
    Pre-curation BUSCO QC (eudicots_odb10), length filtering (≥100 kb),
    bidirectional RagTag curation (minimap2) between haplotypes, final
    chromosome ID assignment, Hi-C and synteny validation of the result.

04_Gap_Closure_and_Final_QC/
    TGS-GapCloser gap closing using the raw HiFi reads; final QUAST +
    BUSCO + Merqury QC pass on the gap-closed assembly.

05_Annotations/
    EDTA transposable-element discovery and soft-masking; NUMT/NUPT
    detection and hard-masking; STAR RNA-seq alignment; BRAKER3 + ANNEVO
    gene prediction; GFF3 harmonization with PASA UTR annotation; eggNOG-
    mapper + InterProScan functional annotation; final combined,
    submission-ready GFF3 (ChrUn purging, EMBLmyGFF3).

06_Organellar_Genomes/
    TIPPo and HiMT independent chloroplast/mitochondrial assembly methods
    (cross-validated against each other); chloroplast genome rotation to
    the conventional start boundary; OGDraw/pyCirclize mitochondrial map
    rendering.

07_Comparative_Analysis/
    BUSCO completeness (multiple lineage databases), PacBio HiFi coverage
    mapping, NOR/rDNA confirmation, Hi-C contact-map visualization,
    comparison against *H. p.* chinense (synteny, BUSCO, QUAST, telomeres),
    MCScanX/Circos synteny plots, comparative gene-structure metrics across
    6 species.

08_Manuscript_Figures/
    Telomeric-repeat detection and per-chromosome ideograms.
```

See `PIPELINE_OVERVIEW.md` for the full per-script breakdown of every folder.

## Data availability

Sequencing data and assembled genome FASTAs referenced by these scripts are
deposited at the European Nucleotide Archive (ENA). **All FASTA files below
— nuclear, chloroplast, and mitochondrion — belong to the same genotype:
*Hypericum perforatum* HyPR01.**

| Accession | What it holds |
|---|---|
| [PRJEB123953](https://www.ebi.ac.uk/ena/browser/view/PRJEB123953) | *H. perforatum* HyPR01 **nuclear genome** — the PacBio HiFi, Hi-C, and RNA-seq reads generated for this study, and the haplotype-resolved nuclear genome assembly itself (Hap1 and Hap2 FASTA, ChrUn scaffolds excluded). Also the accession this pipeline's ENA submission scripts target (`05_Annotations/Gene_Annotations/Final_Combined_GFF3/`). Not yet publicly browsable at time of writing — pending release alongside the manuscript; check the link for current status. |
| [PRJEB124985](https://www.ebi.ac.uk/ena/browser/view/PRJEB124985) | *H. perforatum* HyPR01 **chloroplast genome** assembly. Not yet publicly browsable at time of writing — pending release alongside the manuscript; check the link for current status. |
| [PRJEB124986](https://www.ebi.ac.uk/ena/browser/view/PRJEB124986) | *H. perforatum* HyPR01 **mitochondrion genome** assembly. Not yet publicly browsable at time of writing — pending release alongside the manuscript; check the link for current status. |
| [PRJEB21346](https://www.ebi.ac.uk/ena/browser/view/PRJEB21346) | *H. perforatum* genotype HyPR01 — *"Identification of key regulators of dark gland development and hypericin biosynthesis in Hypericum perforatum"*: an earlier, already-public RNA-seq dataset (petal-rim dark-nodule vs. petal-center tissue, paired-end HiSeq2000) submitted by IPK Gatersleben and reused here as gene-prediction evidence (`00_Raw_Data_QC/RNAseq_Evidence_QC_and_Trimming/`, `05_Annotations/Gene_Annotations/Evidence_RNAseq_STAR/`). |

## Before running any `.sbatch` script

Every SLURM script in this collection writes its `--output`/`--error` logs
to a relative `logs/` subdirectory (so it works unedited under any
submission directory). **Create that directory first**, from inside
whichever stage folder you're submitting from:

```bash
mkdir -p logs
sbatch path/to/script.sbatch
```

SLURM will fail to submit the job (with no log file written to explain why)
if `logs/` doesn't already exist.

## Editing paths before running a script

Absolute filesystem paths have been replaced with placeholders for
publication (see the note at the end of `PIPELINE_OVERVIEW.md`):

- `/path/to/input/directory` — data a script reads, produced by an earlier stage.
- `/path/to/output/directory` — a script's own results.
- `/path/to/your/directory` — a shared project root, or a tool/environment install location.

Real filenames were kept, so the variable name plus filename together still
tell you exactly what a script expects as input or produces as output —
edit the placeholder prefix to match your own environment before running.

## Software

Five different provisioning methods are used across these scripts (HPC
environment modules, named conda/mamba/micromamba environments, Apptainer/
Singularity containers, Python venvs, and standalone downloaded binaries) —
see `SOFTWARE_VERSIONS.md` for the full tool-by-tool breakdown, including an
explicit list of what is *not* version-pinned in the scripts as written.

## Citation

This repository is archived on Zenodo with a permanent DOI:
[10.5281/zenodo.22258703](https://doi.org/10.5281/zenodo.22258703). Cite
this DOI (rather than a GitHub URL alone) if you use these scripts, so
your citation points to the exact archived version.

## License

This code is released under the Creative Commons Attribution 4.0
International License (CC-BY-4.0) — see `LICENSE`.
