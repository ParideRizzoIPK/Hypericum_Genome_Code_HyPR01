# HyPR01 Pipeline — Script Overview

This collection accompanies the *Hypericum perforatum* (HyPR01) genome
assembly and annotation. Folders `00`–`08` are grouped **by function**
(mirroring how a genome paper's Methods section reads) — each stage
generally consumes the previous stage's output. Paths inside the scripts
have been generalized for
publication (see the note at the end of this file); tool names, parameters,
and logic are unchanged from what was actually run.

Each stage folder also has its own short `README.md`. For software
names/versions, see **`SOFTWARE_VERSIONS.md`**.

---

## 00_Raw_Data_QC — input data QC, trimming, and k-mer profiling

Raw PacBio HiFi and Illumina Hi-C and Illumina RNA-seq inputs, their initial QC, and genome-size
estimation, before any assembly step.

- `Genome_Read_QC/HyPR01_nanoplot.sbatch` — NanoPlot QC report (read-length/quality plots) on the raw PacBio HiFi CCS reads.
- `HiC_Read_QC/HyPR01_HiC_fastqc.sbatch` — FastQC quality report on the raw Illumina Hi-C paired-end reads (R1/R2).
- `Kmer_Profiling_KMC/HyPR01_RAW_KMC.sbatch` — k-mer counting (KMC, k=21) on the raw HiFi reads and histogram export, for genome-size/heterozygosity estimation (e.g. via GenomeScope).
- `RNAseq_Evidence_QC_and_Trimming/HyPR01_Root_Transcriptome_fastp_QC.sbatch` — fastp adapter/quality trimming of the 2022 root RNA-seq reads (later used as BRAKER3/PASA evidence).
- `RNAseq_Evidence_QC_and_Trimming/trim_Hperforatum_H01.sh` — legacy in-house trimming pipeline (cutadapt + CLC Assembly Cell) for the 2022 leaf RNA-seq lanes; unzips, pairs, quality-trims, and re-splits reads per sample.
- `RNAseq_Evidence_QC_and_Trimming/fastp_leaves_RNAseq_evaluation.sbatch` — fastp QC-report-only run (all trimming disabled) on the leaf RNA-seq samples; merges each sample's lanes and generates per-sample HTML/JSON reports.
- `RNAseq_Evidence_QC_and_Trimming/fastp_roots_RNAseq_evaluation.sbatch` — same QC-report-only run for the root RNA-seq samples.

## 01_Hifiasm_Assembly — primary haplotype-resolved assembly

- `01_HyPR01RAW_hifiasm_HiC_script.sbatch` — hifiasm assembly integrating Hi-C reads (`--h1`/`--h2`) for haplotype phasing; produces the `hic.hap1`/`hic.hap2` contig graphs used downstream.
- `Initial_Assembly_QUAST/01_QUAST_HyPR01_hap1_script.sbatch` / `02_QUAST_HyPR01_hap2_script.sbatch` — converts each haplotype's hifiasm GFA to FASTA (gfatools) and runs QUAST on the raw, pre-scaffolding haplotigs; this is the earliest QC checkpoint in the pipeline, run before any YAHS scaffolding.

## 02_Scaffolding_YAHS — Hi-C scaffolding and manual curation

YAHS + 3D-DNA/Juicer scaffolding, with Hi-C contact maps manually reviewed
and corrected in Juicebox.

- `scripts/01_YAHS_scaffolding_HyPR01.sh` — telomere-aware, idempotent YAHS + bwa-mem2 (AVX2) Hi-C scaffolding run on the hifiasm hap1/hap2 contigs; auto-detects the YAHS `juicer`/`juicer_pre` converter, checkpoints alignment/indexing steps, and scaffolds using the `CCCTAAA` telomere motif.
- `scripts/02_yahs_juicer_pre_assembly_mode.sh` — YAHS scaffolding in Juicer "assembly mode," generating per-haplotype JBAT files for Juicebox review.
- `scripts/03_juicer_pre_MANUALFIX_part2.sh` — regenerates Juicer `.hic`/`.assembly` files after a manual Juicebox fix, using embedded Perl helper scripts (identity AGP, etc.).

## 03_Curation_RagTag — reference-guided curation and validation

- `scripts/01_PreRagTag_QC_Eudicots.sh` — pre-curation BUSCO (eudicots_odb10) QC on the YAHS scaffolds.
- `scripts/02_Length_Filter_100Kb.sh` — filters scaffolds to ≥100 kb (embedded Python streaming filter) ahead of RagTag.
- `scripts/03_RagTag_Curation.sh` — bidirectional RagTag scaffolding/curation between the two haplotypes (minimap2-based).
- `scripts/04_Define_Chromosomes.sh` — embedded Python engine that ranks scaffolds by length, assigns final chromosome IDs, and inserts the 10 kb N-spacer padding for the "UnChr" bin.
- `scripts/05_Rename_And_Verify_Hap2.sh` — renames/re-orders Hap2 chromosomes to match Hap1's numbering and verifies the mapping.
- `scripts/06_HiC_QC.sh` — Hi-C contact-map validation of the curated assembly (bwa-mem2 + pairtools).
- `scripts/07_PostCuration_Synteny_Assessment.sh` — whole-genome Hap1-vs-Hap2 minimap2 alignment + dotplot for a curation sanity check.
- `Validation_Dotplots/Pre_RagTag_Curation/pre_ragtag_HyPR01_dotplot.py` — synteny dot-plot from the pre-curation Hap1-vs-Hap2 PAF.
- `Validation_Dotplots/Post_RagTag_Curation/Post_RagTag_Synteny_plot.py` — same dot-plot, post-curation, for before/after comparison.
- `Validation_Dotplots/Post_RagTag_Curation/Merqury_HyPR01_barchart.py` — per-chromosome Merqury QV bar chart for both haplotypes.

## 04_Gap_Closure_and_Final_QC — long-read gap closing and whole-assembly QC

- `Gap_Closing/scripts/01_PacBio_HiFi_GapClosing.sh` — TGS-GapCloser run using the raw HiFi reads to close remaining assembly gaps.
- `Final_QC/scripts/02_Final_QC_Pipeline.sh` — final QUAST + BUSCO + Merqury QC pass on the gap-closed assembly.
- `Final_QC/Merqury_Analysis/03_Run_Matplotlib_Interpretation.sh` — per-chromosome Merqury QV bar chart (final-assembly version of the plot in stage 03); a self-contained SLURM wrapper that writes out and runs the plotting script (the former standalone `Matplotlib_Merqury_Interpretation.py` is now embedded via heredoc rather than kept as a separate file). Parses per-chromosome QV directly from Merqury's own `.qv` output files (`HyPR01_Final_QC.Hap{1,2}_curated_chr_un_gapclosed.qv`, expected alongside the script) rather than hand-typed values.

## 05_Annotations — repeat/TE masking, gene prediction, functional annotation

**EDTA — transposable-element discovery and visualization**
- `EDTA/scripts/01_Run_EDTA_Both_Haps.sh` — EDTA transposable-element discovery + LAI scoring, run in parallel per haplotype.
- `EDTA/EDTA_circos_plot/HyPR01_EDTA_Circos_Plot_Hap1_and_Hap2.sh` — computes TE-density and gene/TE collision tracks (bedtools) for both haplotypes and writes the R circos plotting scripts below.
- `EDTA/EDTA_circos_plot/Hap1/plot_circos.R` — renders the Hap1 TE-landscape circos plot (circlize).
- `EDTA/EDTA_circos_plot/Hap2/plot_circos.R` — renders the Hap2 TE-landscape circos plot.

**Masking — soft-masking + NUMT/NUPT hard-masking** (sequential steps: soft-mask everything EDTA found, then identify and hard-mask organellar-derived insertions in the nuclear genome)
NOTE: This section connects to the data generated in section 06_Organellar_Genomes for the organellar assembly
- `Masking/Softmasking/scripts/01_Generate_Softmasked_TEs_Genomes.sh` — generates the soft-masked genome FASTAs from the EDTA output.
- `Masking/NUMT_NUPT_Detection/Strict/scripts/02_Detect_NUMT_NUPT_HyPR01_BothHaps_STRICT.sh` — strict-mode minimap2 detection of NUMT/NUPT sequence using in-house HyPR01 organellar references.
- `Masking/NUMT_NUPT_Detection/Relaxed/scripts/03_Detect_NUMT_NUPT_HyPR01_BothHaps_RELAXED.sh` — relaxed-sensitivity minimap2 detection of NUMT/NUPT sequence using in-house HyPR01 organellar references.
- `Masking/NUMT_NUPT_Detection/Relaxed_K19/scripts/04_Detect_NUMT_NUPT_Relaxed_HyPR01_REF_K19.sh` — control track with inverted minimap2 seeding parameters (k=19, w=10) for a sensitivity comparison.
- `Masking/NUMT_NUPT_Detection/scripts/05_Plot_Organelle_Synteny.sh` — per-contig orientation/synteny plot between HyPR01 and NCBI organellar references.
- `Masking/NUMT_NUPT_Detection/scripts/06_NUMT_NUPT_density_mapping.sh` — genome-wide NUMT/NUPT density mapping and plotting.
- `Masking/NUMT_NUPT_Detection/Hard_Masked_Final/scripts/07_Hard_Masking_NUMT_NUPT.sh` — final hard-masking of the detected NUMT/NUPT regions in both haplotypes.

**Gene_density_vs_TE_density — combined gene/TE annotation visualization**
- `Gene_density_vs_TE_density/HyPR01_Local_Gene_vs_TE_Density.sh` — bins gene and TE annotations into 500 kb windows per haplotype for density co-plotting.
- `Gene_density_vs_TE_density/scripts/render_density_local.R` — renders the paired gene-density/TE-density tracks.

**Gene_Annotations / Evidence_RNAseq_STAR**
- `Gene_Annotations/Evidence_RNAseq_STAR/scripts/01_STAR_indexer.sh` — builds the combined (nuclear + organellar) STAR genome index.
- `Gene_Annotations/Evidence_RNAseq_STAR/scripts/02_STAR_aligner.sh` — competitive alignment of RNA-seq reads (2022 + newer sets) against the combined index.
- `Gene_Annotations/Evidence_RNAseq_STAR/scripts/03_BAM_header_purge.sh` — splits/fixes BAM headers so alignments are haplotype-specific.
- `Gene_Annotations/Evidence_RNAseq_STAR/scripts/04_BAM_evidence_audit.sh` — QC report on the final per-haplotype RNA-seq BAMs.

**Gene_Annotations / Gene_Prediction**
- `Gene_Annotations/Gene_Prediction/ANNEVO/scripts/05_run_ANNEVO_tracks.sh` — ANNEVO (deep-learning gene predictor) run on both haplotypes as a second, independent prediction track.
- `Gene_Annotations/Gene_Prediction/BRAKER3/scripts/06_BRAKER3_Hap1_CDS_ONLY.sh` — BRAKER3 gene prediction (RNA-seq + protein evidence) for Hap1.
- `Gene_Annotations/Gene_Prediction/BRAKER3/scripts/07_BRAKER3_Hap2_CDS_ONLY.sh` — same for Hap2.
- `Gene_Annotations/Gene_Prediction/BRAKER3/scripts/08_ExtractAugustusSeq_HyPR01_arabidopsis_hap2_script.sbatch` — extracts protein/CDS/transcript FASTAs from the AUGUSTUS-only subset of a BRAKER3 GFF3.

**Gene_Annotations / GFF3_Harmonization**
- `Gene_Annotations/GFF3_Harmonization/scripts/09_GFF3_Harmonization.sh` — merges/harmonizes the BRAKER3 and ANNEVO gene sets into one consensus GFF3 per haplotype.
- `Gene_Annotations/GFF3_Harmonization/scripts/10_PASA_UTR_Annotation.sh` — PASA UTR annotation on top of the harmonized gene models (via containerized PASApipeline).
- `Gene_Annotations/GFF3_Harmonization/GFF3_Cleanup_process/restore_haplotype_tiers.py` — restores original BRAKER/ANNEVO confidence-tier annotations lost during PASA processing.
- `Gene_Annotations/GFF3_Harmonization/scripts/GFF3toGTF_and_Stats.sbatch` — converts the AUGUSTUS-only GFF3 to GTF and runs GenomeTools (`gt stat`) summary statistics.
- `Gene_Annotations/GFF3_Harmonization/scripts/run_gt_stats_python_reformatted_gff3.sbatch` — same GenomeTools stats pass, on the Python-reformatted GFF3.

**Gene_Annotations / Functional_Annotation**
- `Gene_Annotations/Functional_Annotation/EggNOG/HyPR01_Hap1_Hap2_EGGNOG.sbatch` — eggNOG-mapper functional annotation for both haplotypes, merged back into GFF3 attributes.
- `Gene_Annotations/Functional_Annotation/InterProScan/HyPR01_Hap1_Hap2_Interproscan.sbatch` — InterProScan domain/family annotation for both haplotypes.

**Gene_Annotations / Final_Combined_GFF3**
- `Gene_Annotations/Final_Combined_GFF3/11_Annotation_Merge_and_Rename_combined.sbatch` — canonical final-merge pipeline: collapses redundant transcript models, removes internal-stop/missing-start transcripts, merges eggNOG + InterProScan into GFF3 column 9, renames/reindexes locus tags, recomputes CDS phase, does a post-hoc verification pass, purges ChrUn scaffolds into separate submission FASTAs/GFF3s with a full audit trail, and runs EMBLmyGFF3 to produce ENA `.embl.gz` flatfiles.

## 06_Organellar_Genomes — chloroplast and mitochondrial genomes

**Assembly**
- `Assembly/TIPPo/scripts/HyPR01_TIPPo.sbatch` — TIPPo organellar assembly (downsampled HiFi reads → chloroplast + mitochondrial contigs); the primary organellar assembly method used.
- `Assembly/HiMT/HyPR01_HiMT_organellar_assembly.sbatch` — HiMT organellar assembly directly from raw HiFi reads (independent method, for cross-validation against TIPPo).
- `Assembly/Assembly_Method_Comparison/scripts/build_organelle_summary_tables.py` — builds publication summary tables (size, GC, gene counts) for both organelle genomes from the annotation files directly (no hand-typed numbers); covers both TIPPo and HiMT output in one script.
- `Assembly/Assembly_Method_Comparison/scripts/run_table4_organelle_stats.sbatch` — environment wrapper (modules + a one-time biopython venv) for the script below.
- `Assembly/Assembly_Method_Comparison/scripts/build_table4_organelle_assembly_stats.py` — computes the TIPPo-vs-HiMT assembly-statistics table (reads used, N50/N90, assembly length/fragments, and an independently-remapped mean coverage) directly from the classified-reads and assembly files, rather than from an assembler's self-reported summary line. For HiMT, which pools both organelles' candidate reads into one file and never reports an organelle-specific count, the mitochondrion-specific read set is recovered by re-aligning against a combined mitochondrion+chloroplast reference and classifying by best hit. **Fully verified against a real run, all three assemblies.**

**Annotation / Chloroplast**
- `Annotation/Chloroplast/scripts/build_rotated_chloroplast_genome.py` — rotates the curated chloroplast genome to the conventional trnH-GUG/LSC start boundary.

**Annotation / Mitochondria**
- `Annotation/Mitochondria/scripts/01_setup_ogdraw_env.sh` — builds a local `ogdraw_env` conda environment running the real OGDraw (GeneMap) Perl codebase, for local/offline mitochondrial map rendering.
- `Annotation/Mitochondria/scripts/02_render_flattened_with_gc.pl` — OGDraw wrapper that enables the native GC-content ring (not exposed by OGDraw's stock CLI).
- `Annotation/Mitochondria/scripts/03_render_mito_contigs_ogdraw.sh` — renders each of the 3 mitochondrial contigs as an individual OGDraw circular map.
- `Annotation/Mitochondria/scripts/04_plot_mito_contigs_pycirclize.py` — alternative circular mitochondrial maps rendered with pyCirclize.
- `Annotation/Mitochondria/scripts/05_build_combined_map_with_contig_ring.py` — overlays a contig-boundary ring onto the combined/flattened OGDraw mitochondrial map.
- `Annotation/Mitochondria/scripts/run_mito_map_pipeline.sh` — orchestrator that runs the mitochondrial-map stages above in order, skipping completed steps.

---

## 07_Comparative_Analysis — genome-wide QC, coverage, synteny, and cross-species comparison

**Genome BUSCO completeness**
- `Genome_BUSCO/HyPR01_Final_comparison/01_BUSCO_HyPR01_Genome_vs_Protein_comprehesive_comparison.sbatch` — final-assembly BUSCO in both genome and protein mode, across 4 lineage databases.
- `Genome_BUSCO/HyPR01_Final_comparison/generate_final_plots.py` — aggregates and plots the above BUSCO results.
- `Genome_BUSCO/HyPR01_Final_comparison_noChrUn/02_BUSCO_HyPR01_Genome_vs_Protein_comprehesive_comparison_without_ChrUn.sbatch` — same BUSCO comparison, repeated with unplaced ("ChrUn") sequence excluded.
- `Genome_BUSCO/HyPR01_Final_comparison_noChrUn/generate_final_plots.py` — plots for the noChrUn BUSCO run.
- `Genome_BUSCO/scripts/03_BUSCO_HyPR01_bothhap_script.sbatch` — earliest BUSCO pass, on the raw hifiasm both-haplotype contig graph (converted from GFA).

**PacBio HiFi coverage**
- `Coverage_Map/01_Pacbio_Hifi_Coverage_Map_of_HyPR01_and_plotting.sbatch` — consolidated v10 pipeline: mapping + mosdepth + plotting in one job.
- `Coverage_Map/generate_coverage_v10.py` — the v10 plotting/QC-table engine invoked by the script above.
- `Coverage_Map/02_Pacbio_Hifi_Coverage_Map_plotting.sbatch` — re-plot of the v10 results with centered figure titles.
- `Coverage_Map/replot_centered_titles/generate_coverage_v10_centered.py` — the centered-titles variant of the v10 plotting engine.

**Nuclear organizing region / rDNA (local macOS analysis)**
- `Nuclear_Organizing_Region/01_run_NOR_confirmation_pipeline.sh` — confirms the Chr08 read-depth anomaly is a 45S rDNA/NOR via 4 independent lines of evidence (barrnap-HGV, BLASTn, nucmer self-dotplot, TRF).
- `Nuclear_Organizing_Region/02_run_genomewide_repeat_scan.sh` — extends that scan genome-wide (all chromosomes, both haplotypes) via barrnap-HGV + TRF.
- `Nuclear_Organizing_Region/03_run_blast_confirmation.sh` — same-species BLASTn confirmation for the Chr07 rDNA unit and the Chr01/Chr06 5S loci found by the genome-wide scan.
- `Nuclear_Organizing_Region/04_run_5S_depth_analysis.sh` — combines existing mosdepth coverage with barrnap tandem-copy counts to estimate absolute 5S rDNA copy number.

**Hi-C contact maps**
- `HiC_Visualization/Hap1_Hap2_Visualization_HiC_Data_v3.py` — Hi-C contact-map heatmaps (Juicebox color scheme) from `.mcool` files, with Mb-scale axis ticks and ICE-balanced/raw matrix options.

**Comparison against *H. perforatum* subsp. *chinense* (GWHEUWF00000000.1)**
- `Comparison_vs_Hp_chinense/Hypericum_perforatum_chinense_BUSCO/03_BUSCO_Hypericum_perforatum_subsp_chinense_Genome_comprehesive_comparison.sbatch` — BUSCO on the 4 chinense subgenomes (A–D) for comparison against HyPR01.
- `Comparison_vs_Hp_chinense/Hypericum_perforatum_chinense_BUSCO/generate_final_plots.py` — plots for the chinense subgenome BUSCO run.
- `Comparison_vs_Hp_chinense/Hypericum_perforatum_chinense_BUSCO/04_Telomere_Hypericum_perforatum_subsp_chinense_Genome_plot.sbatch` — telomeric-repeat detection/plot on the 4 chinense subgenomes; writes the R engine referenced below.
- `Comparison_vs_Hp_chinense/Hypericum_perforatum_chinense_Telomeres/scripts/render_telomeres_unified.R` — the telomere-scanning/plotting R engine for the subgenome comparison.
- `QUAST_Comparison_vs_Hp_chinense/01_QUAST_Comparison_HyPR01_and_Hypericum_tetraploid_chinense.sbatch` — QUAST comparison of HyPR01 against the tetraploid chinense assembly.
- `QUAST_Comparison_vs_Hp_chinense/02_QUAST_Comparison_HyPR01_and_Hypericum_tetraploid_chinense_v2.sbatch` — repeat QUAST run using the prepared subgenome-split FASTAs.

**Synteny circos plots**
- `Synteny_Circos/HyPR01_Hap1_Hap2_Synteny_Circos_plot_pipeline_v2.sh` — Hap1-vs-Hap2 synteny circos plot pipeline (MCScanX + Circos), mirrored hemispheres and balanced gaps.
- `Synteny_Circos/HyPR01_Hap1_Hap2_circos_plot/scripts/extract_haplotype_tracks.py` — extracts per-haplotype protein/BED tracks for the circos inputs.
- `Synteny_Circos/HyPR01_Hap1_Hap2_circos_plot/scripts/generate_circos_assets.py` — builds the karyotype/links Circos config files for the Hap1-vs-Hap2 plot.
- `Synteny_Circos/HyPR01_Hap2_1_vs_1_Synteny_Circos_plot_pipeline_v2.sh` — Hap2 all-vs-all pairwise master synteny pipeline (multiple reference species), with larger ideogram labels.
- `Synteny_Circos/HyPR01_Hap2_all_vs_all/scripts/build_pairwise_configs.py` — builds per-species-pair karyotype/link configs for the pairwise circos plots.
- `Synteny_Circos/Malpighiales_HyPR01_Hap2_Synteny_Circos_plot_pipeline_v1.sh` — Hap2-vs-Malpighiales (5 reference species) synteny circos pipeline; builds the multi-species karyotype/links Circos config itself via two inline `python3 -c` blocks (not a separate script — see the note on `build_circos_inputs.py` below).
- `Synteny_Circos/scripts/build_circos_inputs.py` — not used by the kept pipeline script above, which builds its karyotype/links config inline instead; do not treat this file as the source for Figure 16.
- `Synteny_Circos/scripts/extract_clean_proteins.py` — extracts/translates clean protein + BED tracks per species for MCScanX input.
- `Synteny_Circos/scripts/get_chrom_lengths.py` — records chromosome lengths per species from FASTA into a shared length table.
- `Synteny_Circos/scripts/rename_bed_chromosomes.py` — harmonizes chromosome naming across species' GFF3/BED files (e.g. `HpChr01`, `AtChr01`) and records the ID mapping.

**Comparative gene structure**
- `Gene_Structure_plotting/HyPR01_GeneStructure_MultiSpecies_v2.sh` — comparative gene-structure analysis (gene/CDS/exon/intron length distributions) across HyPR01 Hap2 and 5 reference species, plus a summary-statistics TSV export.
- `Gene_Structure_plotting/scripts/parse_gff3.R` — parses gene/CDS/exon features from a GFF3 into per-species metric tables.
- `Gene_Structure_plotting/scripts/plot_gene_structure.R` — renders the multi-species comparative gene-structure figures (ggplot2/patchwork).

## 08_Manuscript_Figures — remaining publication figures

All scripts in this stage were run locally on macOS, not on the HPC
cluster (see `SOFTWARE_VERSIONS.md` for the local-environment notes).

- `Telomeric_repeats_plot/HyPR01_Telomeric_Repeats_Plot_v3.sh` — scans both haplotypes for telomeric repeat motifs and renders per-chromosome ideograms in a faceted publication layout.
- `Telomeric_repeats_plot/scripts/render_telomeres_unified.R` — the telomere-scanning/plotting R engine written out by the wrapper script above.

---

## A note on the generalized paths

Every script originally contained absolute filesystem paths pointing at
a SLURM HPC cluster or a personal Mac (usernames, personal
cloud-storage sync paths, institution-specific mount points). Those have
been replaced with:

- `/path/to/input/directory` / `/path/to/input/directory/<filename>` — data this script reads, produced by an earlier stage.
- `/path/to/output/directory` / `/path/to/output/directory/<filename>` — data this script writes.
- `/path/to/your/directory` — a shared project root used for both reading and writing in the same script, or a local tool/environment install location (e.g. a conda env, a downloaded software package).
- `#SBATCH --output=`/`--error=` were changed to relative `logs/...` paths (works unedited under any SLURM submission directory).
- `#SBATCH --mail-user=`/`--mail-type=` lines were removed.
