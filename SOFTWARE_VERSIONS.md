# HyPR01 Pipeline — Software & Dependencies

This is a manifest, not an installable environment file — five different
provisioning methods are used across the scripts in this collection, so one
combined `environment.yml` doesn't fit (see `PIPELINE_OVERVIEW.md` for the
per-script breakdown).

**Every version number below is taken directly from what's written in the
scripts** (a `module load tool/X.Y.Z` line, a container image path, a
directory name, a code comment), **except entries marked "confirmed
live"** (checked directly against the HPC modules / conda environments),
**"confirmed directly from the installed software itself"** (e.g.
GeneMark-ETP, Bandage), **or "confirmed by the author"** (e.g. YAHS's
standalone build) — none of which were stated in any script. The last
section lists everything with no version anywhere, by any of these routes.

Stage numbers below refer to `PIPELINE_OVERVIEW.md`'s folders (`00`–`08`,
all function-based).

## How software is provisioned across this pipeline

| Method | Where used |
|---|---|
| HPC environment modules (`module load tool/version`) | Nearly everything run on the SLURM cluster (stages 00–06, 07) |
| Cluster `micromamba`/`conda` named environments | NanoPlot QC (00), the synteny-circos pipelines (07), the EMBL submission tooling (05) |
| Apptainer/Singularity containers | EDTA, RagTag, PASA, ANNEVO, TIPPo, BRAKER3, BUSCO (all bound via `APPTAINER_BIND`/`SINGULARITY_BIND`) |
| Python `venv` + `pip install` at job runtime | BUSCO/coverage-map/synteny plotting steps (04, 07); organelle assembly-stats script (06, biopython) |
| Local macOS `mamba`/`conda` named environments | Organellar figure rendering (06), Nuclear Organizing Region analysis (07), EDTA-circos tool binaries (05), Hi-C contact-map visualization (07) |
| Standalone downloaded software (no module/container) | YAHS (in some script variants), GeneMark-ETP, GeneMap/OGDraw, barrnap-HGV, Juicer Tools |
| R package auto-install from CRAN at runtime (`install.packages(...)`, no version pin) | All local macOS R plotting scripts (05, 07, 08) |
| Local macOS, plain `python3` — no environment/package needed (pure standard library) | Chloroplast genome rotation script (06) |

---

## HPC environment modules (versioned in-script)

| Tool | Version | Stage(s) | Notes |
|---|---|---|---|
| ANNEVO | 2.3.1 | 05 | |
| augustus | 3.3.1 | 07 | earliest/legacy BUSCO script only |
| BRAKER | 3.0.8 | 05 | |
| BUSCO | 5.8.2 / 5.2.2 | 03, 04, 07 | 5.2.2 used in one early/legacy script only |
| bwa-mem2 | 2.2.1 | 02, 03 | |
| cutadapt | 3.3 | 00 | from a code comment ("default cutadapt/3.3"), not an explicit module-load pin |
| clc-assembly-cell | 5.1.1 | 00 | from a code comment ("default clc-assembly-cell/5.1.1"), not an explicit module-load pin |
| EDTA | 2.2.2 | 05 | |
| eggnog-mapper | 2.1.12 | 05 | |
| fastp | 1.0.1 | 00 | |
| fastqc | 0.12.1 | 00 | Hi-C raw-read QC |
| gfatools | 0.5-r234 | 01, 06, 07 | confirmed live; not stated in any script |
| genometools (`gt`) | 1.6.2 | 05 | confirmed live; not stated in any script |
| gffread | 0.12.6 | 05 | pinned in-script for 3 of 5 scripts (EggNOG, InterProScan, ANNEVO); confirmed live the other 2 load the same module |
| gmap | 2019-09-12 | 05 | |
| hifiasm | 0.25.0 | 01 | confirmed live; not stated in any script |
| interproscan | 5.57-90.0 | 05 | |
| kmc | 3.1.1 | 00 | confirmed live; not stated in any script |
| jdk / openjdk | 1.8.0_191 / 11 | 02 | 1.8.0_191 for the 3D-DNA post-review step, 11 for later YAHS variants — different tools, not version drift |
| lastz | 1.03.73 | 02 | |
| matplotlib (HPC module) | 3.7.1 | 03, 04, 05 | distinct from the pip-installed matplotlib listed separately below |
| merqury | 1.3 | 04 | |
| minimap2 | 2.24 | 02, 03, 05, 07 | |
| mosdepth | 0.2.6 | 07 | |
| pairtools | 1.0.3 | 03 | |
| parallel (GNU Parallel) | 20210222 | 02, 05, 07 | confirmed live; not stated in any script |
| python (HPC module) | 3.12.0 | throughout 01–07 | |
| quast | 5.2.0 | 01, 04, 07 | pinned in-script for 04/07; confirmed live the stage-01 pair loads the same module — an inconsistent comment there, not a real version difference |
| R | 4.5.1 | 07 | chinense telomere plot |
| RagTag | 2.1.0 | 03 | |
| samtools | 1.23.1 | throughout 02–07 | |
| STAR | 2.7.9a | 05 | |
| TGS-GapCloser | 1.2.1 | 04 | |
| TIPP (TIPPo) | 1.3.0 | 06 | |
| 3d-dna | 180922 | 02 | date-based version tag |
| yahs | 1.2 | 02 | module, most scripts; the kept `01_YAHS_scaffolding_HyPR01.sh` instead builds YAHS outside the module system — see the Standalone table below |

## Named conda/mamba/venv environments

| Environment | Where built | Packages named in-script | Version | Notes |
|---|---|---|---|---|
| `nanoplot` | cluster `micromamba` | NanoPlot | 1.47.1 | confirmed live |
| `synteny_circos` | cluster `micromamba` | Circos | 0.69.9 | confirmed live |
| `last_env` | cluster `micromamba` | LAST (the aligner) | 1651 | confirmed live; **not MCScanX** — see the MCScanX note below |
| `busco_plot_env` | cluster `python -m venv` + `pip install matplotlib numpy` | matplotlib, numpy | not pinned | a plain venv, not a micromamba env |
| `.venv_coverage_plotting` / `py_env` / `py_plot_env` | cluster/local `python -m venv` + `pip install pandas numpy matplotlib tabulate` | pandas, numpy, matplotlib, seaborn, tabulate | not pinned | plain venvs, not micromamba envs; `seaborn` is also used (see `generate_coverage_v10.py` imports) though not in the explicit `pip install` line |
| `emblmygff3` | cluster `micromamba` (`$MAMBA_ROOT_PREFIX/envs/emblmygff3/bin/EMBLmyGFF3`) | EMBLmyGFF3 | 2.4 | confirmed live |
| `ogdraw_env` | local macOS `mamba create` (`osx-64`/Rosetta 2) | perl-bioperl, perl-gd, perl-postscript-simple, imagemagick, clang_osx-64/clangxx_osx-64 | not pinned | osx-64/Rosetta 2 build chosen because arm64 builds aren't available for this exact BioPerl/GD combination; bioconda/conda-forge solve at build time |
| `pycirclize_env` | local macOS `mamba` | python=3.11, pycirclize, biopython/Pillow (via matplotlib's dependency chain) | python 3.11 (pinned); pycirclize/biopython not pinned | |
| `nor_test` | local macOS `conda` | barrnap-HGV, BLAST+, MUMmer (nucmer), TRF | not pinned | |
| `bedtools_env`, `samtools-env` | local macOS `miniforge3` | bedtools, samtools | not pinned | referenced via `$PATH`, not activated in-script |
| `hictk_vis` | local macOS `miniforge3` mamba env | python, hictk, cooler, matplotlib, numpy | python 3.12.13, hictk 2.2.0, cooler 0.10.4, matplotlib 3.10.9, numpy 2.4.6 | confirmed live (`mamba list -n hictk_vis`); activated manually before running `Hap1_Hap2_Visualization_HiC_Data_v3.py` (not activated in-script); also provides the `hictk` CLI used for the one-off `.hic`→`.mcool` conversion ahead of plotting |
| R user library `~/Library/R/hypr01-lib` | local macOS, auto-installed by each R script on first run | ggplot2, dplyr, tidyr, patchwork, svglite, data.table, scales, ggnewscale, circlize | not pinned | installed from whatever is current on CRAN at run time |

**MCScanX, not resolved:** the Synteny_Circos pipelines (07) reference
`synteny_circos` and `last_env`. Neither env contains an MCScanX-named
conda package, and `last_env` actually holds the LAST aligner, not
MCScanX. If installed, it's likely a standalone binary elsewhere on
`$PATH`, not a conda package.

## Containers (Apptainer/Singularity)

Referenced via `APPTAINER_BIND`/`SINGULARITY_BIND`; the container engine's
own version is never stated in any script.

| Tool | Container/image reference | Version |
|---|---|---|
| PASApipeline | `/opt/Bio/PASApipeline/2.5.2/bin/pasapipeline.sif` | 2.5.2 |
| ANNEVO | `/opt/Bio/ANNEVO/2.3.1/bin/container` | 2.3.1 |
| EDTA, RagTag, BRAKER3, TIPPo, BUSCO | run through the `EDTA`/`RagTag`/`BRAKER`/`TIPP`/`BUSCO` HPC modules (see module table above) with Apptainer bind paths | as listed above |

## Standalone / downloaded software (no module or container)

| Tool | Version | Notes |
|---|---|---|
| Juicer Tools | 1.22.01 | `juicer_tools_1.22.01.jar`, used directly (not a module) |
| YAHS | 1.2.2 | Most scripts instead use the `yahs` module (version 1.2 — see the HPC modules table); the kept `01_YAHS_scaffolding_HyPR01.sh` builds YAHS outside the module system via a standalone `YAHS_ROOT` instead — [github.com/c-zhou/yahs](https://github.com/c-zhou/yahs) release 1.2.2 (confirmed by the author; not stated in-script, where it's only referenced via the sanitized `${YAHS_ROOT}` path). **This is a genuine version difference from the module's 1.2, not a typo** — the standalone build used for the final scaffolding run is one patch version ahead of what the module provides. |
| GeneMark-ETP | 1.02 | shipped as `gmetp_linux_64.tar.gz` + a license key file (`gm_key.gz`), required by BRAKER3; version confirmed directly from the bundle (`bin/gmetp.pl`'s `my $version = "1.02"`, printed as "ETP version 1.02" at runtime; files dated 2023-11-07). The bundled GeneMark-ES Suite component in `bin/gmes/` is a separate, unrelated version (4.x). |
| GeneMap / OGDraw | 1.1.1 | `GeneMap-1.1.1`, downloaded from https://chlorobox.mpimp-golm.mpg.de/OGDraw-Downloads.html |
| barrnap-HGV | not stated | a specific fork/build (not stock barrnap), referenced by path |
| phyloFlash | not stated | a dependency of the local `nor_test` environment; the tool's own source tree was excluded from this scripts collection as a bundled third-party download (see `PATH_PII_AUDIT.md`) |
| TRF (Tandem Repeat Finder) | not stated | via `nor_test` env |
| BLASTn / makeblastdb | not stated | via `nor_test` env |
| MUMmer (nucmer) | not stated | via `nor_test` env |
| Bandage | 0.8.1 | manual GUI step, not scripted (confirmed from the app's own About dialog); used once to manually resolve the chloroplast inverted-repeat graph (see the chloroplast docstrings in stage 06) |

## Python packages (pip-installed at job runtime, not version-pinned)

pandas, numpy, matplotlib, seaborn, tabulate, biopython (`Bio.SeqIO`),
Pillow (`PIL`), pycirclize, argparse (stdlib). (`cooler`, along with the
`matplotlib`/`numpy` used by the Hi-C visualization script specifically,
is **not** in this bucket — see the `hictk_vis` row above; it's pinned via
a local conda env, not installed unpinned at job runtime.)

## R packages (CRAN, auto-installed at runtime, not version-pinned)

ggplot2, dplyr, tidyr, patchwork, svglite, data.table, scales, ggnewscale,
circlize.

## Perl modules (bundled with the local `ogdraw_env` / GeneMap-1.1.1 install, not version-pinned)

BioPerl (`Bio::Perl`), `GD`, `GeneMap::Plastome`, `PostScript::Simple`.

## A note on `PROJECT_ID`

The `Final_Combined_GFF3` scripts (05) set `PROJECT_ID="PRJEB123953"` — this
is a real, public ENA/EBI study accession, kept intentionally as submission
metadata (same treatment as the `AUTHORS=` field elsewhere: legitimate
provenance data, not something to scrub).

---

## Everything not version-pinned in the scripts

Split by which machine actually ran it, checked against each script's own
`#SBATCH` directives — HPC entries below check via `module avail`/`module
spider`; local macOS entries check via `mamba`/`conda list -n <env>`,
`brew list`, or the app's own version info.

### HPC (SLURM cluster)

- **MCScanX** — see the note above; not resolved as a conda package in
  either `synteny_circos` or `last_env`.
- **`busco_plot_env`, `.venv_coverage_plotting`/`py_env`/`py_plot_env`** —
  plain Python `venv`s, not micromamba envs.

### Local macOS

- **Named environments:** `ogdraw_env`, `pycirclize_env` (except its
  `python=3.11`, which is pinned), `nor_test`, `bedtools_env`,
  `samtools-env`.
- **Standalone tools:** barrnap-HGV, phyloFlash, TRF, BLASTn/makeblastdb,
  MUMmer/nucmer (all via `nor_test`).
- **Language-level packages:** every remaining Python pip package and
  every R CRAN package listed in the sections above this one.