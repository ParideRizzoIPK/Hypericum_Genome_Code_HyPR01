#!/bin/bash
# =============================================================================
# run_mito_map_pipeline.sh  --  end-to-end driver for the HyPR01 mitochondrial
# combined gene-map figure (OGDraw base + native GC graph + contig ring/legend).
#
# This is a thin, idempotent ORCHESTRATOR. The real work lives in the stage
# scripts in this directory; this just runs them in the correct order, skips
# work whose output already exists (unless --force), and fails fast.
#
# Stages:
#   1  ensure the ogdraw_env mamba environment            (setup_ogdraw_env.sh)
#   2  species-fix GenBank copies -> OGDraw_local_maps/fixed_gbk/   (sed)
#   3  render base map WITH native GC graph              (render_flattened_with_gc.pl)
#   4  [opt] per-contig OGDraw maps + montage           (render_mito_contigs_ogdraw.sh)
#   5  [opt] pyCirclize alternative maps                (plot_mito_contigs_pycirclize.py)
#   6  verify coordinates + overlay contig ring/legend  (build_combined_map_with_contig_ring.py)
#
# Usage:
#   ./run_mito_map_pipeline.sh [options]
#     --force            re-run every stage even if its output exists
#     --with-per-contig  also run stage 4 (3 separate per-contig maps + montage)
#     --with-pycirclize  also run stage 5 (pyCirclize alternative maps)
#     --skip-setup       skip stage 1 (assume ogdraw_env already exists)
#     --only <n>         run only stage <n> (1-6)
#     -h | --help        show this help
#
#   GENEMAP_DIR=/path/to/GeneMap-1.1.1 may be set in the environment to point
#   at a non-default OGDraw install.
# =============================================================================
set -Eeuo pipefail

# ----------------------------------------------------------------------------
# Paths (resolved relative to this script so the pipeline is relocatable)
# ----------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MITO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUT_DIR="${MITO_DIR}/OGDraw_local_maps"
FIXED_GBK="${OUT_DIR}/fixed_gbk"
PMGA_DIR="${MITO_DIR}/PMGA_Annotation_of_Mitochondria"
PMGA_ANNOT="${PMGA_DIR}/01.Annotation"
PMGA_JOB_ID="20250917150619818593"
GENEMAP_DIR="${GENEMAP_DIR:-/path/to/your/directory/GeneMap-1.1.1}"
OGDRAW_ENV="ogdraw_env"
PYC_ENV="pycirclize_env"
SPECIES="Hypericum perforatum"

GC_BASE="${OUT_DIR}/combined_flattened_ogdraw_GC.png"
FINAL_PNG="${OUT_DIR}/mito_combined_flattened_with_contig_ring_FINAL.png"

# ----------------------------------------------------------------------------
# Options
# ----------------------------------------------------------------------------
FORCE=0; WITH_PER_CONTIG=0; WITH_PYC=0; SKIP_SETUP=0; ONLY=""
while [ $# -gt 0 ]; do
  case "$1" in
    --force)           FORCE=1 ;;
    --with-per-contig) WITH_PER_CONTIG=1 ;;
    --with-pycirclize) WITH_PYC=1 ;;
    --skip-setup)      SKIP_SETUP=1 ;;
    --only)            ONLY="${2:-}"; shift ;;
    -h|--help)         sed -n '2,40p' "$0"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
  shift
done

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
log()  { printf '\n\033[1;34m==== %s\033[0m\n' "$*"; }
info() { printf '     %s\n' "$*"; }
die()  { printf '\033[1;31mERROR: %s\033[0m\n' "$*" >&2; exit 1; }
run_stage() { [ -z "$ONLY" ] || [ "$ONLY" = "$1" ]; }   # true if stage should run
# fresh <output-file>: true if we must (re)build it
fresh() { [ "$FORCE" = 1 ] || [ ! -s "$1" ]; }

mamba_has_env() { mamba env list 2>/dev/null | awk '{print $1}' | grep -qx "$1"; }

# ----------------------------------------------------------------------------
# Preflight
# ----------------------------------------------------------------------------
log "Preflight checks"
command -v mamba >/dev/null 2>&1 || die "mamba not found on PATH"
[ -d "$GENEMAP_DIR" ] || die "GeneMap-1.1.1 not found at: $GENEMAP_DIR (set GENEMAP_DIR=...)"
[ -f "${PMGA_DIR}/${PMGA_JOB_ID}.gbf" ] || die "PMGA flattened record missing: ${PMGA_DIR}/${PMGA_JOB_ID}.gbf"
for c in contig1 contig2 contig3; do
  [ -f "${PMGA_ANNOT}/${PMGA_JOB_ID}.${c}.gb" ] || die "PMGA per-contig record missing: ${c}"
done
mkdir -p "$OUT_DIR" "$FIXED_GBK"
info "MITO_DIR    = $MITO_DIR"
info "GENEMAP_DIR = $GENEMAP_DIR"
info "OUT_DIR     = $OUT_DIR"

# ----------------------------------------------------------------------------
# Stage 1 - ensure ogdraw_env
# ----------------------------------------------------------------------------
if run_stage 1 && [ "$SKIP_SETUP" = 0 ]; then
  log "Stage 1: ensure '$OGDRAW_ENV' mamba environment"
  if mamba_has_env "$OGDRAW_ENV" && [ "$FORCE" = 0 ]; then
    info "'$OGDRAW_ENV' already exists -- skipping (use --force to rebuild)"
  else
    info "building '$OGDRAW_ENV' via setup_ogdraw_env.sh ..."
    bash "${SCRIPT_DIR}/setup_ogdraw_env.sh"
  fi
fi
# pycirclize_env is required from stage 5/6 onward; check it exists
mamba_has_env "$PYC_ENV" || die "mamba env '$PYC_ENV' not found (needed for the Python stages)"

# ----------------------------------------------------------------------------
# Stage 2 - species-fix GenBank copies (YOUR_SPECIES -> Hypericum perforatum)
#           on COPIES only; PMGA originals are never modified.
# ----------------------------------------------------------------------------
if run_stage 2; then
  log "Stage 2: species-fixed GenBank copies -> ${FIXED_GBK}"
  # flattened record (the one OGDraw renders for the combined map)
  if fresh "${FIXED_GBK}/combined_flattened.gb"; then
    sed "s/YOUR_SPECIES/${SPECIES}/g" "${PMGA_DIR}/${PMGA_JOB_ID}.gbf" \
        > "${FIXED_GBK}/combined_flattened.gb"
    info "wrote fixed_gbk/combined_flattened.gb"
  else
    info "fixed_gbk/combined_flattened.gb exists -- skipping"
  fi
  # per-contig records (used for verification lengths + per-contig renders)
  for c in contig1 contig2 contig3; do
    if fresh "${FIXED_GBK}/${c}.gb"; then
      sed "s/YOUR_SPECIES/${SPECIES}/g" "${PMGA_ANNOT}/${PMGA_JOB_ID}.${c}.gb" \
          > "${FIXED_GBK}/${c}.gb"
      info "wrote fixed_gbk/${c}.gb"
    else
      info "fixed_gbk/${c}.gb exists -- skipping"
    fi
  done
fi

# ----------------------------------------------------------------------------
# Stage 3 - render base map WITH OGDraw's native GC graph
# ----------------------------------------------------------------------------
if run_stage 3; then
  log "Stage 3: OGDraw base map WITH native GC graph"
  if fresh "$GC_BASE"; then
    ( cd "$GENEMAP_DIR" && \
      mamba run -n "$OGDRAW_ENV" perl -Ilib \
        "${SCRIPT_DIR}/render_flattened_with_gc.pl" \
        "${FIXED_GBK}/combined_flattened.gb" \
        "$GC_BASE" ) 2>&1 | grep -viE "deprecated|redefined|uninitialized" || true
    [ -s "$GC_BASE" ] || die "GC base image was not produced: $GC_BASE"
    info "wrote $(basename "$GC_BASE")"
  else
    info "$(basename "$GC_BASE") exists -- skipping (use --force to re-render)"
  fi
fi

# ----------------------------------------------------------------------------
# Stage 4 (optional) - per-contig OGDraw maps + 3-panel montage
# ----------------------------------------------------------------------------
if run_stage 4 && { [ "$WITH_PER_CONTIG" = 1 ] || [ "$ONLY" = 4 ]; }; then
  log "Stage 4: per-contig OGDraw maps + montage"
  bash "${SCRIPT_DIR}/render_mito_contigs_ogdraw.sh"
fi

# ----------------------------------------------------------------------------
# Stage 5 (optional) - pyCirclize alternative maps
# ----------------------------------------------------------------------------
if run_stage 5 && { [ "$WITH_PYC" = 1 ] || [ "$ONLY" = 5 ]; }; then
  log "Stage 5: pyCirclize alternative maps"
  mamba run -n "$PYC_ENV" python3 "${SCRIPT_DIR}/plot_mito_contigs_pycirclize.py"
fi

# ----------------------------------------------------------------------------
# Stage 6 - verify coordinates + overlay contig ring/legend -> FINAL PNG + SVG
# ----------------------------------------------------------------------------
if run_stage 6; then
  log "Stage 6: verify + overlay contig ring/legend (FINAL PNG + SVG)"
  [ -s "$GC_BASE" ] || die "GC base image missing; run stage 3 first ($GC_BASE)"
  if fresh "$FINAL_PNG"; then
    mamba run -n "$PYC_ENV" python3 \
      "${SCRIPT_DIR}/build_combined_map_with_contig_ring.py"
  else
    info "$(basename "$FINAL_PNG") exists -- skipping (use --force to rebuild)"
  fi
fi

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
log "Done. Key outputs in ${OUT_DIR}:"
for f in "combined_flattened_ogdraw_GC.png" \
         "mito_combined_flattened_with_contig_ring_FINAL.png" \
         "mito_combined_flattened_with_contig_ring_FINAL.svg" \
         "build_combined_map_with_contig_ring.log"; do
  if [ -s "${OUT_DIR}/${f}" ]; then
    printf '     [ok]   %s\n' "$f"
  else
    printf '     [--]   %s (not generated in this run)\n' "$f"
  fi
done
[ "$WITH_PER_CONTIG" = 1 ] && printf '     [ok]   contig{1,2,3}_ogdraw.png + mito_contigs_combined_ogdraw.png\n' || true
[ "$WITH_PYC" = 1 ]        && printf '     [ok]   pyCirclize_maps/*.png\n' || true
echo
