#!/bin/bash
# ==============================================================================
# Render the 3 HyPR01 mitochondrial contigs as circular gene maps using the
# real, local OGDraw engine (GeneMap-1.1.1). Requires `ogdraw_env` -- see
# setup_ogdraw_env.sh in this same directory if it doesn't exist yet.
#
# Note: bin/drawgenemap's circular-map code path is hardcoded to instantiate
# GeneMap::Plastome regardless of organelle type (there is no CLI flag to
# select GeneMap::Chondriome instead). In practice this produced correct
# mitochondrial gene-category colors/legend anyway (categories are derived
# from the actual gene names in the GenBank file, not hardcoded per module),
# confirmed by inspecting the rendered output directly -- not assumed.
# ==============================================================================
set -e -u -o pipefail

GENEMAP_DIR="/path/to/your/directory/GeneMap-1.1.1"
SRC_GBK_DIR="/path/to/input/directory"
OUT_DIR="/path/to/output/directory"
FIXED_GBK_DIR="${OUT_DIR}/fixed_gbk"
PMGA_JOB_ID="20250917150619818593"

mkdir -p "$OUT_DIR" "$FIXED_GBK_DIR"

# PMGA's GenBank template leaves DEFINITION/ORGANISM/SOURCE as "YOUR_SPECIES"
# placeholder text -- fixed here on COPIES only; PMGA's originals are never
# modified. (This placeholder text was previously found baked directly into
# a rendered figure title before this fix was applied.)
echo "=== Fixing species placeholder text (copies only, originals untouched) ==="
for c in contig1 contig2 contig3; do
    sed -e 's/YOUR_SPECIES/Hypericum perforatum/g' \
        "${SRC_GBK_DIR}/${PMGA_JOB_ID}.${c}.gb" > "${FIXED_GBK_DIR}/${c}.gb"
done

echo "=== Rendering all 3 contigs ==="
cd "$GENEMAP_DIR"
for c in contig1 contig2 contig3; do
    echo "  $c"
    mamba run -n ogdraw_env perl -Ilib bin/drawgenemap \
        --infile "${FIXED_GBK_DIR}/${c}.gb" \
        --format png \
        --outfile "${OUT_DIR}/${c}_ogdraw.png" \
        --density 300x300 \
        2>&1 | grep -v "deprecated" || true
done

echo
echo "=== Building combined multi-panel figure (montage of the 3 separate renders) ==="
echo "    Deliberately NOT a single concatenated-sequence render -- see the"
echo "    module docstring above and Figure_Legends.md for why: joining the"
echo "    3 contigs into one artificial sequence would recreate the exact"
echo "    fabricated-adjacency problem already identified and rejected for"
echo "    PMGA's own flattened annotation (see the mitochondria audit doc, §8)."
mamba run -n ogdraw_env montage \
    -label "A" "${OUT_DIR}/contig1_ogdraw.png" \
    -label "B" "${OUT_DIR}/contig2_ogdraw.png" \
    -label "C" "${OUT_DIR}/contig3_ogdraw.png" \
    -tile 3x1 -geometry 1000x1000+20+20 -background white \
    -pointsize 48 -font Helvetica \
    "${OUT_DIR}/mito_contigs_combined_ogdraw.png"

echo
echo "Done. Output in: $OUT_DIR"
ls -la "$OUT_DIR"/*.png
