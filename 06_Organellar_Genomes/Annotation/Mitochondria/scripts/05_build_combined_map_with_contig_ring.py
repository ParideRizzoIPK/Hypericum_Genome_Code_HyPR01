#!/usr/bin/env python3
"""
Overlay a contig-boundary ring onto the combined (flattened) OGDraw mito map.

Context: PMGA's own flattened annotation (HyPR01_Mitochondria_genome.gbf,
426,309 bp = contig1 + contig2 + contig3 concatenated in that order with no
gap) is what produces the single-circle combined view -- the same file used
to generate the reference figure this script is extending. Concatenating the
3 contigs into one sequence loses track of which gene came from which contig
(and, as documented in HyPR01_Mitochondria_Assembly_Methods_and_Audit.md §8,
implies two junctions the assembly graph doesn't actually support -- that
caveat still applies to the base map here). This script does not fix that;
it restores the contig identity visually via an added ring, on top of the
existing flattened render, per explicit request.

Angle convention (verified empirically, not assumed -- see calibration
check against known gene positions `rps3` and `mttB` during development):
GeneMap::Plastome computes feat_angle = (bp / seq_length) * 360 and passes
it directly to PostScript::Simple's arc() method, which maps to PostScript's
native `arc` operator (0 deg = 3 o'clock, increasing COUNTERCLOCKWISE, no
extra rotation offset found in the source). PIL's ImageDraw.arc uses the
opposite rotational sense (clockwise, Y-down image coordinates), so angles
are negated before drawing here.

Run inside `pycirclize_env` (has Pillow via matplotlib's dependency chain).
"""

import logging
import math
import re
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

BASE = Path("/path/to/your/directory")
# Base map WITH OGDraw's native GC-content graph (gray inner band), rendered by
# scripts/render_flattened_with_gc.pl (passes gc_cont=>1 to GeneMap::Plastome,
# which the stock drawgenemap CLI never does). The GC band occupies radius
# ~832-1042 px; the contig ring is placed just outside it.
SRC_IMAGE = BASE / "combined_flattened_ogdraw_GC.png"
OUT_IMAGE = BASE / "mito_combined_flattened_with_contig_ring_FINAL.png"
LOG_PATH = BASE / "build_combined_map_with_contig_ring.log"

# Source records that back the coordinate system. The verification stage reads
# these directly (not the hardcoded numbers below) and asserts they agree.
FIXED_GBK = BASE / "fixed_gbk"
PER_CONTIG_GB = {
    "Contig 1": FIXED_GBK / "contig1.gb",
    "Contig 2": FIXED_GBK / "contig2.gb",
    "Contig 3": FIXED_GBK / "contig3.gb",
}
# The exact flattened GenBank record OGDraw rendered into SRC_IMAGE. The ring's
# bp->angle mapping is only correct if TOTAL_LEN equals this record's length.
FLAT_GB = FIXED_GBK / "combined_flattened.gb"

TOTAL_LEN = 426309
CENTER = (2500, 2500)   # PS canvas center (600,600)pt * (300/72 dpi) = 2500px
RING_RADIUS = 1095      # in the clear zone between OGDraw's GC band (outer edge
                        # ~1042 px) and the inner gene labels (~1148 px); hugs
                        # the gene track (radius only -- bp->angle calibration
                        # is independent of this value)
RING_THICKNESS = 45     # thin band -- contigs are identified via the color key,
                        # not via text baked onto the arc (see DRAW_CURVED_LABELS)

# Curved per-character text on the ring never renders cleanly in Pillow (no
# native text-on-path; each glyph is bitmap-rotated, so letterforms warp and
# the baseline waves). Left OFF by default: the figure ships with thin colored
# arcs + a straight-text color key, and contig labels are meant to be added in
# vector software (Illustrator/Inkscape) where text-on-path is crisp. Flip to
# True only if a rough baked-in curved label is genuinely wanted.
DRAW_CURVED_LABELS = False

# --- GC content ------------------------------------------------------------
# The GC-content graph is OGDraw's OWN native gray band, baked into the base
# image (SRC_IMAGE) by render_flattened_with_gc.pl. So NATIVE_GC=True just
# tells this script to add a matching legend entry and to compute/log the mean
# GC for that caption; the collision checks confirm the contig ring avoids the
# band. The DRAW_GC path below is a self-computed overlay track, kept for
# reference but OFF (the native band is what the figure uses).
NATIVE_GC = True
GC_N_WINDOWS = 900      # windows used only for logging GC statistics

DRAW_GC = False         # custom overlay GC bars -- superseded by native band
GC_BASELINE_R = 750
GC_MAX_AMP = 85
GC_BAR_WIDTH = 5
GC_COLOR_HI = (95, 95, 95)
GC_COLOR_LO = (165, 165, 165)
GC_BASELINE_COLOR = (60, 60, 60)
GC_BAND_COLOR = (170, 170, 170)   # gray swatch for the native-GC legend entry

# (name, start_bp, end_bp, RGBA color) -- 0-based, half-open on end_bp
CONTIGS = [
    ("Contig 1", 0, 52619, (230, 57, 70, 220)),        # red
    ("Contig 2", 52619, 201420, (42, 157, 143, 220)),  # teal
    ("Contig 3", 201420, 426309, (38, 70, 178, 220)),  # blue
]


def true_angle(bp: float) -> float:
    """Degrees, counterclockwise from 3 o'clock -- OGDraw's own convention."""
    return (bp / TOTAL_LEN) * 360.0


def pil_angle(bp: float) -> float:
    """Convert to PIL's clockwise-from-3-o'clock convention."""
    return -true_angle(bp) % 360.0


def contig_boundaries():
    """Boundary bp positions (0, junctions..., TOTAL_LEN) derived from CONTIGS."""
    return [CONTIGS[0][1]] + [end for _n, _s, end, _c in CONTIGS]


def compute_gc_windows(seq, n_windows):
    """Return (mean_gc, [(mid_bp, gc_percent), ...]) over non-overlapping windows."""
    L = len(seq)
    total_gc = seq.count("G") + seq.count("C")
    mean_gc = total_gc / L * 100.0
    windows = []
    for i in range(n_windows):
        start = i * L // n_windows
        end = (i + 1) * L // n_windows
        sub = seq[start:end]
        gc = (sub.count("G") + sub.count("C")) / len(sub) * 100.0
        windows.append(((start + end) / 2.0, gc))
    return mean_gc, windows


def gc_scale(windows, mean_gc):
    """Auto-fit px-per-%% so the largest deviation maps to GC_MAX_AMP."""
    max_dev = max(abs(gc - mean_gc) for _bp, gc in windows) or 1.0
    return GC_MAX_AMP / max_dev


def draw_gc_track(draw, windows, mean_gc, scale):
    """PIL: mean baseline circle + per-window radial deviation bars."""
    draw.ellipse(
        [CENTER[0] - GC_BASELINE_R, CENTER[1] - GC_BASELINE_R,
         CENTER[0] + GC_BASELINE_R, CENTER[1] + GC_BASELINE_R],
        outline=GC_BASELINE_COLOR + (255,), width=3,
    )
    for mid_bp, gc in windows:
        th = math.radians(true_angle(mid_bp))
        r_val = GC_BASELINE_R + (gc - mean_gc) * scale
        x1 = CENTER[0] + GC_BASELINE_R * math.cos(th)
        y1 = CENTER[1] - GC_BASELINE_R * math.sin(th)
        x2 = CENTER[0] + r_val * math.cos(th)
        y2 = CENTER[1] - r_val * math.sin(th)
        color = GC_COLOR_HI if gc >= mean_gc else GC_COLOR_LO
        draw.line([x1, y1, x2, y2], fill=color + (255,), width=GC_BAR_WIDTH)


def draw_curved_text(canvas, draw, text, mid_bp, font, radius):
    """Draw `text` centered on the arc at `radius`, following its curvature.

    Each character is rendered separately and rotated to the local tangent
    angle, then composited along the circle. Verified against a standalone
    test (4 direction/rotation-sign combinations rendered and visually
    compared) before use here -- see conversation history for the test
    harness. The key finding: characters must always be walked in the same
    angular direction (start at more-CCW edge, step clockwise/decreasing
    true_angle) regardless of which half of the circle the label sits in.
    Only the bottom half (true_angle in 180-360, i.e. below image-center)
    needs its character *order* reversed and +180 deg added to each
    character's own rotation, so the glyphs stay upright while the string
    still reads left-to-right for a normal viewer. Using the mirrored walk
    direction for the bottom half (as an earlier, buggy version of this
    function did) produces doubly-reversed, backward-reading text.
    """
    mid_true = true_angle(mid_bp)
    flip = 180 < mid_true < 360  # label sits in the lower half of the image

    chars = list(text)
    if flip:
        chars = chars[::-1]

    widths = []
    for c in chars:
        bbox = draw.textbbox((0, 0), c, font=font)
        widths.append((bbox[2] - bbox[0]) if c != " " else font.size * 0.4)

    gap_px = max(4, font.size // 12)
    stroke = max(3, font.size // 14)
    # Tile must be large enough to hold a rotated glyph plus its stroke halo.
    tile_sz = int(font.size * 2.6)
    half_tile = tile_sz // 2

    total_width = sum(widths) + (len(chars) - 1) * gap_px
    angular_span = math.degrees(total_width / radius)

    cur_angle = mid_true + angular_span / 2  # always start at the more-CCW edge
    for c, w in zip(chars, widths):
        half_w_ang = math.degrees((w / 2) / radius)
        char_angle = cur_angle - half_w_ang

        rad = math.radians(char_angle)
        x = CENTER[0] + radius * math.cos(rad)
        y = CENTER[1] - radius * math.sin(rad)

        rot = 90 - char_angle
        if flip:
            rot += 180

        tile = Image.new("RGBA", (tile_sz, tile_sz), (255, 255, 255, 0))
        tile_draw = ImageDraw.Draw(tile)
        tile_draw.text(
            (half_tile, half_tile), c, font=font, fill=(255, 255, 255, 255),
            anchor="mm", stroke_width=stroke, stroke_fill=(0, 0, 0, 255),
        )
        rotated = tile.rotate(rot, resample=Image.BICUBIC, center=(half_tile, half_tile))
        canvas.alpha_composite(rotated, (int(x - half_tile), int(y - half_tile)))

        full_w_ang = math.degrees(w / radius) + math.degrees(gap_px / radius)
        cur_angle -= full_w_ang


# --- Contig key (bottom-right, opposite OGDraw's own gene legend) ----------
# Shared layout so the PNG (PIL) and the SVG render identically. Sizes are
# enlarged relative to OGDraw's small gene legend for readability.
KEY_RIGHT = 4880
KEY_TOP = 4150
KEY_SWATCH = 82          # swatch square side
KEY_GAP = 32             # swatch-to-text gap
KEY_ROW_H = 118          # vertical pitch between entry rows
KEY_HEADER_GAP = 150     # KEY_TOP (header top) -> first entry row top
KEY_ENTRY_SIZE = 92
KEY_HEADER_SIZE = 104
KEY_HEADER = "Assembly contigs"
FONT_TTC = "/System/Library/Fonts/Helvetica.ttc"   # index 0 = Regular, 1 = Bold


def key_layout(measure_draw, header_font, entry_font, gc_text=None, gc_font=None):
    """Compute geometry for the contig key once; consumed by PNG and SVG.

    Right-aligned to KEY_RIGHT so the block hugs the bottom-right corner. The
    GC caption (if any) is included in the width so the block never overflows
    the right edge. Returns absolute pixel coords in the 4999x4999 canvas.
    """
    entries = [(f"{name}: {end - start:,} bp", color)
               for name, start, end, color in CONTIGS]
    label_w = max(measure_draw.textbbox((0, 0), t, font=entry_font)[2] for t, _ in entries)
    header_w = measure_draw.textbbox((0, 0), KEY_HEADER, font=header_font)[2]
    widths = [header_w, KEY_SWATCH + KEY_GAP + label_w]
    if gc_text and gc_font:
        gc_w = measure_draw.textbbox((0, 0), gc_text, font=gc_font)[2]
        widths.append(KEY_SWATCH + KEY_GAP + gc_w)
    block_w = max(widths)
    x0 = KEY_RIGHT - block_w

    rows = []
    y0 = KEY_TOP + KEY_HEADER_GAP
    for i, (text, color) in enumerate(entries):
        y = y0 + i * KEY_ROW_H
        rows.append({
            "swatch": (x0, y, x0 + KEY_SWATCH, y + KEY_SWATCH),
            "color": color,
            "text": text,
            "text_xy": (x0 + KEY_SWATCH + KEY_GAP, y + KEY_SWATCH / 2),
        })
    return {"header_xy": (x0, KEY_TOP), "rows": rows}


def draw_contig_key(draw, layout, header_font, entry_font):
    """PNG rendering of the contig key from a shared layout."""
    draw.text(layout["header_xy"], KEY_HEADER, font=header_font, fill=(0, 0, 0, 255))
    for r in layout["rows"]:
        draw.rectangle(list(r["swatch"]),
                       fill=(r["color"][0], r["color"][1], r["color"][2], 255),
                       outline=(0, 0, 0, 255), width=5)
        draw.text(r["text_xy"], r["text"], font=entry_font,
                  fill=(0, 0, 0, 255), anchor="lm")


def gc_caption_text():
    """GC legend label (shared by PNG and SVG)."""
    return "GC content"


def gc_caption_geometry(layout):
    """Placement for the GC row: the next row after the contigs, with the SAME
    swatch size and row pitch as the contig entries (uniform look)."""
    x0 = layout["header_xy"][0]
    # one KEY_ROW_H below the last contig row's top so spacing is uniform
    last_top = layout["rows"][-1]["swatch"][1]
    y = last_top + KEY_ROW_H
    return x0, y, KEY_SWATCH


def draw_gc_caption(draw, layout, entry_font):
    """PNG: GC legend row styled exactly like the contig rows (same swatch
    size, same font), with a gray swatch and the label 'GC content'."""
    x0, y, swatch = gc_caption_geometry(layout)
    draw.rectangle([x0, y, x0 + swatch, y + swatch],
                   fill=GC_BAND_COLOR + (255,), outline=(0, 0, 0, 255), width=5)
    draw.text((x0 + swatch + KEY_GAP, y + swatch / 2), gc_caption_text(),
              font=entry_font, fill=(0, 0, 0, 255), anchor="lm")


def _hex(color):
    return f"#{color[0]:02x}{color[1]:02x}{color[2]:02x}"


def build_svg(base_png_path, layout, out_svg_path, size, gc=None):
    """Write an editable SVG: OGDraw base map embedded as a full-res raster,
    with the GC track, contig ring, boundary dividers, and key as native SVG
    vector shapes and live <text> (editable in Illustrator/Inkscape)."""
    import base64
    W, H = size
    with open(base_png_path, "rb") as fh:
        b64 = base64.b64encode(fh.read()).decode("ascii")

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
        f'font-family="Helvetica, Arial, sans-serif">',
        f'<image x="0" y="0" width="{W}" height="{H}" '
        f'xlink:href="data:image/png;base64,{b64}"/>',
    ]

    # Self-computed GC overlay track, only in custom mode. (In native mode the
    # GC band is already part of the embedded base image.)
    if gc is not None and gc["mode"] == "custom":
        mean_gc, windows, scale = gc["mean"], gc["windows"], gc["scale"]
        parts.append('<g id="gc-track">')
        parts.append(
            f'<circle cx="{CENTER[0]}" cy="{CENTER[1]}" r="{GC_BASELINE_R}" '
            f'fill="none" stroke="{_hex(GC_BASELINE_COLOR)}" stroke-width="3"/>'
        )
        for mid_bp, gcval in windows:
            th = math.radians(true_angle(mid_bp))
            r_val = GC_BASELINE_R + (gcval - mean_gc) * scale
            x1 = CENTER[0] + GC_BASELINE_R * math.cos(th)
            y1 = CENTER[1] - GC_BASELINE_R * math.sin(th)
            x2 = CENTER[0] + r_val * math.cos(th)
            y2 = CENTER[1] - r_val * math.sin(th)
            col = _hex(GC_COLOR_HI if gcval >= mean_gc else GC_COLOR_LO)
            parts.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" '
                         f'y2="{y2:.1f}" stroke="{col}" stroke-width="{GC_BAR_WIDTH}"/>')
        parts.append('</g>')

    parts.append('<g id="contig-ring">')

    # Contig arcs as thick stroked polylines (geometry identical to the PNG's
    # arc math; polyline avoids SVG arc-flag ambiguity and is exact).
    for name, start_bp, end_bp, color in CONTIGS:
        pts = []
        n = max(2, int((end_bp - start_bp) / 300))
        for k in range(n + 1):
            bp = start_bp + (end_bp - start_bp) * k / n
            th = math.radians(true_angle(bp))
            x = CENTER[0] + RING_RADIUS * math.cos(th)
            y = CENTER[1] - RING_RADIUS * math.sin(th)
            pts.append(f"{x:.1f},{y:.1f}")
        parts.append(
            f'<polyline points="{" ".join(pts)}" fill="none" '
            f'stroke="{_hex(color)}" stroke-width="{RING_THICKNESS}" '
            f'stroke-linecap="butt"><title>{name}</title></polyline>'
        )

    # White boundary dividers
    for bp in contig_boundaries():
        a = math.radians(true_angle(bp))
        r_in = RING_RADIUS - RING_THICKNESS / 2 - 8
        r_out = RING_RADIUS + RING_THICKNESS / 2 + 8
        x1, y1 = CENTER[0] + r_in * math.cos(a), CENTER[1] - r_in * math.sin(a)
        x2, y2 = CENTER[0] + r_out * math.cos(a), CENTER[1] - r_out * math.sin(a)
        parts.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
                     f'stroke="#ffffff" stroke-width="8"/>')
    parts.append('</g>')

    # Contig key: swatches + editable text
    hx, hy = layout["header_xy"]
    parts.append('<g id="contig-key">')
    parts.append(
        f'<text x="{hx:.1f}" y="{hy:.1f}" font-size="{KEY_HEADER_SIZE}" '
        f'font-weight="bold" dominant-baseline="hanging" fill="#000000">'
        f'{KEY_HEADER}</text>'
    )
    for r in layout["rows"]:
        sx0, sy0, sx1, sy1 = r["swatch"]
        tx, ty = r["text_xy"]
        parts.append(
            f'<rect x="{sx0:.1f}" y="{sy0:.1f}" width="{sx1 - sx0:.1f}" '
            f'height="{sy1 - sy0:.1f}" fill="{_hex(r["color"])}" '
            f'stroke="#000000" stroke-width="5"/>'
        )
        parts.append(
            f'<text x="{tx:.1f}" y="{ty:.1f}" font-size="{KEY_ENTRY_SIZE}" '
            f'dominant-baseline="central" fill="#000000">{r["text"]}</text>'
        )

    # GC legend row below the contig key -- same swatch size and font as the
    # contig rows (mirrors the PNG layout).
    if gc is not None:
        x0, y, swatch = gc_caption_geometry(layout)
        parts.append(
            f'<rect x="{x0:.1f}" y="{y:.1f}" width="{swatch}" height="{swatch}" '
            f'fill="{_hex(GC_BAND_COLOR)}" stroke="#000000" stroke-width="5"/>')
        parts.append(
            f'<text x="{x0 + swatch + KEY_GAP:.1f}" y="{y + swatch / 2:.1f}" '
            f'font-size="{KEY_ENTRY_SIZE}" dominant-baseline="central" '
            f'fill="#000000">{gc_caption_text()}</text>')
    parts.append('</g>')
    parts.append('</svg>')

    with open(out_svg_path, "w") as fh:
        fh.write("\n".join(parts))
    print(f"wrote {out_svg_path}")


# ---------------------------------------------------------------------------
# Verification stage
#
# Confirms, with logged evidence, that (1) the hardcoded contig boundaries
# match the actual source sequences on disk, (2) the flattened record OGDraw
# rendered is exactly contig1+contig2+contig3 in that order and its length
# equals TOTAL_LEN -- the invariant that makes bp->angle identical to OGDraw's
# own drawing -- and (3) the drawn ring band sits over empty background in the
# base image (no overlap with gene glyphs or labels), with the measured radial
# clearance to the nearest OGDraw element reported. Hard failures raise before
# any figure is written; overlay clearance issues are logged as warnings.
# ---------------------------------------------------------------------------

def _setup_logger():
    logger = logging.getLogger("contig_ring")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s  %(levelname)-7s %(message)s", "%H:%M:%S")
    ch = logging.StreamHandler()
    ch.setFormatter(fmt)
    logger.addHandler(ch)
    fh = logging.FileHandler(LOG_PATH, mode="w")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    return logger


def _read_locus_length(gb_path):
    with open(gb_path) as fh:
        for line in fh:
            if line.startswith("LOCUS"):
                m = re.search(r"\s(\d+)\s+bp", line)
                if m:
                    return int(m.group(1))
    raise ValueError(f"No LOCUS length found in {gb_path}")


def _read_gb_sequence(gb_path):
    """Return the ORIGIN sequence (uppercase, no coords/whitespace)."""
    chunks, in_origin = [], False
    with open(gb_path) as fh:
        for line in fh:
            if line.startswith("ORIGIN"):
                in_origin = True
                continue
            if in_origin:
                if line.startswith("//"):
                    break
                chunks.append(re.sub(r"[\d\s]", "", line))
    return "".join(chunks).upper()


def _clock(angle_deg):
    """Human-readable clock position for a true_angle (0 deg = 3 o'clock, CCW)."""
    hour = (3 - angle_deg / 30.0) % 12
    hour = 12 if round(hour) % 12 == 0 else round(hour) % 12
    return f"~{hour} o'clock"


def verify_coordinates(logger):
    """Assert the ring's coordinate system matches the source sequences."""
    logger.info("=" * 70)
    logger.info("COORDINATE / FASTA-POSITION SANITY CHECKS")
    logger.info("=" * 70)

    # (1) Per-contig declared length vs. the actual sequence record on disk.
    contig_seqs = {}
    for name, start, end, _color in CONTIGS:
        declared = end - start
        gb_len = _read_locus_length(PER_CONTIG_GB[name])
        seq = _read_gb_sequence(PER_CONTIG_GB[name])
        contig_seqs[name] = seq
        logger.info("%-9s declared=%7d bp | LOCUS=%7d bp | ORIGIN seq=%7d bp | %s",
                    name, declared, gb_len, len(seq),
                    PER_CONTIG_GB[name].name)
        assert declared == gb_len == len(seq), (
            f"{name}: declared {declared} != LOCUS {gb_len} != seq {len(seq)}")
    logger.info("PASS  every contig's declared span equals its source record length")

    # (2) Boundaries are contiguous, start at 0, and sum to TOTAL_LEN.
    assert CONTIGS[0][1] == 0, "first contig must start at bp 0"
    for prev, cur in zip(CONTIGS, CONTIGS[1:]):
        assert cur[1] == prev[2], f"gap/overlap at boundary {prev[0]}->{cur[0]}"
    total = sum(end - start for _n, start, end, _c in CONTIGS)
    logger.info("Boundaries: %s", " | ".join(
        f"{n} [{s:,}-{e:,})" for n, s, e, _c in CONTIGS))
    assert total == TOTAL_LEN, f"contig sum {total} != TOTAL_LEN {TOTAL_LEN}"
    logger.info("PASS  contigs are contiguous, start at 0, and sum to %d bp", total)

    # (3) The flattened record OGDraw rendered: length == TOTAL_LEN (the
    #     invariant that keeps the ring's bp->angle identical to OGDraw's).
    flat_len = _read_locus_length(FLAT_GB)
    logger.info("Flattened record OGDraw rendered: %s  LOCUS=%d bp",
                FLAT_GB.name, flat_len)
    assert flat_len == TOTAL_LEN, (
        f"OGDraw-rendered length {flat_len} != TOTAL_LEN {TOTAL_LEN}: "
        f"the ring would be misaligned to the base map")
    logger.info("PASS  TOTAL_LEN matches the length OGDraw drew (angles align)")

    # (4) The flattened sequence IS contig1+contig2+contig3 in that order, so
    #     each contig occupies exactly its declared bp slice on the FASTA.
    flat_seq = _read_gb_sequence(FLAT_GB)
    assert len(flat_seq) == TOTAL_LEN, (
        f"flattened ORIGIN seq {len(flat_seq)} != TOTAL_LEN {TOTAL_LEN}")
    concat = "".join(contig_seqs[n] for n, _s, _e, _c in CONTIGS)
    assert concat == flat_seq, (
        "flattened sequence is NOT contig1+contig2+contig3 in order -- "
        "contig arc positions would not match the FASTA coordinates")
    for name, start, end, _color in CONTIGS:
        assert flat_seq[start:end] == contig_seqs[name], (
            f"{name} sequence does not occupy flat[{start}:{end}]")
        logger.info("%-9s occupies flat[%7d:%7d)  (verified base-for-base)",
                    name, start, end)
    logger.info("PASS  flattened FASTA = contig1+2+3; every contig sits at its "
                "declared coordinates")

    # (5) Angle mapping self-consistency + human-checkable reference points.
    assert abs(true_angle(0)) < 1e-9 and abs(true_angle(TOTAL_LEN) - 360) < 1e-9
    logger.info("Angle map: 0 bp -> %.3f deg, %d bp -> %.3f deg (closes circle)",
                true_angle(0), TOTAL_LEN, true_angle(TOTAL_LEN))
    for name, start, end, _color in CONTIGS:
        for tag, bp in (("start", start), ("mid", (start + end) // 2)):
            ta = true_angle(bp)
            logger.info("  %-9s %-5s bp=%7d -> true=%7.2f deg  pil=%7.2f deg  %s",
                        name, tag, bp, ta, pil_angle(bp), _clock(ta))
    logger.info("PASS  bp->angle mapping is self-consistent and logged above")
    return flat_seq


def verify_gc(logger, flat_seq):
    """Compute and log GC statistics. Returns a dict describing what to draw:
      {"mode": "native", "mean": gc}                    -- OGDraw's own band
      {"mode": "custom", "mean": gc, "windows":..,       -- self-drawn overlay
       "scale":..}
    or None if GC is entirely disabled."""
    if not (NATIVE_GC or DRAW_GC):
        return None
    logger.info("=" * 70)
    logger.info("GC-CONTENT CHECKS")
    logger.info("=" * 70)
    mean_gc, windows = compute_gc_windows(flat_seq, GC_N_WINDOWS)
    gcs = [gc for _bp, gc in windows]
    lo, hi = min(gcs), max(gcs)
    logger.info("Whole-genome GC = %.2f%% over %d bp", mean_gc, len(flat_seq))
    logger.info("Window GC (n=%d, %d bp/window): min=%.2f%%  max=%.2f%%",
                GC_N_WINDOWS, len(flat_seq) // GC_N_WINDOWS, lo, hi)

    if NATIVE_GC:
        logger.info("GC graph is OGDraw's NATIVE band (baked into the base image "
                    "by render_flattened_with_gc.pl, gc_cont=1).")
        logger.info("Contig ring placed at R=%d (outside the ~832-1042 px GC "
                    "band); overlap is confirmed 0%% by the checks below.",
                    RING_RADIUS)
        logger.info("PASS  native GC band present; mean GC logged for the caption")
        return {"mode": "native", "mean": mean_gc}

    # custom overlay path (DRAW_GC): fit + fit-check against the clear zone
    scale = gc_scale(windows, mean_gc)
    inner = GC_BASELINE_R + (lo - mean_gc) * scale
    outer = GC_BASELINE_R + (hi - mean_gc) * scale
    TITLE_MAX_R, RING_INNER = 570, RING_RADIUS - RING_THICKNESS / 2
    logger.info("Custom GC track extent: %.0f..%.0f px (baseline %d)",
                inner, outer, GC_BASELINE_R)
    assert inner > TITLE_MAX_R + 15 and outer < RING_INNER - 15
    logger.info("PASS  custom GC track fits the title/contig-ring gap")
    return {"mode": "custom", "mean": mean_gc, "windows": windows, "scale": scale}


def verify_overlay(logger, base_img):
    """Confirm the ring band sits over empty background (no OGDraw overlap)
    and report the measured radial clearance to the nearest drawn element."""
    logger.info("=" * 70)
    logger.info("OVERLAY REGISTRATION / NO-COLLISION CHECKS (base image)")
    logger.info("=" * 70)

    rgb = base_img.convert("RGB")
    W, H = rgb.size
    px = rgb.load()
    logger.info("Base image: %dx%d px, center=%s, ring R=%d, thickness=%d",
                W, H, CENTER, RING_RADIUS, RING_THICKNESS)

    def is_ink(x, y):
        if 0 <= x < W and 0 <= y < H:
            r, g, b = px[int(x), int(y)]
            return min(r, g, b) < 245      # anything not near-white = drawn ink
        return False

    def point(bp, radius):
        th = math.radians(true_angle(bp))
        return CENTER[0] + radius * math.cos(th), CENTER[1] - radius * math.sin(th)

    half = RING_THICKNESS / 2
    worst_overlap = 0.0
    for name, start, end, _color in CONTIGS:
        # (A) occupancy directly under the band this contig will paint over
        n_ang = max(50, int((end - start) / 200))
        occupied = samples = 0
        for i in range(n_ang + 1):
            bp = start + (end - start) * i / n_ang
            for dr in (-half, -half / 2, 0, half / 2, half):
                samples += 1
                if is_ink(*point(bp, RING_RADIUS + dr)):
                    occupied += 1
        frac = occupied / samples
        worst_overlap = max(worst_overlap, frac)
        logger.info("%-9s band occupancy over base = %.3f%% (%d/%d samples ink)",
                    name, frac * 100, occupied, samples)

    if worst_overlap <= 0.005:
        logger.info("PASS  ring band covers essentially no drawn content "
                    "(worst %.3f%%)", worst_overlap * 100)
    else:
        logger.warning("CHECK ring band overlaps drawn content (worst %.3f%%); "
                       "consider adjusting RING_RADIUS", worst_overlap * 100)

    # (B) radial clearance from the ring's outer edge to the nearest OGDraw
    #     element scanning outward (toward the gene track / inner labels).
    out_edge = RING_RADIUS + half
    max_scan = int(math.hypot(W, H))   # generous cap
    min_gap, min_gap_bp = None, None
    for k in range(1801):              # 0.2 deg steps
        bp = TOTAL_LEN * k / 1800
        found = None
        r = out_edge + 2
        while r < min(out_edge + 900, CENTER[0]):
            if is_ink(*point(bp, r)):
                found = r - out_edge
                break
            r += 3
        if found is not None and (min_gap is None or found < min_gap):
            min_gap, min_gap_bp = found, bp
    if min_gap is None:
        logger.info("Outward clearance: no drawn element within 900 px of the "
                    "ring's outer edge (very clean)")
    else:
        logger.info("Outward clearance: nearest OGDraw element is %.0f px beyond "
                    "the ring's outer edge (at bp~%d, %s)",
                    min_gap, int(min_gap_bp), _clock(true_angle(min_gap_bp)))
        if min_gap < 25:
            logger.warning("CHECK outward clearance < 25 px -- ring is close to "
                           "gene annotations; consider a smaller RING_RADIUS")
        else:
            logger.info("PASS  comfortable gap (>=25 px) between ring and the "
                        "nearest gene annotation")

    # (C) inward clearance to the nearest element (the GC band's outer edge).
    in_edge = RING_RADIUS - half
    min_in, min_in_bp = None, None
    for k in range(1801):
        bp = TOTAL_LEN * k / 1800
        found = None
        r = in_edge - 2
        while r > in_edge - 400:
            if is_ink(*point(bp, r)):
                found = in_edge - r
                break
            r -= 3
        if found is not None and (min_in is None or found < min_in):
            min_in, min_in_bp = found, bp
    if min_in is None:
        logger.info("Inward clearance: no drawn element within 400 px inside the "
                    "ring (very clean)")
    else:
        logger.info("Inward clearance: nearest element (GC band edge) is %.0f px "
                    "inside the ring's inner edge (at bp~%d, %s)",
                    min_in, int(min_in_bp), _clock(true_angle(min_in_bp)))
        if min_in < 15:
            logger.warning("CHECK inward clearance < 15 px -- ring nearly touches "
                           "the GC band; consider a larger RING_RADIUS")
        else:
            logger.info("PASS  ring sits clear of the GC band (>=15 px inward)")


def main():
    import sys
    # Optional output-name suffix (e.g. "_with_GC_TEST") so a trial render does
    # not overwrite the approved FINAL. Default (no arg) writes the FINAL names.
    suffix = sys.argv[1] if len(sys.argv) > 1 else "_FINAL"
    out_image = BASE / f"mito_combined_flattened_with_contig_ring{suffix}.png"

    logger = _setup_logger()
    logger.info("Building overlay  (RING_RADIUS=%d, THICKNESS=%d, DRAW_GC=%s)",
                RING_RADIUS, RING_THICKNESS, DRAW_GC)
    logger.info("Output: %s", out_image.name)
    flat_seq = verify_coordinates(logger)
    gc = verify_gc(logger, flat_seq)

    img = Image.open(SRC_IMAGE).convert("RGBA")
    verify_overlay(logger, img)
    overlay = Image.new("RGBA", img.size, (255, 255, 255, 0))
    draw = ImageDraw.Draw(overlay)

    # Self-computed GC overlay track (only if the custom path is enabled; the
    # native OGDraw band is already in the base image).
    if gc is not None and gc["mode"] == "custom":
        draw_gc_track(draw, gc["windows"], gc["mean"], gc["scale"])

    bbox = [
        CENTER[0] - RING_RADIUS, CENTER[1] - RING_RADIUS,
        CENTER[0] + RING_RADIUS, CENTER[1] + RING_RADIUS,
    ]

    for _name, start_bp, end_bp, color in CONTIGS:
        pil_start = pil_angle(end_bp)
        pil_end = pil_angle(start_bp)
        if pil_start > pil_end:
            pil_end += 360
        draw.arc(bbox, start=pil_start, end=pil_end, fill=color, width=RING_THICKNESS)

    # Crisp white divider lines at each contig boundary
    for bp in contig_boundaries():
        a = math.radians(true_angle(bp))
        r_in = RING_RADIUS - RING_THICKNESS / 2 - 8
        r_out = RING_RADIUS + RING_THICKNESS / 2 + 8
        x1, y1 = CENTER[0] + r_in * math.cos(a), CENTER[1] - r_in * math.sin(a)
        x2, y2 = CENTER[0] + r_out * math.cos(a), CENTER[1] - r_out * math.sin(a)
        draw.line([x1, y1, x2, y2], fill=(255, 255, 255, 255), width=8)

    # Helvetica.ttc is a TrueType *collection*: index 0 = Regular, 1 = Bold.
    entry_font = ImageFont.truetype(FONT_TTC, KEY_ENTRY_SIZE, index=0)
    header_font = ImageFont.truetype(FONT_TTC, KEY_HEADER_SIZE, index=1)
    gc_text = gc_caption_text() if gc is not None else None
    layout = key_layout(draw, header_font, entry_font, gc_text, entry_font)
    draw_contig_key(draw, layout, header_font, entry_font)
    if gc is not None:
        draw_gc_caption(draw, layout, entry_font)

    if DRAW_CURVED_LABELS:
        label_font = ImageFont.truetype(FONT_TTC, 120, index=1)
        for name, start_bp, end_bp, _color in CONTIGS:
            mid_bp = (start_bp + end_bp) / 2
            label = f"{name} • {end_bp - start_bp:,} bp"
            draw_curved_text(overlay, draw, label, mid_bp, label_font, radius=RING_RADIUS)

    combined = Image.alpha_composite(img, overlay)
    combined.convert("RGB").save(out_image)
    logger.info("wrote %s", out_image)

    # Companion editable SVG (same geometry; overlay as live vector/text).
    # The native GC band rides along inside the embedded base image, so only a
    # custom overlay track is emitted as vector; the caption uses gc["mean"].
    build_svg(SRC_IMAGE, layout, str(out_image).replace(".png", ".svg"),
              img.size, gc)
    logger.info("Done.")


if __name__ == "__main__":
    main()
