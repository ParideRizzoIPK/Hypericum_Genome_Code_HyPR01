#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --partition=cpu        
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8      
#SBATCH --mem=16G              
#SBATCH --time=01:00:00        
#SBATCH --job-name=organelle_synteny
#SBATCH --output=logs/synteny_%j.out
#SBATCH --error=logs/synteny_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  

echo "########################################################################"
echo "Pipeline: Independent Per-Contig Orientation & Curved Synteny Engine"
echo "Date: $(date)"
echo "########################################################################"

# --- 1. ENVIRONMENT INITIALIZATION ---
source ~/.bashrc
module load minimap2/2.24
set -u            

export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# --- 2. CONFIGURATION AND DIRECTORIES ---
MASTER_OUT_DIR="/path/to/output/directory"
mkdir -p "${MASTER_OUT_DIR}"

HYPR01_CP="/path/to/input/directory/HyPR01_Chloroplast_Genome.fasta"
HYPR01_MT="/path/to/input/directory/HyPR01_Mitochondria_genome_collapsed.fasta"

NCBI_CP="/path/to/input/directory/hypericum_cp.fasta"
NCBI_MT="/path/to/input/directory/hypericum_mt.fasta"

THREADS=${SLURM_CPUS_PER_TASK}

cd "${MASTER_OUT_DIR}"

# ==============================================================================
# MODULE 1: GLOBAL PER-CONTIG STRAND VOTING & CURATION ENGINE
# ==============================================================================
echo ">>>> EXECUTING PRELIMINARY STRUCTURAL SCANS FOR STRAND CHECKING... <<<<"
minimap2 -x asm5 -c -t "${THREADS}" "${NCBI_CP}" "${HYPR01_CP}" > pre_cp.paf 2>/dev/null
minimap2 -x asm5 -c -t "${THREADS}" "${NCBI_MT}" "${HYPR01_MT}" > pre_mt.paf 2>/dev/null

cat << 'EOF' > preprocess_contigs_globally.py
import sys
import os
from collections import defaultdict

def read_fasta(path):
    seqs = {}
    current_id = None
    current_seq = []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_id: seqs[current_id] = "".join(current_seq)
                current_id = line.strip().split()[0][1:]
                current_seq = []
            else:
                current_seq.append(line.strip())
    if current_id: seqs[current_id] = "".join(current_seq)
    return seqs

def process_organelle_tracks(paf_path, fasta_path, out_fasta_path, is_cp=False):
    seqs = read_fasta(fasta_path)
    contig_alns = defaultdict(list)
    
    if os.path.exists(paf_path) and os.path.getsize(paf_path) > 0:
        with open(paf_path, 'r') as f:
            for line in f:
                p = line.strip().split('\t')
                if len(p) < 12: continue
                contig_alns[p[0]].append({
                    'qlen': int(p[1]), 'qstart': int(p[2]), 'qend': int(p[3]),
                    'strand': p[4], 'tname': p[5], 'tlen': int(p[6]),
                    'tstart': int(p[7]), 'tend': int(p[8])
                })
                
    processed_seqs = {}
    contig_target_starts = {}
    comp_map = str.maketrans("ATCGNatsgn", "TAGCNtagcn")
    
    for cid, seq in seqs.items():
        alns = contig_alns[cid]
        fwd_bases = sum((a['tend'] - a['tstart']) for a in alns if a['strand'] == '+')
        rev_bases = sum((a['tend'] - a['tstart']) for a in alns if a['strand'] == '-')
        
        # FIX 1: Bulletproof Independent Contig Curation Check
        if rev_bases > fwd_bases:
            print(f"  [Reorientation] Contig '{cid}' is predominantly inverted. Flipping to forward strand.")
            seq = seq.translate(comp_map)[::-1]
            for a in alns:
                a['strand'] = '+' if a['strand'] == '-' else '-'
                old_qs = a['qstart']
                a['qstart'] = a['qlen'] - a['qend']
                a['qend'] = a['qlen'] - old_qs
                
        min_ts = min([a['tstart'] for a in alns], default=float('inf'))
        contig_target_starts[cid] = min_ts
        
        # Circular Phasing Optimization (Applied to Chloroplast)
        if is_cp and alns:
            best_aln = min(alns, key=lambda x: x['tstart'])
            rotation_point = best_aln['qstart'] if best_aln['strand'] == '+' else best_aln['qend']
            print(f"  [Phasing] Rotating circular molecule at coordinate breakpoint: {rotation_point} bp")
            seq = seq[rotation_point:] + seq[:rotation_point]
            
        processed_seqs[cid] = seq
        
    # Syntenic Sort Order Array
    sorted_cids = sorted(processed_seqs.keys(), key=lambda k: contig_target_starts.get(k, float('inf')))
    
    with open(out_fasta_path, 'w') as out_f:
        for cid in sorted_cids:
            out_f.write(f">{cid}\n")
            s = processed_seqs[cid]
            for i in range(0, len(s), 80):
                out_f.write(s[i:i+80] + "\n")

print("Executing Chloroplast curation...")
process_organelle_tracks('pre_cp.paf', sys.argv[1], 'HyPR01_Chloroplast_Genome_curated.fasta', is_cp=True)
print("Executing Mitochondrial multi-contig curation...")
process_organelle_tracks('pre_mt.paf', sys.argv[2], 'HyPR01_Mitochondria_genome_curated.fasta', is_cp=False)
EOF

python3 preprocess_contigs_globally.py "${HYPR01_CP}" "${HYPR01_MT}"
rm -f pre_cp.paf pre_mt.paf preprocess_contigs_globally.py

# ==============================================================================
# MODULE 2: PRODUCTION STRUCTURAL ALIGNMENTS
# ==============================================================================
echo ">>>> COMPUTING PRODUCTION MOLECULAR MAPS... <<<<"
minimap2 -x asm5 -c -t "${THREADS}" "${NCBI_CP}" "HyPR01_Chloroplast_Genome_curated.fasta" > final_cp.paf 2> cp_final.log
minimap2 -x asm5 -c -t "${THREADS}" "${NCBI_MT}" "HyPR01_Mitochondria_genome_curated.fasta" > final_mt.paf 2> mt_final.log

# ==============================================================================
# MODULE 3: RE-SCALED PUBLICATION GRAPHICS RENDERER (SVG ENGINE)
# ==============================================================================
echo ">>>> RENDER SEPARATE VECTOR IMAGES WITH RE-SCALED FONTS... <<<<"

cat << 'EOF' > render_publication_plots.py
import sys
import os

def read_fasta_lens(path):
    lens = {}
    current_id = None
    climb = 0
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_id: lens[current_id] = climb
                current_id = line.strip().split()[0][1:]
                climb = 0
            else:
                climb += len(line.strip())
    if current_id: lens[current_id] = climb
    return lens

def generate_svg_plot(paf_path, fasta_path, output_svg, title_text):
    query_lens = read_fasta_lens(fasta_path)
    query_offsets = {}
    contig_boundaries = []
    current_offset = 0
    gap_size_bp = 4000  
    
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                cid = line.strip().split()[0][1:]
                length = query_lens[cid]
                query_offsets[cid] = current_offset
                contig_boundaries.append((current_offset, current_offset + length, cid))
                current_offset += length + gap_size_bp
                
    total_q_span = current_offset - gap_size_bp
    
    alignments = []
    tlen = 0
    with open(paf_path, 'r') as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) < 12: continue
            alignments.append({
                'qname': p[0], 'qstart': int(p[2]), 'qend': int(p[3]), 'strand': p[4],
                'tlen': int(p[6]), 'tstart': int(p[7]), 'tend': int(p[8])
            })
            tlen = max(tlen, int(p[6]))

    # FIX 3: Expanded canvas properties to guarantee margin clearance
    margin_left, margin_right = 260, 80
    width, height = 1500, 420
    plot_width = width - margin_left - margin_right
    
    max_total_len = max(tlen, total_q_span)
    scale = plot_width / max_total_len
    
    y_ncbi, y_hypr = 120, 270
    track_h = 18
    y_mid = (y_ncbi + y_hypr) / 2
    
    with open(output_svg, 'w') as f:
        f.write(f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg" style="background-color:white;">\n')
        
        # FIX 5: Integrated high-contrast font sizes matching reference scales
        f.write('<style>\n')
        f.write('  .title { font-family: sans-serif; font-size: 22px; font-weight: bold; fill: #111; }\n')
        f.write('  .label { font-family: sans-serif; font-size: 14px; font-weight: bold; fill: #222; }\n')
        f.write('  .subtext { font-family: sans-serif; font-size: 12px; fill: #555; }\n')
        f.write('  .tick-text { font-family: sans-serif; font-size: 12px; font-weight: bold; fill: #333; }\n')
        f.write('  .legend-text { font-family: sans-serif; font-size: 12px; font-weight: bold; fill: #222; }\n')
        f.write('</style>\n')
        
        f.write(f'  <text x="{width/2}" y="50" text-anchor="middle" class="title">{title_text}</text>\n')
        
        # Render Curved Ribbon Connectors
        for aln in alignments:
            ts = margin_left + aln['tstart'] * scale
            te = margin_left + aln['tend'] * scale
            
            offset = query_offsets[aln['qname']]
            qs = margin_left + (offset + aln['qstart']) * scale
            qe = margin_left + (offset + aln['qend']) * scale
            
            if aln['strand'] == '+':
                color = '#008080'
                path_d = f"M {ts},{y_ncbi+track_h/2} L {te},{y_ncbi+track_h/2} C {te},{y_mid} {qe},{y_mid} {qe},{y_hypr-track_h/2} L {qs},{y_hypr-track_h/2} C {qs},{y_mid} {ts},{y_mid} {ts},{y_ncbi+track_h/2} Z"
            else:
                color = '#FF4500'
                path_d = f"M {ts},{y_ncbi+track_h/2} L {te},{y_ncbi+track_h/2} C {te},{y_mid} {qs},{y_mid} {qs},{y_hypr-track_h/2} L {qe},{y_hypr-track_h/2} C {qe},{y_mid} {ts},{y_mid} {ts},{y_ncbi+track_h/2} Z"
            
            f.write(f'  <path d="{path_d}" fill="{color}" fill-opacity="0.45" stroke="{color}" stroke-opacity="0.6" stroke-width="0.5"/>\n')
            
        # Draw Master NCBI Reference Bar
        f.write(f'  <rect x="{margin_left}" y="{y_ncbi-track_h/2}" width="{tlen*scale}" height="{track_h}" rx="4" fill="#EAEAEA" stroke="#666" stroke-width="1.5"/>\n')
        
        # Draw Segmented HyPR01 Contig Blocks
        for start, end, name in contig_boundaries:
            b_width = (end - start) * scale
            b_x = margin_left + start * scale
            f.write(f'  <rect x="{b_x}" y="{y_hypr-track_h/2}" width="{b_width}" height="{track_h}" rx="2" fill="#DCDCDA" stroke="#555" stroke-width="1.5"/>\n')
            if len(contig_boundaries) > 1:
                f.write(f'  <line x1="{b_x}" y1="{y_hypr-track_h/2-4}" x2="{b_x}" y2="{y_hypr+track_h/2+4}" stroke="#CC0000" stroke-width="1" stroke-dasharray="2,2"/>\n')

        # FIX 4: Simplified track labels to eliminate redundancy
        f.write(f'  <text x="{margin_left-25}" y="{y_ncbi+2}" text-anchor="end" class="label">NCBI Reference</text>\n')
        f.write(f'  <text x="{margin_left-25}" y="{y_ncbi+18}" text-anchor="end" class="subtext">({tlen/1000.0:.1f} kb)</text>\n')
        
        f.write(f'  <text x="{margin_left-25}" y="{y_hypr+2}" text-anchor="end" class="label">HyPR01</text>\n')
        f.write(f'  <text x="{margin_left-25}" y="{y_hypr+18}" text-anchor="end" class="subtext">({sum(query_lens.values())/1000.0:.1f} kb)</text>\n')
        
        # Dynamic Axis Scale Ticks
        tick_interval = 20000 if max_total_len < 200000 else 50000
        for bp in range(0, int(max_total_len), tick_interval):
            x_pos = margin_left + bp * scale
            if x_pos <= width - margin_right:
                f.write(f'  <line x1="{x_pos}" y1="{y_hypr+track_h/2}" x2="{x_pos}" y2="{y_hypr+track_h/2+6}" stroke="#333" stroke-width="1.5"/>\n')
                f.write(f'  <text x="{x_pos}" y="{y_hypr+track_h/2+20}" text-anchor="middle" class="tick-text">{int(bp/1000)} kb</text>\n')
        
        # Legend Configuration
        f.write(f'  <rect x="{width-260}" y="{height-50}" width="18" height="12" fill="#008080" fill-opacity="0.6"/>\n')
        f.write(f'  <text x="{width-230}" y="{height-40}" class="legend-text">Forward Collinear Match</text>\n')
        f.write(f'  <rect x="{width-260}" y="{height-30}" width="18" height="12" fill="#FF4500" fill-opacity="0.6"/>\n')
        f.write(f'  <text x="{width-230}" y="{height-20}" class="legend-text">Inverted Match</text>\n')
        f.write('</svg>\n')
    print(f"  -> Generated publication vector: {output_svg}")

generate_svg_plot("final_cp.paf", "HyPR01_Chloroplast_Genome_curated.fasta", "Chloroplast_Synteny_NCBI_vs_HyPR01.svg", "Chloroplast Pairwise Genome Synteny Architecture")
generate_svg_plot("final_mt.paf", "HyPR01_Mitochondria_genome_curated.fasta", "Mitochondrion_Synteny_NCBI_vs_HyPR01.svg", "Mitochondrion Pairwise Genome Synteny Architecture")
EOF

python3 render_publication_plots.py
rm -f render_publication_plots.py

echo "########################################################################"
echo "PIPELINE COMPLETED SUCCESSFULLY."
echo "Output files stored at: ${MASTER_OUT_DIR}/"
echo "########################################################################"