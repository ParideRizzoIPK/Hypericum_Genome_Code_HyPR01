#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --partition=cpu        # OPTIMIZED: Target mixed cpu nodes for backfilling
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16     # OPTIMIZED: Scaled down from 32 to trigger immediate scheduling
#SBATCH --mem=32G              # OPTIMIZED: Scaled down from 64G to fit empty memory slots
#SBATCH --time=02:00:00        # OPTIMIZED: Lowered to 2 hours to leverage backfill engine
#SBATCH --job-name=NUMT_NUPT
#SBATCH --output=logs/numt_nupt_%j.out
#SBATCH --error=logs/numt_nupt_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  

echo "#########################################################"
echo "Pipeline: Step 1 (Script R3) - NUMT/NUPT Organellar Detection"
echo "Date: $(date)"
echo "#########################################################"

# --- 1. ENVIRONMENT & MODULE INITIALIZATION ---
source ~/.bashrc

# Load native cluster alignment tools
module load minimap2/2.24

set -u            

export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# --- 2. CONFIGURATION AND DIRECTORIES ---
MASTER_OUT_DIR="/path/to/output/directory"
mkdir -p "${MASTER_OUT_DIR}"

# FIXED: Swapped out legacy NCBI references for in-house HyPR01 assemblies
CP_REF="/path/to/input/directory/HyPR01_Chloroplast_Genome.fasta"
MT_REF="/path/to/input/directory/HyPR01_Mitochondria_genome_collapsed.fasta"

HAP1_SOFTMASKED="/path/to/input/directory/hap1.softmasked.fa"
HAP2_SOFTMASKED="/path/to/input/directory/hap2.softmasked.fa"

THREADS=${SLURM_CPUS_PER_TASK}
HALF_THREADS=$(( THREADS / 2 ))

echo "Allocated total CPUs: ${THREADS} cores."
echo "Splitting pipeline resource footprint: ${HALF_THREADS} threads per haplotype track."

# ==============================================================================
# MODULE R3.1: ORGANELLAR DATABASE CONCATENATION & HEADER STANDARDIZATION
# ==============================================================================
echo ">>>> STANDARDIZING AND CONCATENATING ORGANELLAR FASTA QUERIES... <<<<"
COMBINED_REF="${MASTER_OUT_DIR}/combined_organelles.fasta"

# Use python to safely standardize multi-scaffold or single circular organelle headers
python3 -c "
def process_fasta(infile, outfile, prefix):
    idx = 1
    with open(infile, 'r') as f_in, open(outfile, 'a') as f_out:
        for line in f_in:
            if line.startswith('>'):
                f_out.write(f'>{prefix}_{idx}\n')
                idx += 1
            else:
                f_out.write(line)

import os
if os.path.exists('${COMBINED_REF}'): os.remove('${COMBINED_REF}')
process_fasta('${CP_REF}', '${COMBINED_REF}', 'Plastid')
process_fasta('${MT_REF}', '${COMBINED_REF}', 'Mitochondrion')
"
echo "  -> Combined reference library generated at: ${COMBINED_REF}"

# ==============================================================================
# PIPELINE EXECUTION LOOP: PARALLEL HOMOLOGY STREAMING
# ==============================================================================

# --- RUNTRACK 1: HAPLOTYPE 1 TRACK ---
echo ">>>> DETACHING HAPLOTYPE 1 MAPPING LOOP... <<<<"
HAP1_WSPACE="${MASTER_OUT_DIR}/Hap1_wspace"
mkdir -p "${HAP1_WSPACE}"

(
    cd "${HAP1_WSPACE}"
    echo "[$(date)] Aligning organelle genomes to Hap1 nuclear space via minimap2..." > r3_hap1.log
    
    # Run minimap2 using '-x asm20' to ensure deep structural sensitivity for ancient inserts
    minimap2 -x asm20 -c -t "${HALF_THREADS}" "${HAP1_SOFTMASKED}" "${COMBINED_REF}" > hap1_organelle_hits.paf 2>> r3_hap1.log
    echo "[$(date)] Alignment mapping completed. Engaging GFF3 Parser..." >> r3_hap1.log
) &
HAP1_PID=$!


# --- RUNTRACK 2: HAPLOTYPE 2 TRACK ---
echo ">>>> DETACHING HAPLOTYPE 2 MAPPING LOOP... <<<<"
HAP2_WSPACE="${MASTER_OUT_DIR}/Hap2_wspace"
mkdir -p "${HAP2_WSPACE}"

(
    cd "${HAP2_WSPACE}"
    echo "[$(date)] Aligning organelle genomes to Hap2 nuclear space via minimap2..." > r3_hap2.log
    
    minimap2 -x asm20 -c -t "${HALF_THREADS}" "${HAP2_SOFTMASKED}" "${COMBINED_REF}" > hap2_organelle_hits.paf 2>> r3_hap2.log
    echo "[$(date)] Alignment mapping completed. Engaging GFF3 Parser..." >> r3_hap2.log
) &
HAP2_PID=$!

# ==============================================================================
# WAITING BARRIER & PARSING ENGINE
# ==============================================================================
echo ">>>> Synchronization barrier engaged: Awaiting structural maps... <<<<"
wait "${HAP1_PID}"
echo "  -> Haplotype 1 tracking thread terminated."
wait "${HAP2_PID}"
echo "  -> Haplotype 2 tracking thread terminated."

echo ">>>> INITIALIZING CUSTOM GFF3 PARSING BLOCKS... <<<<"

# Embedded Python routine handles length filtering, sequence identity calculations, 
# flat structure normalization, and 0-to-1 based coordinate translations natively.
cat << 'EOF' > parse_paf_to_flat_gff3.py
import sys

paf_file = sys.argv[1]
gff_file = sys.argv[2]
hap_tag = sys.argv[3]

records = []

with open(paf_file, 'r') as infile:
    for line in infile:
        if not line.strip():
            continue
        p = line.strip().split('\t')
        
        qname = p[0]       # Standardized Query Header (Plastid_X / Mitochondrion_X)
        strand = p[4]      # Relative strand orientation (+ / -)
        tname = p[5]       # Target chromosome/nuclear scaffold identifier
        tstart = int(p[7]) # 0-based coordinate start position
        tend = int(p[8])   # 0-based coordinate end position
        n_match = int(p[9])# Number of residue matches
        aln_len = int(p[10])# Total alignment block length
        mapq = p[11]       # Alignment Mapping Quality Score

        # Criterion 1: Absolute Length Threshold
        if aln_len < 500:
            continue
            
        # Criterion 2: Absolute Sequence Identity Calculation
        identity_pct = (n_match / aln_len) * 100
        if identity_pct < 90.0:
            continue
            
        # Separate structural feature types based on database prefix mapping
        if "Plastid" in qname:
            gff_type = "NUPT"
            organelle_source = "Plastid"
        elif "Mitochondrion" in qname:
            gff_type = "NUMT"
            organelle_source = "Mitochondrion"
        else:
            gff_type = "Organelle_Insertion"
            organelle_source = "Unknown"
            
        # Criterion 3: Coordinate system translation (0-based half-open -> 1-based closed)
        gff_start = tstart + 1
        gff_end = tend
        
        records.append({
            'seqid': tname,
            'source': 'NUMT_NUPT',
            'type': gff_type,
            'start': gff_start,
            'end': gff_end,
            'score': mapq,
            'strand': strand,
            'phase': '.',
            'organelle': organelle_source,
            'identity': f"{identity_pct:.2f}"
        })

# Perform standard positional sorting by chromosome and then start coordinates
records.sort(key=lambda x: (x['seqid'], x['start']))

with open(gff_file, 'w') as outfile:
    outfile.write("##gff-version 3\n")
    feature_counter = 1
    for r in records:
        # Build clean, flat attribute string lacking redundant parent records
        attributes = f"ID={hap_tag}_{r['type']}_{feature_counter:05d};organelle={r['organelle']};percent_identity={r['identity']}"
        outfile.write(f"{r['seqid']}\t{r['source']}\t{r['type']}\t{r['start']}\t{r['end']}\t{r['score']}\t{r['strand']}\t{r['phase']}\t{attributes}\n")
        feature_counter += 1

print(f"Successfully compiled {feature_counter - 1} records into {gff_file}")
EOF

# Run parsing loops across the haplotype outputs
python3 parse_paf_to_flat_gff3.py "${HAP1_WSPACE}/hap1_organelle_hits.paf" "${MASTER_OUT_DIR}/hap1_NUMT_NUPT.gff3" "HyPR01_Hap1"
python3 parse_paf_to_flat_gff3.py "${HAP2_WSPACE}/hap2_organelle_hits.paf" "${MASTER_OUT_DIR}/hap2_NUMT_NUPT.gff3" "HyPR01_Hap2"

# --- CLEANUP SCRATCH DIRECTORIES ---
rm -f parse_paf_to_flat_gff3.py
rm -rf "${HAP1_WSPACE}" "${HAP2_WSPACE}"

echo "#########################################################"
echo "SCRIPT R3 PIPELINE RUN COMPLETE."
echo "Production Flat GFF3 Tracks safely stored:"
echo "  - Hap1 GFF3 Output: ${MASTER_OUT_DIR}/hap1_NUMT_NUPT.gff3"
echo "  - Hap2 GFF3 Output: ${MASTER_OUT_DIR}/hap2_NUMT_NUPT.gff3"
echo "#########################################################"