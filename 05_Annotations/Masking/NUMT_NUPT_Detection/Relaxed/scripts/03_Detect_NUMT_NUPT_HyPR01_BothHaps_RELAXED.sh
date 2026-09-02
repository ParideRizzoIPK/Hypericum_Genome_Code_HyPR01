#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --partition=cpu        
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16     
#SBATCH --mem=32G              
#SBATCH --time=02:00:00        
#SBATCH --job-name=NUMT_NUPT_relaxed
#SBATCH --output=logs/numt_nupt_relaxed_%j.out
#SBATCH --error=logs/numt_nupt_relaxed_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  

echo "########################################################################"
echo "Pipeline: Step 1 (Script R3b) - Refined Sensitive NUMT/NUPT Mapping"
echo "Date: $(date)"
echo "########################################################################"

# --- 1. ENVIRONMENT & MODULE INITIALIZATION ---
source ~/.bashrc
module load minimap2/2.24
set -u            

export APPTAINER_BIND="/path/to/your/directory"
export SINGULARITY_BIND="/path/to/your/directory"

# --- 2. CONFIGURATION AND DIRECTORIES ---
MASTER_OUT_DIR="/path/to/output/directory"
mkdir -p "${MASTER_OUT_DIR}"

# FIXED: Replaced legacy NCBI dependencies with clean in-house HyPR01 genome targets
CP_REF="/path/to/input/directory/HyPR01_Chloroplast_Genome.fasta"
MT_REF="/path/to/input/directory/HyPR01_Mitochondria_genome_collapsed.fasta"

HAP1_SOFTMASKED="/path/to/input/directory/hap1.softmasked.fa"
HAP2_SOFTMASKED="/path/to/input/directory/hap2.softmasked.fa"

THREADS=${SLURM_CPUS_PER_TASK}
HALF_THREADS=$(( THREADS / 2 ))

# ==============================================================================
# MODULE R3b.1: ORGANELLAR DATABASE CONCATENATION & HEADER STANDARDIZATION
# ==============================================================================
COMBINED_REF="${MASTER_OUT_DIR}/combined_organelles.fasta"

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

# ==============================================================================
# PIPELINE EXECUTION LOOP: PARALLEL INVERTED STREAMING
# ==============================================================================

# --- RUNTRACK 1: HAPLOTYPE 1 TRACK ---
HAP1_WSPACE="${MASTER_OUT_DIR}/Hap1_relaxed_wspace"
mkdir -p "${HAP1_WSPACE}"

(
    cd "${HAP1_WSPACE}"
    # REFINED: Switched to -k 15 -w 4 to suppress micro-homology noise while maintaining divergence sensitivity
    minimap2 -k 15 -w 4 -s 40 -c -t "${HALF_THREADS}" "${COMBINED_REF}" "${HAP1_SOFTMASKED}" > hap1_inverted_hits.paf 2> r3b_hap1.log
) &
HAP1_PID=$!

# --- RUNTRACK 2: HAPLOTYPE 2 TRACK ---
HAP2_WSPACE="${MASTER_OUT_DIR}/Hap2_relaxed_wspace"
mkdir -p "${HAP2_WSPACE}"

(
    cd "${HAP2_WSPACE}"
    minimap2 -k 15 -w 4 -s 40 -c -t "${HALF_THREADS}" "${COMBINED_REF}" "${HAP2_SOFTMASKED}" > hap2_inverted_hits.paf 2> r3b_hap2.log
) &
HAP2_PID=$!

# ==============================================================================
# WAITING BARRIER & FIXED GFF3 PARSING ENGINE
# ==============================================================================
wait "${HAP1_PID}"
wait "${HAP2_PID}"

echo ">>>> INITIALIZING REFINED COALESCENCE PARSING BLOCKS... <<<<"

cat << 'EOF' > parse_inverted_paf_to_coalesced_gff3.py
import sys
from collections import defaultdict

paf_file = sys.argv[1]
gff_file = sys.argv[2]
hap_tag = sys.argv[3]

raw_hits = []

with open(paf_file, 'r') as infile:
    for line in infile:
        if not line.strip():
            continue
        p = line.strip().split('\t')
        if len(p) < 12:
            continue
        
        tname = p[0]        # Nuclear scaffold ID (Query Name)
        tstart = int(p[2])  # Nuclear start position (Query Start)
        tend = int(p[3])    # Nuclear end position (Query End)
        strand = p[4]       
        qname = p[5]        # FIXED: Organelle ID string is Column 6 (Target Name)
        n_match = int(p[9]) 
        aln_len = int(p[10])
        mapq = p[11]

        # REFINED: Enforce standard 500 bp limit to filter out low-complexity SSRs
        if aln_len < 500:
            continue
        identity_pct = (n_match / aln_len) * 100
        if identity_pct < 75.0:
            continue
            
        if "Plastid" in qname:
            gff_type = "NUPT"
            organelle_source = "Plastid"
        elif "Mitochondrion" in qname:
            gff_type = "NUMT"
            organelle_source = "Mitochondrion"
        else:
            gff_type = "Organelle_Insertion"
            organelle_source = "Unknown"
            
        gff_start = tstart + 1
        gff_end = tend
        
        raw_hits.append({
            'seqid': tname,
            'strand': strand,
            'type': gff_type,
            'organelle': organelle_source,
            'start': gff_start,
            'end': gff_end,
            'n_match': n_match,
            'aln_len': aln_len,
            'score': int(mapq)
        })

groups = defaultdict(list)
for hit in raw_hits:
    key = (hit['seqid'], hit['strand'], hit['type'], hit['organelle'])
    groups[key].append(hit)

merged_records = []

for key, hits in groups.items():
    seqid, strand, gff_type, organelle = key
    hits.sort(key=lambda x: x['start'])
    
    current_merged = None
    for h in hits:
        if current_merged is None:
            current_merged = {
                'seqid': seqid, 'strand': strand, 'type': gff_type, 'organelle': organelle,
                'start': h['start'], 'end': h['end'],
                'total_matches': h['n_match'], 'total_aln_len': h['aln_len'],
                'max_score': h['score'], 'gaps': []
            }
        else:
            if h['start'] <= current_merged['end'] + 100:
                if h['start'] > current_merged['end'] + 1:
                    gap_start = current_merged['end'] + 1
                    gap_end = h['start'] - 1
                    current_merged['gaps'].append(f"{gap_start}-{gap_end}")
                
                current_merged['end'] = max(current_merged['end'], h['end'])
                current_merged['total_matches'] += h['n_match']
                current_merged['total_aln_len'] += h['aln_len']
                current_merged['max_score'] = max(current_merged['max_score'], h['score'])
            else:
                merged_records.append(current_merged)
                current_merged = {
                    'seqid': seqid, 'strand': strand, 'type': gff_type, 'organelle': organelle,
                    'start': h['start'], 'end': h['end'],
                    'total_matches': h['n_match'], 'total_aln_len': h['aln_len'],
                    'max_score': h['score'], 'gaps': []
                }
    if current_merged is not None:
        merged_records.append(current_merged)

merged_records.sort(key=lambda x: (x['seqid'], x['start']))

with open(gff_file, 'w') as outfile:
    outfile.write("##gff-version 3\n")
    outfile.write("#seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes\n")
    
    feature_counter = 1
    for r in merged_records:
        weighted_identity = (r['total_matches'] / r['total_aln_len']) * 100
        gap_str = ",".join(r['gaps']) if r['gaps'] else "none"
        
        attributes = f"ID={hap_tag}_{r['type']}_{feature_counter:05d};organelle={r['organelle']};percent_identity={weighted_identity:.2f};gap_intervals={gap_str}"
        outfile.write(f"{r['seqid']}\tNUMT_NUPT_relaxed\t{r['type']}\t{r['start']}\t{r['end']}\t{r['max_score']}\t{r['strand']}\t.\t{attributes}\n")
        feature_counter += 1

print(f"Successfully compiled {feature_counter - 1} records into {gff_file}")
EOF

python3 parse_inverted_paf_to_coalesced_gff3.py "${HAP1_WSPACE}/hap1_inverted_hits.paf" "${MASTER_OUT_DIR}/hap1_NUMT_NUPT.gff3" "HyPR01_Hap1"
python3 parse_inverted_paf_to_coalesced_gff3.py "${HAP2_WSPACE}/hap2_inverted_hits.paf" "${MASTER_OUT_DIR}/hap2_NUMT_NUPT.gff3" "HyPR01_Hap2"

rm -f parse_inverted_paf_to_coalesced_gff3.py
rm -rf "${HAP1_WSPACE}" "${HAP2_WSPACE}"

echo "########################################################################"
echo "SCRIPT R3b PIPELINE RUN COMPLETE."
echo "########################################################################"