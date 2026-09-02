#!/bin/bash

# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --partition=cpu        
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16     
#SBATCH --mem=32G              
#SBATCH --time=02:00:00        
#SBATCH --job-name=Hard_Mask_UnChr
#SBATCH --output=logs/hard_mask_unchr_%j.out
#SBATCH --error=logs/hard_mask_unchr_%j.err

# --- SCRIPT SETUP ---
set -eo pipefail  

# ==============================================================================
#                      DIAGNOSTICS & LOGIC CHECK SECTION
# ==============================================================================
echo "[$(date)] Starting Pre-flight Diagnostic Checks..."

source ~/.bashrc

OUTPUT_DIR="/path/to/output/directory"
INPUT_GFF_DIR="/path/to/input/directory"
INPUT_FA_DIR="/path/to/your/directory"

HAP1_GFF="${INPUT_GFF_DIR}/hap1_NUMT_NUPT.gff3"
HAP2_GFF="${INPUT_GFF_DIR}/hap2_NUMT_NUPT.gff3"
HAP1_FASTA="${INPUT_FA_DIR}/hap1.softmasked.fa"
HAP2_FASTA="${INPUT_FA_DIR}/hap2.softmasked.fa"

# Verify Essential Executables
if ! command -v python3 &> /dev/null; then
    echo "ERROR: Required executable 'python3' is missing from path environment." >&2
    exit 1
fi
echo "✔ Toolchain verification passed."

# Verify Input Files
for file in "$HAP1_GFF" "$HAP2_GFF" "$HAP1_FASTA" "$HAP2_FASTA"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Missing required input component: $file" >&2
        exit 1
    fi
done
echo "✔ All input FASTA and GFF3 files verified."

# Establish Output Workspace
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi
cd "$OUTPUT_DIR" || exit 1
echo "[$(date)] Pre-flight diagnostics passed. Executing masking and reporting payload..."

# ==============================================================================
#             IN-MEMORY COORDINATE-SAFE MASKING & REPORTING ENGINE
# ==============================================================================
python3 - << 'EOF'
import os

def process_haplotype_masking(fasta_in, gff_in, fasta_out, hap_name):
    print(f"[{hap_name}] Extracting coordinates and sequences...")
    
    # 1. Parse strict UnChr masking intervals
    intervals = []
    with open(gff_in, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            if parts[0] == 'UnChr':
                start = int(parts[3])
                end = int(parts[4])
                intervals.append((start, end))
                
    # Merge overlapping intervals to guarantee accurate position tracking
    intervals.sort()
    merged_intervals = []
    for start, end in intervals:
        if not merged_intervals or merged_intervals[-1][1] < start:
            merged_intervals.append([start, end])
        else:
            merged_intervals[-1][1] = max(merged_intervals[-1][1], end)
            
    # 2. Parse FASTA into memory structural arrays
    genome = {}
    headers = {}
    order = []
    
    with open(fasta_in, 'r') as f:
        current_id = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    genome[current_id] = bytearray("".join(seq_lines), 'ascii')
                current_id = line.split()[0][1:]
                headers[current_id] = line
                order.append(current_id)
                seq_lines = []
            else:
                seq_lines.append(line)
        if current_id:
            genome[current_id] = bytearray("".join(seq_lines), 'ascii')

    # 3. Profile baseline sequence parameters
    genome_size_before = sum(len(seq) for seq in genome.values())
    n_before = {chrid: (seq.count(b'N') + seq.count(b'n')) for chrid, seq in genome.items()}
    
    # 4. Inject hard-mask parameters across UnChr segment
    total_positions_masked = 0
    if 'UnChr' in genome:
        unchr_seq = genome['UnChr']
        for start, end in merged_intervals:
            # GFF coordinates (1-based) to Python indices (0-based)
            s_idx = start - 1
            e_idx = end
            total_positions_masked += (end - start + 1)
            unchr_seq[s_idx:e_idx] = b'N' * (end - start + 1)
    else:
        print(f"WARNING: 'UnChr' header not detected in {fasta_in}")

    # 5. Profile post-masking sequence parameters
    genome_size_after = sum(len(seq) for seq in genome.values())
    n_after = {chrid: (seq.count(b'N') + seq.count(b'n')) for chrid, seq in genome.items()}
    
    # 6. Output coordinate-invariant softmasked/hardmasked combined FASTA file
    with open(fasta_out, 'w') as f:
        for chrid in order:
            f.write(headers[chrid] + '\n')
            seq_str = genome[chrid].decode('ascii')
            for i in range(0, len(seq_str), 60):
                f.write(seq_str[i:i+60] + '\n')
                
    return {
        'size_before': genome_size_before,
        'size_after': genome_size_after,
        'masked_positions': total_positions_masked,
        'n_before': n_before,
        'n_after': n_after,
        'chromosomes': order
    }

output_dir = "/path/to/output/directory"

h1_stats = process_haplotype_masking(
    "/path/to/input/directory/hap1.softmasked.fa",
    "/path/to/input/directory/hap1_NUMT_NUPT.gff3",
    os.path.join(output_dir, "hap1.masked_unchr.fa"),
    "Haplotype 1"
)

h2_stats = process_haplotype_masking(
    "/path/to/input/directory/hap2.softmasked.fa",
    "/path/to/input/directory/hap2_NUMT_NUPT.gff3",
    os.path.join(output_dir, "hap2.masked_unchr.fa"),
    "Haplotype 2"
)

# ==============================================================================
#                      GENOMIC STATISTICAL REPORT GENERATOR
# ==============================================================================
report_path = os.path.join(output_dir, "UnChr_organellar_masking_summary.txt")
print(f"Writing integrated consolidation metrics report to: {report_path}")

with open(report_path, 'w') as r:
    r.write("========================================================================\n")
    r.write("        HyPR01 GENOME UNCHR ORGANELLAR HARD-MASKING SUMMARY REPORT      \n")
    r.write("========================================================================\n\n")
    
    for label, stats in [("HAPLOTYPE 1 (Hap1)", h1_stats), ("HAPLOTYPE 2 (Hap2)", h2_stats)]:
        r.write(f"--- {label} ---\n")
        r.write(f"Genome Size Before Masking : {stats['size_before']:,} bp\n")
        r.write(f"Genome Size After Masking  : {stats['size_after']:,} bp\n")
        r.write(f"Total Positions Hard Masked: {stats['masked_positions']:,} bp\n\n")
        
        r.write("Chromosome Breakdown:\n")
        r.write(f"{'Chromosome':<15}{'N-Count Before':<20}{'N-Count After':<20}\n")
        r.write("-" * 55 + "\n")
        for chrid in stats['chromosomes']:
            # Normalize display output targeting core pseudomolecules and UnChr blocks
            if chrid.startswith('Chr') or chrid.startswith('Ch') or chrid == 'UnChr':
                r.write(f"{chrid:<15}{stats['n_before'].get(chrid, 0):<20,}{stats['n_after'].get(chrid, 0):<20,}\n")
        r.write("\n" + "="*72 + "\n\n")

print("✔ Metrics tracking and multi-format extraction completed successfully.")
EOF

echo "[$(date)] Hard-masking pipeline workflow finished execution."