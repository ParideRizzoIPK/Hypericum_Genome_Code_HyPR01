#!/bin/bash
#SBATCH --job-name=Hifiasm_HiC_Assembly
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32  # Hifiasm can use multiple threads
#SBATCH --mem=512G            # Memory requirement depends heavily on genome size, complexity, and Hi-C data.
                              # May need adjustment.
#SBATCH --time=7-00:00:00     # Assembly can be time-consuming. Adjust as needed.
#SBATCH --output=logs/hifiasm_hic_%j.out
#SBATCH --error=logs/hifiasm_hic_%j.err
#SBATCH --partition=cpu # Specify your Slurm partition

# --- Load necessary modules ---
module load hifiasm # Or the specific module name for Hifiasm on your cluster

# --- Define variables ---
HIFI_READS="/path/to/input/directory/hypr01_hifi_reads.fastq"  # Your input PacBio CCS/HiFi fastq file
HIC_R1_READS="/path/to/input/directory/2037511_Hypericum_PR1_HiC166_S68_L002_R1_001.fastq.gz" # Your Hi-C R1 reads
HIC_R2_READS="/path/to/input/directory/2037511_Hypericum_PR1_HiC166_S68_L002_R2_001.fastq.gz" # Your Hi-C R2 reads

OUTPUT_PREFIX="HyPR01_genome_hic_assembly" # Desired prefix for output files
THREADS=32

# --- Run Hifiasm with Hi-C data ---
echo "Starting Hifiasm assembly with Hi-C integration..."
hifiasm -o ${OUTPUT_PREFIX} -t ${THREADS} --h1 ${HIC_R1_READS} --h2 ${HIC_R2_READS} ${HIFI_READS}

# Note: Older versions of hifiasm might use -1 and -2 for Hi-C reads.
# Check your hifiasm version's help message (hifiasm -h) if --h1/--h2 doesn't work.
# For example: hifiasm -o ${OUTPUT_PREFIX} -t ${THREADS} -1 ${HIC_R1_READS} -2 ${HIC_R2_READS} ${HIFI_READS}


echo "Hifiasm assembly with Hi-C finished."
echo "Primary assembly graph: ${OUTPUT_PREFIX}.p_ctg.gfa (and .p_ctg.fa)"
echo "Check hifiasm documentation for specific Hi-C related outputs and logs."
# Hifiasm with Hi-C will produce haplotype-resolved assemblies.
# You will typically get *.hap1.p_ctg.gfa and *.hap2.p_ctg.gfa files,
# or a combination of primary (*.p_ctg.fa) and alternate (*.a_ctg.fa) contigs
# that are better phased and scaffolded thanks to Hi-C.