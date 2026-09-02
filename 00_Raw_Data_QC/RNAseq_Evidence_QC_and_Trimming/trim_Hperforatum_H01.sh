#!/bin/bash
 
#SBATCH --job-name=trim_Hperforatum_H01                              # Job name, will show up in squeue output
#SBATCH --auks=yes                                                   # use Kerberos
#SBATCH --partition=cpu                                              # possible values: cpu, gpu
#SBATCH --ntasks=4                                                   # Number of tasks (default=1)
#SBATCH --cpus-per-task=1                                            # number of cpus for this task (default=1)
#SBATCH --mem=2G                                                     # size of memory (default 10G in cpu and gpu)
#SBATCH --nodes=1                                                    # Ensure that all cores are on one machine
#SBATCH --output=logs/job_%j.out   # File to which standard out will be written
#SBATCH --error=logs/job_%j.err    # File to which standard err will be written
 
# init module environment
. /etc/profile.d/modules.sh
 
# load your module
module load cutadapt           # default cutadapt/3.3(default)
module load clc-assembly-cell  # default clc-assembly-cell/5.1.1
 
# run your program...

# change to working directory
cd /path/to/your/directory

gunzip 2778380/1/2778380_221_S99_L001_R1_001.fastq.gz
gunzip 2778380/1/2778380_221_S99_L001_R2_001.fastq.gz
gunzip 2778380/2/2778380_221_S99_L002_R1_001.fastq.gz
gunzip 2778380/2/2778380_221_S99_L002_R2_001.fastq.gz
gunzip 2778381/1/2778381_222_S106_L001_R1_001.fastq.gz
gunzip 2778381/1/2778381_222_S106_L001_R2_001.fastq.gz
gunzip 2778381/2/2778381_222_S106_L002_R1_001.fastq.gz
gunzip 2778381/2/2778381_222_S106_L002_R2_001.fastq.gz
gunzip 2778382/1/2778382_223_S113_L001_R1_001.fastq.gz
gunzip 2778382/1/2778382_223_S113_L001_R2_001.fastq.gz
gunzip 2778382/2/2778382_223_S113_L002_R1_001.fastq.gz
gunzip 2778382/2/2778382_223_S113_L002_R2_001.fastq.gz
gunzip 2778383/1/2778383_221_E_S64_L001_R1_001.fastq.gz
gunzip 2778383/1/2778383_221_E_S64_L001_R2_001.fastq.gz
gunzip 2778383/2/2778383_221_E_S64_L002_R1_001.fastq.gz
gunzip 2778383/2/2778383_221_E_S64_L002_R2_001.fastq.gz

cat 2778380/1/2778380_221_S99_L001_R1_001.fastq 2778380/2/2778380_221_S99_L002_R1_001.fastq > Hperforatum_H01_S99_R1.fastq
cat 2778381/1/2778381_222_S106_L001_R1_001.fastq 2778381/2/2778381_222_S106_L002_R1_001.fastq > Hperforatum_H01_S106_R1.fastq
cat 2778382/1/2778382_223_S113_L001_R1_001.fastq 2778382/2/2778382_223_S113_L002_R1_001.fastq > Hperforatum_H01_S113_R1.fastq
cat 2778383/1/2778383_221_E_S64_L001_R1_001.fastq 2778383/2/2778383_221_E_S64_L002_R1_001.fastq > Hperforatum_H01_S64_R1.fastq
clc_sequence_info -k -n -r Hperforatum_H01_S99_R1.fastq > trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S106_R1.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S113_R1.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S64_R1.fastq >> trim_Hperforatum_H01.info
cat 2778380/1/2778380_221_S99_L001_R2_001.fastq 2778380/2/2778380_221_S99_L002_R2_001.fastq > Hperforatum_H01_S99_R2.fastq
cat 2778381/1/2778381_222_S106_L001_R2_001.fastq 2778381/2/2778381_222_S106_L002_R2_001.fastq > Hperforatum_H01_S106_R2.fastq
cat 2778382/1/2778382_223_S113_L001_R2_001.fastq 2778382/2/2778382_223_S113_L002_R2_001.fastq > Hperforatum_H01_S113_R2.fastq
cat 2778383/1/2778383_221_E_S64_L001_R2_001.fastq 2778383/2/2778383_221_E_S64_L002_R2_001.fastq > Hperforatum_H01_S64_R2.fastq
clc_sequence_info -k -n -r Hperforatum_H01_S99_R2.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S106_R2.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S113_R2.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S64_R2.fastq >> trim_Hperforatum_H01.info

rm 2778380/1/2778380_221_S99_L001_R1_001.fastq
rm 2778380/1/2778380_221_S99_L001_R2_001.fastq
rm 2778380/2/2778380_221_S99_L002_R1_001.fastq
rm 2778380/2/2778380_221_S99_L002_R2_001.fastq
rm 2778381/1/2778381_222_S106_L001_R1_001.fastq
rm 2778381/1/2778381_222_S106_L001_R2_001.fastq
rm 2778381/2/2778381_222_S106_L002_R1_001.fastq
rm 2778381/2/2778381_222_S106_L002_R2_001.fastq
rm 2778382/1/2778382_223_S113_L001_R1_001.fastq
rm 2778382/1/2778382_223_S113_L001_R2_001.fastq
rm 2778382/2/2778382_223_S113_L002_R1_001.fastq
rm 2778382/2/2778382_223_S113_L002_R2_001.fastq
rm 2778383/1/2778383_221_E_S64_L001_R1_001.fastq
rm 2778383/1/2778383_221_E_S64_L001_R2_001.fastq
rm 2778383/2/2778383_221_E_S64_L002_R1_001.fastq
rm 2778383/2/2778383_221_E_S64_L002_R2_001.fastq

cutadapt -e 0.1 -O 1 -a AGATCGGAAGAGC -o Hperforatum_H01_S99_R1_adapt_trimmed.fastq Hperforatum_H01_S99_R1.fastq >> trim_Hperforatum_H01.info
cutadapt -e 0.1 -O 1 -a AGATCGGAAGAGC -o Hperforatum_H01_S106_R1_adapt_trimmed.fastq Hperforatum_H01_S106_R1.fastq >> trim_Hperforatum_H01.info
cutadapt -e 0.1 -O 1 -a AGATCGGAAGAGC -o Hperforatum_H01_S113_R1_adapt_trimmed.fastq Hperforatum_H01_S113_R1.fastq >> trim_Hperforatum_H01.info
cutadapt -e 0.1 -O 1 -a AGATCGGAAGAGC -o Hperforatum_H01_S64_R1_adapt_trimmed.fastq Hperforatum_H01_S64_R1.fastq >> trim_Hperforatum_H01.info
cutadapt -e 0.1 -O 1 -a AGATCGGAAGAGC -o Hperforatum_H01_S99_R2_adapt_trimmed.fastq Hperforatum_H01_S99_R2.fastq >> trim_Hperforatum_H01.info
cutadapt -e 0.1 -O 1 -a AGATCGGAAGAGC -o Hperforatum_H01_S106_R2_adapt_trimmed.fastq Hperforatum_H01_S106_R2.fastq >> trim_Hperforatum_H01.info
cutadapt -e 0.1 -O 1 -a AGATCGGAAGAGC -o Hperforatum_H01_S113_R2_adapt_trimmed.fastq Hperforatum_H01_S113_R2.fastq >> trim_Hperforatum_H01.info
cutadapt -e 0.1 -O 1 -a AGATCGGAAGAGC -o Hperforatum_H01_S64_R2_adapt_trimmed.fastq Hperforatum_H01_S64_R2.fastq >> trim_Hperforatum_H01.info

rm Hperforatum_H01_S99_R1.fastq
rm Hperforatum_H01_S106_R1.fastq
rm Hperforatum_H01_S113_R1.fastq
rm Hperforatum_H01_S64_R1.fastq
rm Hperforatum_H01_S99_R2.fastq
rm Hperforatum_H01_S106_R2.fastq
rm Hperforatum_H01_S113_R2.fastq
rm Hperforatum_H01_S64_R2.fastq

clc_quality_trim -f 33 -c 20 -b 0.1 -l 0.9 -m 80 -r -i Hperforatum_H01_S99_R1_adapt_trimmed.fastq Hperforatum_H01_S99_R2_adapt_trimmed.fastq -o Hperforatum_H01_S99_trimmed_single.fastq -p Hperforatum_H01_S99_trimmed_pairs.fastq >> trim_Hperforatum_H01.info
clc_quality_trim -f 33 -c 20 -b 0.1 -l 0.9 -m 80 -r -i Hperforatum_H01_S106_R1_adapt_trimmed.fastq Hperforatum_H01_S106_R2_adapt_trimmed.fastq -o Hperforatum_H01_S106_trimmed_single.fastq -p Hperforatum_H01_S106_trimmed_pairs.fastq >> trim_Hperforatum_H01.info
clc_quality_trim -f 33 -c 20 -b 0.1 -l 0.9 -m 80 -r -i Hperforatum_H01_S113_R1_adapt_trimmed.fastq Hperforatum_H01_S113_R2_adapt_trimmed.fastq -o Hperforatum_H01_S113_trimmed_single.fastq -p Hperforatum_H01_S113_trimmed_pairs.fastq >> trim_Hperforatum_H01.info
clc_quality_trim -f 33 -c 20 -b 0.1 -l 0.9 -m 80 -r -i Hperforatum_H01_S64_R1_adapt_trimmed.fastq Hperforatum_H01_S64_R2_adapt_trimmed.fastq -o Hperforatum_H01_S64_trimmed_single.fastq -p Hperforatum_H01_S64_trimmed_pairs.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S99_trimmed_single.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S106_trimmed_single.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S113_trimmed_single.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S64_trimmed_single.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S99_trimmed_pairs.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S106_trimmed_pairs.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S113_trimmed_pairs.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S64_trimmed_pairs.fastq >> trim_Hperforatum_H01.info

rm Hperforatum_H01_S99_R1_adapt_trimmed.fastq
rm Hperforatum_H01_S106_R1_adapt_trimmed.fastq
rm Hperforatum_H01_S113_R1_adapt_trimmed.fastq
rm Hperforatum_H01_S64_R1_adapt_trimmed.fastq
rm Hperforatum_H01_S99_R2_adapt_trimmed.fastq
rm Hperforatum_H01_S106_R2_adapt_trimmed.fastq
rm Hperforatum_H01_S113_R2_adapt_trimmed.fastq
rm Hperforatum_H01_S64_R2_adapt_trimmed.fastq

perl splitInterleavedFastqCL.pl -in=Hperforatum_H01_S99_trimmed_pairs.fastq -out1=Hperforatum_H01_S99_R1_trimmed.fastq -out2=Hperforatum_H01_S99_R2_trimmed.fastq
perl splitInterleavedFastqCL.pl -in=Hperforatum_H01_S106_trimmed_pairs.fastq -out1=Hperforatum_H01_S106_R1_trimmed.fastq -out2=Hperforatum_H01_S106_R2_trimmed.fastq
perl splitInterleavedFastqCL.pl -in=Hperforatum_H01_S113_trimmed_pairs.fastq -out1=Hperforatum_H01_S113_R1_trimmed.fastq -out2=Hperforatum_H01_S113_R2_trimmed.fastq
perl splitInterleavedFastqCL.pl -in=Hperforatum_H01_S64_trimmed_pairs.fastq -out1=Hperforatum_H01_S64_R1_trimmed.fastq -out2=Hperforatum_H01_S64_R2_trimmed.fastq
clc_sequence_info -k -n -r Hperforatum_H01_S99_R1_trimmed.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S106_R1_trimmed.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S113_R1_trimmed.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S64_R1_trimmed.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S99_R2_trimmed.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S106_R2_trimmed.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S113_R2_trimmed.fastq >> trim_Hperforatum_H01.info
clc_sequence_info -k -n -r Hperforatum_H01_S64_R2_trimmed.fastq >> trim_Hperforatum_H01.info

pigz -p 4 Hperforatum_H01_S99_trimmed_single.fastq
pigz -p 4 Hperforatum_H01_S106_trimmed_single.fastq
pigz -p 4 Hperforatum_H01_S113_trimmed_single.fastq
pigz -p 4 Hperforatum_H01_S64_trimmed_single.fastq
pigz -p 4 Hperforatum_H01_S99_trimmed_pairs.fastq
pigz -p 4 Hperforatum_H01_S106_trimmed_pairs.fastq
pigz -p 4 Hperforatum_H01_S113_trimmed_pairs.fastq
pigz -p 4 Hperforatum_H01_S64_trimmed_pairs.fastq
pigz -p 4 Hperforatum_H01_S99_R1_trimmed.fastq
pigz -p 4 Hperforatum_H01_S106_R1_trimmed.fastq
pigz -p 4 Hperforatum_H01_S113_R1_trimmed.fastq
pigz -p 4 Hperforatum_H01_S64_R1_trimmed.fastq
pigz -p 4 Hperforatum_H01_S99_R2_trimmed.fastq
pigz -p 4 Hperforatum_H01_S106_R2_trimmed.fastq
pigz -p 4 Hperforatum_H01_S113_R2_trimmed.fastq
pigz -p 4 Hperforatum_H01_S64_R2_trimmed.fastq

