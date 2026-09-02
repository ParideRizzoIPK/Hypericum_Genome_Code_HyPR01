#!/bin/bash
 
# --- SLURM SUBMISSION DIRECTIVES ---
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G
#SBATCH --time=24:00:00
#SBATCH --job-name=yahs_jbat_final
#SBATCH --output=logs/jbat_final_%j.out
#SBATCH --error=logs/jbat_final_%j.err
 
# --- SCRIPT SETUP ---
set -e
set -o pipefail
 
# --- 1. CONFIGURATION ---
OUT_DIR="/path/to/output/directory"
JUICER_TOOLS_JAR="/path/to/your/directory/juicer_tools_1.22.01.jar"
 
# Modules (Standard tools only)
module load samtools/1.23.1
module load yahs/1.2
 
cd "${OUT_DIR}"
 
# --- 2. AUTO-DETECT YAHS JUICER TOOL ---
YAHS_BIN_DIR=$(dirname $(which yahs))
if [ -f "${YAHS_BIN_DIR}/juicer" ]; then
    JUICER_PRE_CMD="${YAHS_BIN_DIR}/juicer"
elif [ -f "${YAHS_BIN_DIR}/juicer_pre" ]; then
    JUICER_PRE_CMD="${YAHS_BIN_DIR}/juicer_pre"
else
    JUICER_PRE_CMD="/opt/Bio/yahs/1.2/bin/juicer"
fi
echo "Using Juicer conversion tool: ${JUICER_PRE_CMD}"
 
# --- 3. EMBEDDED PERL SCRIPTS ---
 
# Script A: Create Identity AGP (Maps contig -> contig)
cat << 'EOF' > make_identity_agp.pl
#!/usr/bin/perl
use strict;
use warnings;
my $fai_file = $ARGV[0];
open(my $fh, '<', $fai_file) or die "Cannot open $fai_file: $!";
while (<$fh>) {
    chomp;
    my ($ctg, $len) = split(/\t/);
    # AGP: scaffold start end part type ctg start end orient
    print "$ctg\t1\t$len\t1\tW\t$ctg\t1\t$len\t+\n";
}
close($fh);
EOF
 
# Script B: Convert YaHS AGP to valid Juicebox .assembly format
cat << 'EOF' > agp2assembly.pl
#!/usr/bin/perl
use strict;
use warnings;
my ($agp_file, $fai_file, $out_file) = @ARGV;
 
my %ctg_map;   # name -> ID
my %ctg_lens;  # name -> length
my $id_counter = 1;
 
# 1. Read FAI
open(my $fai_fh, '<', $fai_file) or die "Cannot open $fai_file: $!";
while (<$fai_fh>) {
    chomp;
    my ($name, $len) = split(/\t/);
    $ctg_map{$name} = $id_counter;
    $ctg_lens{$name} = $len;
    $id_counter++;
}
close($fai_fh);
 
# 2. Parse AGP
my %scaffolds;      
my @scaffold_order; 
 
open(my $agp_fh, '<', $agp_file) or die "Cannot open $agp_file: $!";
while (<$agp_fh>) {
    next if /^#/;
    next if /^\s*$/;
    my @parts = split(/\s+/);
    next if ($parts[4] eq 'N');
    
    my $scaff = $parts[0];
    my $ctg   = $parts[5];
    my $orient = $parts[8];
    
    if (!exists $scaffolds{$scaff}) {
        push @scaffold_order, $scaff;
        $scaffolds{$scaff} = [];
    }
    
    my $sign = ($orient eq '+') ? "" : "-";
    my $cid = $ctg_map{$ctg};
    if (defined $cid) {
        push @{$scaffolds{$scaff}}, "$sign$cid";
    }
}
close($agp_fh);
 
# 3. Write Output
open(my $out_fh, '>', $out_file) or die "Cannot open $out_file: $!";
# Header: >ctg_name ID length
foreach my $ctg (sort { $ctg_map{$a} <=> $ctg_map{$b} } keys %ctg_map) {
    print $out_fh ">$ctg $ctg_map{$ctg} $ctg_lens{$ctg}\n";
}
# Body: Scaffolds
foreach my $scaff (@scaffold_order) {
    my $comps = $scaffolds{$scaff};
    if (@$comps) {
        print $out_fh join(" ", @$comps) . "\n";
    }
}
close($out_fh);
EOF
 
# --- 4. PROCESSING FUNCTION ---
run_jbat_fix() {
    local HAP=$1
    local PREFIX="HyPR-01_v3.hap${HAP}"
    
    echo "========================================================"
    echo ">>>> MANUAL JBAT FIX FOR HAPLOTYPE ${HAP} <<<<"
    echo "========================================================"
    
    # Inputs
    local BIN_FILE="${PREFIX}.bin"
    local FINAL_AGP="${PREFIX}_scaffolds_final.agp"
    local CTG_FAI="HyPR01_hap${HAP}.fa.fai"
    
    # Check Inputs
    if [ ! -f "${BIN_FILE}" ]; then echo "ERROR: BIN missing"; exit 1; fi
    if [ ! -f "${FINAL_AGP}" ]; then echo "ERROR: AGP missing"; exit 1; fi
    if [ ! -f "${CTG_FAI}" ]; then echo "ERROR: FAI missing"; exit 1; fi
    
    # Outputs
    local IDENTITY_AGP="${PREFIX}.identity.agp"
    local CTG_SIZES="${PREFIX}.contigs.chrom.sizes"
    local JBAT_PREFIX="${PREFIX}_JBAT"
    local RAW_TXT="${JBAT_PREFIX}.raw.txt"
    local SORTED_TXT="${JBAT_PREFIX}.sorted.txt"
    local OUT_HIC="${JBAT_PREFIX}.hic"
    local OUT_ASSEMBLY="${JBAT_PREFIX}.assembly"
    local OUT_LIFTOVER="${JBAT_PREFIX}.liftover.agp"
    
    # 1. Create Helper Files
    echo "[HAP${HAP}] Creating identity AGP and chrom.sizes..."
    cut -f1,2 "${CTG_FAI}" > "${CTG_SIZES}"
    perl make_identity_agp.pl "${CTG_FAI}" > "${IDENTITY_AGP}"
    
    # 2. Generate .assembly File (Correct Format)
    echo "[HAP${HAP}] Converting AGP to .assembly..."
    perl agp2assembly.pl "${FINAL_AGP}" "${CTG_FAI}" "${OUT_ASSEMBLY}"
    cp "${FINAL_AGP}" "${OUT_LIFTOVER}"
    
    # 3. Extract Contig-Level Contacts to File
    # We write to DISK to avoid Java pipe issues
    echo "[HAP${HAP}] Extracting and sorting contacts (this may take time)..."
    
    # Step A: Convert BIN to Text using Identity AGP (forces contig names)
    ${JUICER_PRE_CMD} pre "${BIN_FILE}" "${IDENTITY_AGP}" "${CTG_FAI}" > "${RAW_TXT}"
    
    # Step B: Sort
    sort -k2,2d -k6,6d -T "${OUT_DIR}" --parallel=${SLURM_CPUS_PER_TASK} -S64G "${RAW_TXT}" > "${SORTED_TXT}"
    
    # 4. Generate .hic from File (No pipes!)
    echo "[HAP${HAP}] Building .hic matrix..."
    
    java -Xmx160G -jar "${JUICER_TOOLS_JAR}" pre \
        -j ${SLURM_CPUS_PER_TASK} \
        "${SORTED_TXT}" \
        "${OUT_HIC}" \
        "${CTG_SIZES}"
        
    if [ -f "${OUT_HIC}" ]; then
        echo "[HAP${HAP}] SUCCESS: Created ${OUT_HIC}"
        # Cleanup large intermediates
        rm "${RAW_TXT}" "${SORTED_TXT}" "${IDENTITY_AGP}" "${CTG_SIZES}"
    else
        echo "[HAP${HAP}] ERROR: Failed to create .hic file"
        exit 1
    fi
}
 
# --- EXECUTE ---
run_jbat_fix "1"
run_jbat_fix "2"
 
# Cleanup scripts
rm make_identity_agp.pl agp2assembly.pl
 
echo "========================================================"
echo "MANUAL FIX COMPLETE"
echo "Load *_JBAT.assembly in Juicebox using 'File -> Open Assembly'"
echo "========================================================"