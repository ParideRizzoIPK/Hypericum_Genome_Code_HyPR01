#!/bin/bash
#SBATCH --job-name=Table4_Organelle_Assembly_Stats
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --partition=cpu
#SBATCH --output=logs/table4_organelle_stats_%j.out
#SBATCH --error=logs/table4_organelle_stats_%j.err

set -e -u -o pipefail

# =============================================================================
# Environment wrapper for build_table4_organelle_assembly_stats.py (this
# same directory). Loads the required modules, provisions a biopython venv
# on first run, then executes the script -- see that script's own docstring
# for what it computes and how each row is derived.
# =============================================================================

module purge
module load minimap2/2.24
module load samtools/1.23.1
module load mosdepth/0.2.6
module load gfatools
module load python/3.12.0

export SINGULARITY_BIND="/path/to/your/directory"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="${SCRIPT_DIR}/.venv_table4_organelle_stats"

if [ ! -x "${VENV_DIR}/bin/python" ]; then
    echo "INFO: Creating venv for biopython at ${VENV_DIR}"
    python -m venv "${VENV_DIR}"
    "${VENV_DIR}/bin/pip" install --quiet --upgrade pip
    "${VENV_DIR}/bin/pip" install --quiet biopython
else
    echo "INFO: Reusing existing venv at ${VENV_DIR}"
fi

echo "INFO: Running build_table4_organelle_assembly_stats.py"
"${VENV_DIR}/bin/python" "${SCRIPT_DIR}/build_table4_organelle_assembly_stats.py"
