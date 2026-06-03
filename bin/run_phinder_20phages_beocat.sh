#!/bin/bash
#SBATCH --job-name=phinder_20phages
#SBATCH --partition=batch.q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=phinder_20phages_%j.log
#SBATCH --error=phinder_20phages_%j.err

#==============================================================================
# PHINDER Beocat - 20 Phage Test
#==============================================================================
# Runs PHINDER with 20 diverse phages from SRA
#
# Usage:
#   sbatch bin/run_phinder_20phages_beocat.sh
#
# Includes:
#   - Lambda, T4, T7 (classic model phages)
#   - Various temperate and lytic phages
#   - Different genome sizes (40-200 kb)
#   - Multiple bacterial hosts
#==============================================================================

set -euo pipefail

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------

# Input: SRA accession list (20 phages)
SRR_LIST="test_phages_20_sra.txt"

# Output directory with timestamp
OUTDIR="results_20phages_$(date +%Y%m%d_%H%M%S)"

#------------------------------------------------------------------------------
# Run Pipeline
#------------------------------------------------------------------------------

echo "========================================"
echo "  PHINDER 20-Phage Test Run"
echo "========================================"
echo ""
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "Start time: $(date)"
echo ""
echo "Configuration:"
echo "  SRR List: ${SRR_LIST}"
echo "  Output dir: ${OUTDIR}"
echo "  Total phages: 20"
echo ""

# Load Nextflow
echo "Loading modules..."
if module avail Nextflow 2>&1 | grep -q "24.04"; then
    module load Nextflow/24.04.2
elif module avail Nextflow 2>&1 | grep -q "Nextflow"; then
    module load Nextflow
else
    echo "ERROR: No Nextflow module found"
    echo "Available modules:"
    module avail Nextflow
    exit 1
fi

# Run PHINDER in SRA mode
echo ""
echo "========================================"
echo "  Starting PHINDER Pipeline"
echo "========================================"
echo ""

nextflow run main.nf \
    --input ${SRR_LIST} \
    --input_mode sra \
    --outdir ${OUTDIR} \
    -profile slurm \
    -resume \
    -with-report ${OUTDIR}/phinder_report.html \
    -with-timeline ${OUTDIR}/phinder_timeline.html

EXIT_CODE=$?

echo ""
echo "========================================"
echo "  Pipeline Complete"
echo "========================================"
echo ""
echo "Exit code: ${EXIT_CODE}"
echo "End time: $(date)"
echo ""

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "✓ SUCCESS! All phages processed"
    echo ""
    echo "Results in: ${OUTDIR}/"
    echo ""
    echo "View reports:"
    echo "  - MultiQC: ${OUTDIR}/multiqc/multiqc_report.html"
    echo "  - PHINDER Summary: ${OUTDIR}/summary/phinder_summary.html"
    echo "  - Execution: ${OUTDIR}/phinder_report.html"
    echo ""
    echo "Download report:"
    echo "  scp tylerdoe@beocat.cis.ksu.edu:$(pwd)/${OUTDIR}/summary/phinder_summary.html ."
else
    echo "✗ FAILED with exit code ${EXIT_CODE}"
    echo ""
    echo "Check logs:"
    echo "  - SLURM: phinder_20phages_${SLURM_JOB_ID}.{log,err}"
    echo "  - Nextflow: .nextflow.log"
fi

exit ${EXIT_CODE}
