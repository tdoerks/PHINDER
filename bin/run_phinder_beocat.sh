#!/bin/bash
#SBATCH --job-name=phinder_test
#SBATCH --partition=batch.q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=phinder_%j.log
#SBATCH --error=phinder_%j.err

#==============================================================================
# PHINDER Beocat Submission Script
#==============================================================================
# This script runs PHINDER on Beocat HPC with SLURM scheduler
#
# Usage:
#   sbatch bin/run_phinder_beocat.sh
#
# Before running:
#   1. Edit SAMPLESHEET variable below
#   2. Verify database paths in nextflow.config
#   3. Download test data with: ./bin/download_test_phages.sh
#==============================================================================

set -euo pipefail

#------------------------------------------------------------------------------
# Configuration - EDIT THESE
#------------------------------------------------------------------------------

# Input samplesheet (CSV file)
SAMPLESHEET="test_samplesheet.csv"

# Output directory
OUTDIR="results_$(date +%Y%m%d_%H%M%S)"

# Input mode: 'reads' or 'assembly'
INPUT_MODE="reads"

#------------------------------------------------------------------------------
# DO NOT EDIT BELOW THIS LINE
#------------------------------------------------------------------------------

echo "========================================"
echo "  PHINDER Pipeline - Beocat Run"
echo "========================================"
echo ""
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "Start time: $(date)"
echo ""
echo "Configuration:"
echo "  Samplesheet: ${SAMPLESHEET}"
echo "  Output dir: ${OUTDIR}"
echo "  Input mode: ${INPUT_MODE}"
echo ""

# Load required modules
echo "Loading modules..."
module load Nextflow/24.04.4

# Verify samplesheet exists
if [ ! -f "${SAMPLESHEET}" ]; then
    echo "ERROR: Samplesheet not found: ${SAMPLESHEET}"
    echo "Please create samplesheet or run: ./bin/download_test_phages.sh"
    exit 1
fi

# Run PHINDER
echo ""
echo "========================================"
echo "  Starting PHINDER Pipeline"
echo "========================================"
echo ""

nextflow run main.nf \
    --input ${SAMPLESHEET} \
    --input_mode ${INPUT_MODE} \
    --outdir ${OUTDIR} \
    -profile slurm \
    -resume \
    -with-report ${OUTDIR}/phinder_report.html \
    -with-timeline ${OUTDIR}/phinder_timeline.html \
    -with-dag ${OUTDIR}/phinder_dag.png

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
    echo "✓ SUCCESS!"
    echo ""
    echo "Results available in: ${OUTDIR}/"
    echo ""
    echo "Key outputs:"
    echo "  - MultiQC report: ${OUTDIR}/multiqc/multiqc_report.html"
    echo "  - Assemblies: ${OUTDIR}/assemblies/"
    echo "  - Annotations: ${OUTDIR}/pharokka/"
    echo "  - Quality: ${OUTDIR}/checkv/"
    echo ""
    echo "Nextflow reports:"
    echo "  - Execution report: ${OUTDIR}/phinder_report.html"
    echo "  - Timeline: ${OUTDIR}/phinder_timeline.html"
    echo "  - DAG: ${OUTDIR}/phinder_dag.png"
else
    echo "✗ FAILED with exit code ${EXIT_CODE}"
    echo ""
    echo "Check logs:"
    echo "  - SLURM output: phinder_${SLURM_JOB_ID}.log"
    echo "  - SLURM errors: phinder_${SLURM_JOB_ID}.err"
    echo "  - Nextflow log: .nextflow.log"
    echo ""
    echo "Common issues:"
    echo "  1. Database paths incorrect in nextflow.config"
    echo "  2. Missing input files"
    echo "  3. Insufficient memory/time allocation"
fi

exit ${EXIT_CODE}
