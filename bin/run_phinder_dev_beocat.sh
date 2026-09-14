#!/bin/bash
#SBATCH --job-name=phinder_dev
#SBATCH --partition=batch.q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=phinder_dev_%j.log
#SBATCH --error=phinder_dev_%j.err

#==============================================================================
# PHINDER-dev test run — 6 Fusobacterium phages (assembly mode)
# Run from: /fastscratch/tylerdoe/PHINDER-dev/
# Tests new modules before merging to main pipeline
#==============================================================================

set -euo pipefail

SAMPLESHEET="samplesheet_fn_phages.csv"
OUTDIR="results_dev_$(date +%Y%m%d_%H%M%S)"

echo "========================================"
echo "  PHINDER-dev Test Run (6 phages)"
echo "========================================"
echo ""
echo "Job ID:    ${SLURM_JOB_ID}"
echo "Node:      ${SLURM_NODELIST}"
echo "Start:     $(date)"
echo "Output:    ${OUTDIR}"
echo ""

# Load Nextflow
if module avail Nextflow 2>&1 | grep -q "24.04"; then
    module load Nextflow/24.04.2
elif module avail Nextflow 2>&1 | grep -q "Nextflow"; then
    module load Nextflow
else
    echo "ERROR: No Nextflow module found"; exit 1
fi

echo "========================================"
echo "  Starting PHINDER Pipeline"
echo "========================================"
echo ""

nextflow run main.nf \
    --input "${SAMPLESHEET}" \
    --input_mode assembly \
    --outdir "${OUTDIR}" \
    -profile slurm \
    -resume

EXIT_CODE=$?

echo ""
echo "Exit code: ${EXIT_CODE}"
echo "End time:  $(date)"

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "SUCCESS — results in: ${OUTDIR}/"
    echo "Download dashboard:"
    echo "  scp tylerdoe@beocat.cis.ksu.edu:$(pwd)/${OUTDIR}/summary/phinder_summary.html /mnt/c/Users/tdoerks/Downloads/"
else
    echo "FAILED — check phinder_dev_${SLURM_JOB_ID}.{log,err} and .nextflow.log"
fi

exit ${EXIT_CODE}
