#!/bin/bash
#SBATCH --job-name=phinder_full
#SBATCH --partition=batch.q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --output=phinder_full_%j.log
#SBATCH --error=phinder_full_%j.err

#==============================================================================
# PHINDER Beocat - Full Test (3 Phages)
#==============================================================================
# Runs PHINDER with Lambda, T4, and T7 phages
#
# Usage:
#   sbatch bin/run_phinder_full_test_beocat.sh
#
# Downloads and analyzes:
#   - Lambda phage (SRR5131134) - 48.5 kb, temperate
#   - T4 phage (SRR5131135) - 169 kb, lytic
#   - T7 phage (SRR5131136) - 40 kb, lytic
#==============================================================================

set -euo pipefail

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------

# Input: SRA accession list (3 phages)
SRR_LIST="test_phages_sra.txt"

# Output directory with timestamp
OUTDIR="results_full_test_$(date +%Y%m%d_%H%M%S)"

#------------------------------------------------------------------------------
# Run Pipeline
#------------------------------------------------------------------------------

echo "========================================"
echo "  PHINDER Full Test - 3 Phages"
echo "========================================"
echo ""
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "Start time: $(date)"
echo ""
echo "Configuration:"
echo "  SRR List: ${SRR_LIST}"
echo "  Output dir: ${OUTDIR}"
echo ""
echo "Phages to analyze:"
echo "  1. Lambda (SRR5131134) - ~48.5 kb, temperate"
echo "  2. T4 (SRR5131135) - ~169 kb, lytic"
echo "  3. T7 (SRR5131136) - ~40 kb, lytic"
echo ""

# Load Nextflow
echo "Loading modules..."
module load Nextflow/24.04.4

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
    echo "✓ SUCCESS! All 3 phages analyzed"
    echo ""
    echo "Results in: ${OUTDIR}/"
    echo ""
    echo "Assemblies:"
    echo "  - Lambda: ${OUTDIR}/assemblies/SRR5131134_assembly.fasta"
    echo "  - T4: ${OUTDIR}/assemblies/SRR5131135_assembly.fasta"
    echo "  - T7: ${OUTDIR}/assemblies/SRR5131136_assembly.fasta"
    echo ""
    echo "View MultiQC report:"
    echo "  scp tylerdoe@beocat.cis.ksu.edu:$(pwd)/${OUTDIR}/multiqc/multiqc_report.html ."
else
    echo "✗ FAILED with exit code ${EXIT_CODE}"
    echo ""
    echo "Check logs:"
    echo "  - SLURM: phinder_full_${SLURM_JOB_ID}.{log,err}"
    echo "  - Nextflow: .nextflow.log"
fi

exit ${EXIT_CODE}
