#!/bin/bash
#SBATCH --job-name=phinder_sra
#SBATCH --partition=ksu-gen-bio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=phinder_sra_%j.log
#SBATCH --error=phinder_sra_%j.err

#==============================================================================
# PHINDER Beocat - SRA Mode (Quick Test)
#==============================================================================
# Runs PHINDER with direct SRA download (like COMPASS)
#
# Usage:
#   sbatch bin/run_phinder_sra_beocat.sh
#
# This will download Lambda phage (SRR5131134) and run full pipeline
#==============================================================================

set -euo pipefail

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------

# Input: SRA accession list
SRR_LIST="test_lambda_sra.txt"

# Output directory with timestamp
OUTDIR="results_lambda_$(date +%Y%m%d_%H%M%S)"

#------------------------------------------------------------------------------
# Run Pipeline
#------------------------------------------------------------------------------

echo "========================================"
echo "  PHINDER Pipeline - SRA Mode"
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

# Load Nextflow
echo "Loading modules..."
module load Nextflow/24.04.4

# Run PHINDER in SRA mode
echo ""
echo "========================================"
echo "  Starting PHINDER (SRA Mode)"
echo "========================================"
echo ""
echo "Pipeline will:"
echo "  1. Download SRR5131134 (Lambda phage)"
echo "  2. Run FastQC and fastp"
echo "  3. Assemble with Unicycler"
echo "  4. Annotate with Pharokka"
echo "  5. Assess quality with CheckV"
echo "  6. Predict lifestyle with VIBRANT"
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
    echo "✓ SUCCESS!"
    echo ""
    echo "Results in: ${OUTDIR}/"
    echo ""
    echo "Key outputs:"
    echo "  - Assembly: ${OUTDIR}/assemblies/SRR5131134_assembly.fasta"
    echo "  - CheckV: ${OUTDIR}/checkv/SRR5131134_checkv/"
    echo "  - Pharokka: ${OUTDIR}/pharokka/SRR5131134_pharokka/"
    echo "  - MultiQC: ${OUTDIR}/multiqc/multiqc_report.html"
    echo ""
    echo "Download MultiQC report to view:"
    echo "  scp tylerdoe@beocat.cis.ksu.edu:$(pwd)/${OUTDIR}/multiqc/multiqc_report.html ."
else
    echo "✗ FAILED with exit code ${EXIT_CODE}"
    echo ""
    echo "Check logs:"
    echo "  - SLURM: phinder_sra_${SLURM_JOB_ID}.{log,err}"
    echo "  - Nextflow: .nextflow.log"
    echo "  - Work dir: work/"
fi

exit ${EXIT_CODE}
