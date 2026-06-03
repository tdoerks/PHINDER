#!/bin/bash
#SBATCH --job-name=phinder_3phages
#SBATCH --partition=batch.q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=phinder_3phages_%j.log
#SBATCH --error=phinder_3phages_%j.err

#==============================================================================
# PHINDER Beocat - 3 Phage Test (Conservative Resources)
#==============================================================================
# Runs PHINDER with Lambda, T4, and T7 phages
# Uses conservative resource settings to avoid SLURM termination
#==============================================================================

set -euo pipefail

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------

# Input: SRA accession list (3 classic phages)
SRR_LIST="test_phages_sra.txt"

# Output directory with timestamp
OUTDIR="results_3phages_$(date +%Y%m%d_%H%M%S)"

#------------------------------------------------------------------------------
# Run Pipeline
#------------------------------------------------------------------------------

echo "========================================"
echo "  PHINDER 3-Phage Test (Conservative)"
echo "========================================"
echo ""
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "Start time: $(date)"
echo ""
echo "Configuration:"
echo "  SRR List: ${SRR_LIST}"
echo "  Output dir: ${OUTDIR}"
echo "  Phages: Lambda, T4, T7"
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
    -work-dir work_3phages \
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
    echo "Key outputs:"
    echo "  - Lambda: ${OUTDIR}/assemblies/SRR5131134_assembly.fasta"
    echo "  - T4: ${OUTDIR}/assemblies/SRR5131135_assembly.fasta"
    echo "  - T7: ${OUTDIR}/assemblies/SRR5131136_assembly.fasta"
    echo ""
    echo "View reports:"
    echo "  - PHINDER Summary: ${OUTDIR}/summary/phinder_summary.html"
    echo "  - MultiQC: ${OUTDIR}/multiqc/multiqc_report.html"
else
    echo "✗ FAILED with exit code ${EXIT_CODE}"
    echo ""
    echo "Check logs:"
    echo "  - SLURM: phinder_3phages_${SLURM_JOB_ID}.{log,err}"
    echo "  - Nextflow: .nextflow.log"
fi

exit ${EXIT_CODE}
