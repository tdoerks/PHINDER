#!/bin/bash
#SBATCH --job-name=phinder_fuso
#SBATCH --partition=batch.q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --output=phinder_fuso_%j.log
#SBATCH --error=phinder_fuso_%j.err

#==============================================================================
# PHINDER Beocat - All Fusobacterium Phages (Assembly Mode)
#==============================================================================
# Annotates 15 Fusobacterium phage genomes from NCBI GenBank
# (pre-assembled FASTAs — skips QC/assembly steps)
#
# Usage:
#   bash bin/download_fuso_phages_new.sh   # download 9 new FASTAs first
#   sbatch bin/run_phinder_fuso_all_beocat.sh
#
# Phages (15 total):
#   Original 6 (OR492271-276): phiFN37, phiFN38, phiBB, phiPaco, phiHugo, phiKSUM
#   New 9: Fnu1, FNU2, FNU3, TCUFN6, TCUFN3, JDFnp1, JDFnp6, JDFnp4, JDFnp7
#==============================================================================

set -euo pipefail

SAMPLESHEET="samplesheet_fuso_all.csv"
OUTDIR="results_fuso_all_$(date +%Y%m%d_%H%M%S)"

echo "========================================"
echo "  PHINDER - All Fusobacterium Phages"
echo "========================================"
echo ""
echo "Job ID:    ${SLURM_JOB_ID}"
echo "Node:      ${SLURM_NODELIST}"
echo "Start:     $(date)"
echo "Samplesheet: ${SAMPLESHEET}"
echo "Output:    ${OUTDIR}"
echo ""

# Verify FASTAs exist
echo "Checking input FASTAs..."
MISSING=0
while IFS=, read -r sample f1 f2 assembly mode; do
    [[ "$sample" == "sample" ]] && continue
    if [[ ! -f "$assembly" ]]; then
        echo "  MISSING: ${sample} — ${assembly}"
        MISSING=$((MISSING+1))
    else
        echo "  OK: ${sample} ($(wc -c < "$assembly" | numfmt --to=iec))"
    fi
done < "${SAMPLESHEET}"

if [[ $MISSING -gt 0 ]]; then
    echo ""
    echo "ERROR: $MISSING FASTA(s) missing. Run:"
    echo "  bash bin/download_fuso_phages_new.sh"
    exit 1
fi
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
    -profile slurm,beocat \
    -resume \
    -with-report "${OUTDIR}/phinder_report.html" \
    -with-timeline "${OUTDIR}/phinder_timeline.html"

EXIT_CODE=$?

echo ""
echo "Exit code: ${EXIT_CODE}"
echo "End time:  $(date)"

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "SUCCESS - results in: ${OUTDIR}/"
    echo "Download dashboard:"
    echo "  scp tylerdoe@beocat.cis.ksu.edu:$(pwd)/${OUTDIR}/summary/phinder_summary.html /mnt/c/Users/tdoerks/Downloads/"
else
    echo "FAILED - check phinder_fuso_${SLURM_JOB_ID}.{log,err} and .nextflow.log"
fi

exit ${EXIT_CODE}
