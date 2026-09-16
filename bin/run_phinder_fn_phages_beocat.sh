#!/bin/bash
#SBATCH --job-name=phinder_fn
#SBATCH --partition=batch.q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=phinder_fn_%j.log
#SBATCH --error=phinder_fn_%j.err

#==============================================================================
# PHINDER Beocat - Fusobacterium Phages (Assembly Mode)
#==============================================================================
# Annotates 6 Fusobacterium phage genomes from NCBI GenBank
# (pre-assembled FASTAs — skips QC/assembly steps)
#
# Usage:
#   bash bin/download_fn_phages.sh   # download FASTAs first
#   sbatch bin/run_phinder_fn_phages_beocat.sh
#
# Phages:
#   phiFN37  (OR492271) - 85 kb
#   phiRTG5  (OR492272) - 37 kb
#   phiBB    (OR492273) - 36 kb
#   phiPaco  (OR492274) - 43 kb
#   phiHugo  (OR492275) - 35 kb
#   phiKSUM  (OR492276) - 111 kb
#==============================================================================

set -euo pipefail

SAMPLESHEET="samplesheet_fn_phages.csv"
OUTDIR="results_fn_phages_$(date +%Y%m%d_%H%M%S)"

echo "========================================"
echo "  PHINDER - Fusobacterium Phages"
echo "========================================"
echo ""
echo "Job ID:    ${SLURM_JOB_ID}"
echo "Node:      ${SLURM_NODELIST}"
echo "Start:     $(date)"
echo "Samplesheet: ${SAMPLESHEET}"
echo "Output:    ${OUTDIR}"
echo ""

# Verify FASTAs exist before submitting
echo "Checking input FASTAs..."
while IFS=, read -r sample assembly; do
    [[ "$sample" == "sample" ]] && continue
    if [[ ! -f "$assembly" ]]; then
        echo "ERROR: Missing FASTA for ${sample}: ${assembly}"
        echo "Run: bash bin/download_fn_phages.sh"
        exit 1
    fi
    echo "  OK: ${sample} ($(wc -c < "$assembly" | numfmt --to=iec))"
done < "${SAMPLESHEET}"
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
    echo "Download report:"
    echo "  scp tylerdoe@beocat.cis.ksu.edu:$(pwd)/${OUTDIR}/phinder_report.html ."
else
    echo "FAILED - check phinder_fn_${SLURM_JOB_ID}.{log,err} and .nextflow.log"
fi

exit ${EXIT_CODE}
