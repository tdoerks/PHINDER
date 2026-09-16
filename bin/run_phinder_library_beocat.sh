#!/bin/bash
#SBATCH --job-name=phinder_lib
#SBATCH --partition=batch.q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=phinder_lib_%j.log
#SBATCH --error=phinder_lib_%j.err

#==============================================================================
# PHINDER Beocat - Full Phage Library (86 phages, assembly mode)
#==============================================================================
# 86 phages across 9 host groups:
#   Fusobacterium (15), Porphyromonas (4), C. difficile (17),
#   Klebsiella (10), Pseudomonas (13), Staphylococcus (12),
#   Salmonella (7), Classic E. coli (7), Acinetobacter (1)
#
# Usage:
#   bash bin/download_phage_library.sh   # downloads 71 new phages
#   bash bin/download_fuso_phages_new.sh # downloads 9 new Fusobacterium
#   sbatch bin/run_phinder_library_beocat.sh
#==============================================================================

set -euo pipefail

SAMPLESHEET="samplesheet_phage_library.csv"
OUTDIR="results_phage_library_$(date +%Y%m%d_%H%M%S)"

echo "========================================"
echo "  PHINDER - Full Phage Library (86)"
echo "========================================"
echo ""
echo "Job ID:    ${SLURM_JOB_ID}"
echo "Node:      ${SLURM_NODELIST}"
echo "Start:     $(date)"
echo "Output:    ${OUTDIR}"
echo ""

# Verify all FASTAs exist
echo "Checking input FASTAs..."
MISSING=0
while IFS=, read -r sample f1 f2 assembly mode; do
    [[ "$sample" == "sample" ]] && continue
    if [[ ! -f "$assembly" ]]; then
        echo "  MISSING: ${sample} — ${assembly}"
        MISSING=$((MISSING+1))
    fi
done < "${SAMPLESHEET}"

if [[ $MISSING -gt 0 ]]; then
    echo ""
    echo "ERROR: $MISSING FASTA(s) missing. Run:"
    echo "  bash bin/download_phage_library.sh"
    echo "  bash bin/download_fuso_phages_new.sh"
    exit 1
fi
echo "  All FASTAs present ($(grep -c '^[^s]' "$SAMPLESHEET") phages)"
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
    echo "SUCCESS — results in: ${OUTDIR}/"
    echo ""
    echo "Download dashboard:"
    echo "  scp tylerdoe@beocat.cis.ksu.edu:$(pwd)/${OUTDIR}/summary/phinder_summary.html /mnt/c/Users/tdoerks/Downloads/"
else
    echo "FAILED — check phinder_lib_${SLURM_JOB_ID}.{log,err} and .nextflow.log"
fi

exit ${EXIT_CODE}
