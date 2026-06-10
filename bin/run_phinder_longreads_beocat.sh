#!/bin/bash
#SBATCH --job-name=phinder_lr_smoke
#SBATCH --partition=batch.q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=phinder_lr_smoke_%j.log
#SBATCH --error=phinder_lr_smoke_%j.err
#SBATCH --mail-user=tylerdoe@ksu.edu
#SBATCH --mail-type=END,FAIL
#==============================================================================
# PHINDER long-read (Phase 1) test on Beocat
#==============================================================================
# Tier-1 data: simulated Nanopore reads from the Lambda reference (J02459),
# so there is a known ground truth (1 contig, ~48.5 kb) to assemble against.
#
# By default this runs a SMOKE TEST of the new path only:
#   NanoPlot -> Flye -> Medaka -> a clean assembly
# (downstream annotation skipped). Flip SKIP_DOWNSTREAM=false for the full run.
#
# Run from the repo root:
#   sbatch bin/run_phinder_longreads_beocat.sh
#==============================================================================

REF=lambda_ref.fasta
SAMPLE=Lambda
READS=${SAMPLE}_nanopore.fastq.gz
SHEET=sheet_lr_lambda.csv
OUTDIR=results_lr_smoke
SKIP_DOWNSTREAM=true   # set to false to run CheckV/Pharokka/VIBRANT/DIAMOND/Phanotate too

module load Nextflow/24.04.2 2>/dev/null || module load Nextflow

# --- Tier-1 test data (idempotent: only builds what is missing) ---
if [ ! -s "${REF}" ]; then
    echo "Fetching Lambda reference (J02459)..."
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=J02459&rettype=fasta&retmode=text" > "${REF}"
fi
if [ ! -s "${READS}" ]; then
    echo "Simulating ~50x Nanopore reads with Badread..."
    bash bin/simulate_longreads.sh "${REF}" "${SAMPLE}"
fi
printf 'sample,long_reads,platform\n%s,%s,nanopore\n' "${SAMPLE}" "${READS}" > "${SHEET}"

# --- run the long-read pipeline ---
SKIP_ARGS=""
if [ "${SKIP_DOWNSTREAM}" = "true" ]; then
    SKIP_ARGS="--skip_checkv --skip_pharokka --skip_vibrant --skip_diamond --skip_phanotate"
fi

nextflow run main.nf -profile local \
    --input_mode long_reads --input "${SHEET}" \
    ${SKIP_ARGS} \
    --outdir "${OUTDIR}" \
    -resume   # reuse cached steps (e.g. a prior successful Flye) across reruns
