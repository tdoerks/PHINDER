#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=16G
#SBATCH -t 3-00:00:00
#SBATCH -p batch.q
#SBATCH --gres=killable:0
#SBATCH -J phinder-spades
#SBATCH -o "phinder_spades_%j.log"
#SBATCH -e "phinder_spades_%j.log"
#SBATCH --mail-user=tylerdoe@ksu.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load apptainer 2>/dev/null || module load Apptainer 2>/dev/null || true
module load nextflow  2>/dev/null || module load Nextflow  2>/dev/null || true

export NXF_OPTS='-Xms2g -Xmx8g'

PIPELINE="/fastscratch/tylerdoe/PHINDER-spades-test/main_compare.nf"
SAMPLESHEET="/fastscratch/tylerdoe/PHINDER-spades-test/samplesheets/samplesheet_spades_compare.csv"
OUTDIR="/fastscratch/tylerdoe/PHINDER-spades-test/results_spades_compare"

nextflow run "$PIPELINE" \
    -profile slurm,beocat \
    --compare_input "$SAMPLESHEET" \
    --outdir "$OUTDIR" \
    -resume
