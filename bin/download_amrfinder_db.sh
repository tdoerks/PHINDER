#!/usr/bin/env bash
# Download AMRFinder Plus database to Beocat
# Run once from the PHINDER-dev login node before submitting jobs
# Uses the same container as the pipeline module

set -euo pipefail

DB_DIR="/fastscratch/tylerdoe/databases/amrfinder_db"
mkdir -p "$DB_DIR"

echo "Downloading AMRFinder Plus database to: $DB_DIR"
echo "This may take a few minutes (~300 MB)..."
echo ""

apptainer exec \
    --bind "$DB_DIR:$DB_DIR" \
    docker://staphb/ncbi-amrfinderplus:4.2.7-2026-08-07.1 \
    amrfinder --update --database "$DB_DIR"

echo ""
echo "Done. Database contents:"
ls -lh "$DB_DIR"
echo ""
echo "Add to nextflow.config or run with:"
echo "  --amrfinder_db $DB_DIR"
