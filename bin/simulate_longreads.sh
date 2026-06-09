#!/bin/bash
#==============================================================================
# Simulate Nanopore long reads from phage references for PHINDER dev/testing
#==============================================================================
# Tier-1 test data (see docs/LONG_READS_PLAN.md): generate controlled long
# reads with known ground truth, so the long_reads pipeline path can be wired
# and validated before sourcing real Nanopore data.
#
# Usage:  bash bin/simulate_longreads.sh <reference.fasta> <sample_id> [depth]
# Output: <sample_id>_nanopore.fastq.gz  (place in the long_reads samplesheet)
#==============================================================================
set -euo pipefail

REF="${1:?usage: simulate_longreads.sh <reference.fasta> <sample_id> [depth]}"
SAMPLE="${2:?provide a sample_id}"
DEPTH="${3:-50}"
# VERIFY container tag resolves in registry before first run
IMG="docker://quay.io/biocontainers/badread:0.4.1--pyhdfd78af_0"

echo "Simulating ~${DEPTH}x Nanopore reads from ${REF} -> ${SAMPLE}_nanopore.fastq.gz"
apptainer exec "${IMG}" badread simulate \
    --reference "${REF}" \
    --quantity "${DEPTH}x" \
    | gzip > "${SAMPLE}_nanopore.fastq.gz"

echo "Done: ${SAMPLE}_nanopore.fastq.gz"
