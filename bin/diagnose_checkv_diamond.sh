#!/bin/bash
#==============================================================================
# Diagnose CheckV DIAMOND completeness failure
#==============================================================================
# CheckV's completeness step fails at "[3/8] Running DIAMOND blastp search...
# DIAMOND task failed. Program should be rerun." while contamination (HMMER)
# succeeds. DB is readable (dbinfo OK) and DIAMOND is 2.1.8, so the leading
# hypothesis is that DIAMOND's memory-mapped temp files fail on the NFS-mounted
# work dir (/fastscratch). This script reproduces the exact blastp by hand and
# A/B tests the temp location: local /tmp vs NFS /fastscratch.
#
# Run from the PHINDER repo root (e.g. /fastscratch/$USER/PHINDER):
#   bash bin/diagnose_checkv_diamond.sh [WORKDIR]
# WORKDIR defaults to work_3phages.
#==============================================================================
set -u

CHECKV_DB="${CHECKV_DB:-/fastscratch/${USER}/databases/checkv-db-v1.5}"
DB="${CHECKV_DB}/genome_db/checkv_reps.dmnd"
IMG="docker://quay.io/biocontainers/checkv:1.0.2--pyhdfd78af_0"
WORKDIR="${1:-work_3phages}"

echo "Searching for a CheckV proteins.faa under ${WORKDIR} ..."
Q=$(find "${WORKDIR}" -name proteins.faa -path "*checkv*" 2>/dev/null | head -1)
if [[ -z "${Q}" ]]; then
  echo "ERROR: no CheckV proteins.faa found under ${WORKDIR}"
  echo "       (run after at least one CheckV task has reached completeness)"
  exit 1
fi
echo "Query proteins : ${Q}"
echo "Database       : ${DB}"
if [[ ! -f "${DB}" ]]; then echo "ERROR: db not found at ${DB}"; exit 1; fi
echo

run_test () {
  local label="$1" tmp="$2"
  echo "================ ${label} (tmpdir=${tmp}) ================"
  rm -rf "${tmp}"; mkdir -p "${tmp}"
  apptainer exec "${IMG}" diamond blastp \
    --query "${Q}" --db "${DB}" \
    --out "${tmp}/aai.tsv" --outfmt 6 --threads 8 \
    --tmpdir "${tmp}" 2>&1 | tail -15
  local rc=${PIPESTATUS[0]}
  echo ">>> EXIT ${label}: ${rc}"
  if [[ ${rc} -eq 0 && -f "${tmp}/aai.tsv" ]]; then
    echo ">>> alignments produced: $(wc -l < "${tmp}/aai.tsv") lines"
  fi
  echo
}

run_test "TEST_B_local_tmp"     "/tmp/dtest_local"
run_test "TEST_A_nfs_fastscratch" "/fastscratch/${USER}/dtest_nfs"

echo "=========================================================="
echo "Interpretation:"
echo "  B(local)=0  A(nfs)!=0  -> NFS temp is the cause; fix = scratch=true on CHECKV"
echo "  B=0  A=0               -> rebuild fixed it; just resubmit the pipeline"
echo "  B!=0 A!=0              -> real DIAMOND error is now visible above; fix that"
echo "=========================================================="
