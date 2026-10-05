#!/bin/bash
#SBATCH --job-name=phinder_stress
#SBATCH --partition=batch.q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=72:00:00
#SBATCH --output=phinder_stress_%j.log
#SBATCH --error=phinder_stress_%j.err
#SBATCH --mail-type=END,FAIL

#==============================================================================
# PHINDER stress test — Beocat
#==============================================================================
# Runs PHINDER in all three input modes against simulated + real phage data with
# known answers, then scores every sample x module.
#
# One-time setup (from the repo root, stress-test branch):
#   python3 bin/stress_fetch_refs.py --outdir stress_data/refs
#   python3 bin/stress_simulate.py --refs stress_data/refs --outdir stress_data \
#       --real-reads /fastscratch/tylerdoe/PHINDER-spades-test/samplesheets/samplesheet_spades_compare.csv
#
# Then:  sbatch bin/run_phinder_stress_beocat.sh
#
# Breadth tier (1000 diverse RefSeq phages, clean 100x) — one-time setup, then submit separately:
#   python3 bin/stress_fetch_refs.py --panel assets/stress_breadth_panel.tsv --outdir stress_data/breadth_refs
#   python3 bin/stress_simulate.py --breadth-panel assets/stress_breadth_panel.tsv \
#       --refs stress_data/breadth_refs --outdir stress_data
#   MODES=breadth sbatch bin/run_phinder_stress_beocat.sh
# ~12k SLURM tasks; the slurm profile submits 10/min, so expect 2-3 days. If the head job hits its
# time limit, resubmit the same command — -resume picks up where it stopped.
#
# Each mode runs from its own launch dir so -resume tracks each run separately.
# Failed tasks are retried twice, then ignored and recorded (conf/stress.config).
#==============================================================================

set -uo pipefail   # no -e: a failed mode must not stop the next mode or the scoring

REPO="$(pwd)"
DATA="${REPO}/stress_data"
RUNS="${REPO}/stress_runs"
MODES="${MODES:-reads assembly sra}"   # override: MODES="assembly" sbatch ...

module load Nextflow/24.04.2 2>/dev/null || module load Nextflow || { echo "ERROR: no Nextflow module"; exit 1; }
# Cap the Nextflow JVM: by default its heap is sized to the whole node, not this job's memory,
# and SLURM kills the head job (OUT_OF_MEMORY) — job 11552233 died this way after 13 h.
export NXF_OPTS="${NXF_OPTS:--Xms1g -Xmx6g}"

declare -A INPUT=(
    [reads]="${DATA}/samplesheet_stress_reads.csv"
    [assembly]="${DATA}/samplesheet_stress_assembly.csv"
    [sra]="${DATA}/stress_sra.txt"
    [breadth]="${DATA}/samplesheet_stress_breadth.csv"
)
declare -A NF_MODE=([reads]=reads [assembly]=assembly [sra]=sra [breadth]=reads)

# Pull any missing containers first, with retries (Nextflow aborts on a failed pull)
bash "${REPO}/bin/prefetch_containers.sh" || echo "WARNING: some container pulls failed — Nextflow will retry them"

for mode in $MODES; do
    [ -f "${INPUT[$mode]:-/nonexistent}" ] || { echo "ERROR: input for '${mode}' missing — run the one-time setup (see header)"; exit 1; }
done

echo "PHINDER stress test — job ${SLURM_JOB_ID:-local} on ${SLURM_NODELIST:-$(hostname)} — $(date)"
echo "Commit: $(git -C "$REPO" rev-parse --short HEAD) ($(git -C "$REPO" branch --show-current))"

for mode in $MODES; do
    outdir="${RUNS}/${mode}/results"
    mkdir -p "${RUNS}/${mode}"
    echo ""
    echo "======== ${mode} mode — $(date) ========"
    ( cd "${RUNS}/${mode}" && nextflow run "${REPO}/main.nf" \
        -profile slurm,beocat \
        -c "${REPO}/conf/stress.config" \
        --input "${INPUT[$mode]}" \
        --input_mode "${NF_MODE[$mode]}" \
        --outdir "${outdir}" \
        -resume \
        -with-trace "${outdir}/pipeline_trace.txt" \
        -with-report "${outdir}/phinder_report.html" )
    echo "${mode} mode exit code: $?"
done

echo ""
echo "======== Scoring — $(date) ========"
# Scores every run that exists so far (stress and/or breadth); missing runs are skipped.
python3 "${REPO}/bin/score_stress_test.py" \
    --truth "${DATA}/stress_truth.tsv" --truth "${DATA}/stress_truth_breadth.tsv" \
    --refs "${DATA}/refs" --refs "${DATA}/breadth_refs" \
    --run reads="${RUNS}/reads/results" \
    --run assembly="${RUNS}/assembly/results" \
    --run sra="${RUNS}/sra/results" \
    --run breadth="${RUNS}/breadth/results" \
    --out "${RUNS}/stress_scorecard"

echo "Scorecard: ${RUNS}/stress_scorecard.html  |  ${RUNS}/stress_scorecard.tsv"
