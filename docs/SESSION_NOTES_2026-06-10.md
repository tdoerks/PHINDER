# PHINDER Development Session - June 10, 2026

**Branches:** `main` (short-read fixes) + `feature/long-reads` (Phase 1 validated)
**Status:** Long-read Phase 1 assembler chain proven end-to-end on simulated data;
short-read pipeline pushed past its two long-standing walls (CheckV, Pharokka).

---

## 🎯 Session Goals

1. Diagnose why the latest 3-phage run (job 9429958) failed
2. Get the short-read pipeline past CheckV and Pharokka to a full end-to-end run
3. Stand up and test the long-read (Phase 1) path on real (simulated) data
4. Unify the SLURM notification email to the new `@ksu.edu` address

---

## 🐛 Short-read fixes (on `main`)

### 1. CheckV — concurrent DIAMOND failure (`cc5c5c0`)
**Symptom (job 9429958):** CheckV `contamination` (HMMER) passed, then `completeness`
died at `[3/8] Running DIAMOND blastp search... DIAMOND task failed. Program should be
rerun.`

**Ruled out:** the DB (the `ba109da` rebuild was already in place) and the filesystem —
`bin/diagnose_checkv_diamond.sh` reproduced the exact blastp by hand and it **passed on
both local `/tmp` and NFS `/fastscratch`** (B=0, A=0), 28,138 alignments.

**Root cause:** concurrency. Under the `local` executor all 3 CheckV tasks ran at once
(`executor > local (7)`), each launching a DIAMOND that sizes its memory block from
*total* node RAM (128 G), not its share. The simultaneous allocations collided
(job MaxRSS only ~50 G of 128 G — a transient spike the 30 s sacct sampler missed). A
single DIAMOND has the whole node and succeeds.

**Fix:** `maxForks = 1` on `CHECKV` and `DIAMOND_PROPHAGE` — one DIAMOND-backed task at a
time, matching the verified-good single-process condition.

### 2. Pharokka — dependency-check crash, NOT the contig count (`a0e8a13`)
**Symptom:** `PHAROKKA (SRR5131136)` exited 1. The actual error (swallowed by Nextflow's
retry line) was at **"Checking dependencies"**, *before the assembly is ever read*:
```
File "input_commands.py", line 386, in check_dependencies
    phanotate_major_version = int(phanotate_version.split(".")[0])
ValueError: invalid literal for int() with base 10: '/usr/local/lib/python3'
```
A numpy `UserWarning: smallest subnormal ... is zero` (hardware-triggered on the compute
node) polluted the string Pharokka parses for PHANOTATE's version.

**Key realization:** this is a **container bug, not a data problem** — it fails on *every*
sample (the crash precedes input). The contig filter would NOT have fixed it.

**Smoke-tested 1.7.x builds with a real `pharokka.py --fast` run on a Beocat compute node**
(the warning is hardware-triggered, so it had to be tested on-node) against the existing
`pharokka_db`:

| Build | Result |
|---|---|
| `1.9.1` / `1.8.2` | ❌ DB incompatible (want a newer DB) |
| `1.7.5--pyhdfd78af_2` | ❌ phanotate `ValueError`; +`PYTHONWARNINGS=ignore` gets past it, then MMseqs2 is **v18** (pharokka 1.7.5 needs **v13.45111**) |
| `1.7.5--pyhdfd78af_1` | ❌ phanotate `ValueError` |
| **`1.7.5--pyhdfd78af_0`** | ✅ phanotate v1.5.1 ok, **MMseqs2 v13.45111 ok**, tRNAscan/MinCED/ARAGORN/mash ok, full run exit 0 |
| `1.7.4--pyhdfd78af_0` | ✅ (also clean) |

**Fix:** pin `pharokka:1.7.5--pyhdfd78af_0` — the *original* 1.7.5 build (numpy-1.x, correct
MMseqs2). Same version → existing DB untouched. The 06-08 notes had chosen the *newest*
build `_2`, but that rebuild drifted to the wrong MMseqs2 + a numpy-2 stack.

### 3. MultiQC — missing-output failure with a title set (`d8c1b21`)
**Symptom:** `Error executing process > MULTIQC` even though `multiqc` exited 0 and logged
`MultiQC complete`.

**Root cause:** the module declared fixed output names (`multiqc_report.html` /
`multiqc_data`), but `multiqc_title='PHINDER Report'` makes MultiQC **prefix** the outputs
(`PHINDER-Report_multiqc_report.html` / `..._multiqc_report_data`). The declared names
don't exist → Nextflow reports a missing output.

**Fix:** match the title-prefixed names with globs (`*multiqc_report.html` / `*_data`).
Surfaced by the long-read smoke test; pre-emptively fixed on both branches.

---

## 🧬 Long-read Phase 1 — IMPLEMENTED → VALIDATED (on `feature/long-reads`)

Stack: **Flye → Medaka** (deliberately not Unicycler; see the rejection rationale in
`docs/LONG_READS_PLAN.md`). Brought from "implemented, untested" to "proven end-to-end".

### Container tag pre-flight (`c252b8d`, `1ed86bb`)
Verified all four long-read container tags against quay.io *before* running (the recurring
lesson):
- `nanoplot` / `badread` — original tags OK
- `flye:2.9.5--py39hdd1f253_0` → **MISSING**, corrected to `2.9.5--py39hdf45acc_0`
- `medaka:1.11.3--py39h05d5c5e_2` → **MISSING**, corrected to `..._0`
- `nanoplot:1.42.0` → pulled fine but **crashed at runtime** with the numpy 1.x/2.x ABI break
  (`_ARRAY_API not found`). Smoke-tested builds in-container; pinned **`1.47.0--pyhdfd78af_0`**
  (clean numpy stack). Same "pull ≠ works" lesson as Pharokka.

### Test data + run script (`cbd4964`)
`bin/run_phinder_longreads_beocat.sh` — idempotently builds Tier-1 data (Lambda ref J02459
→ Badread-simulated Nanopore reads → samplesheet) and runs the `long_reads` path. Defaults
to a smoke test (NanoPlot → Flye → Medaka, downstream skipped); `SKIP_DOWNSTREAM=false` for
the full run.

### Result
Badread: 187 reads / 2.43 Mb (~50×), N50 22.6 kb, 95% identity. Pipeline output:
```
>contig_1
total assembled: 48524 bp     (Lambda reference = 48,502 bp — 0.04% off)
```
**One clean contig.** NanoPlot ✅, Flye ✅, Medaka ✅ (no numpy gremlin). The MULTIQC bug
(above) was the only failure and is now fixed.

---

## ✉️ Email unification (`214011c`, `a3023ad`)
`tdoerks@vet.k-state.edu` → `tylerdoe@ksu.edu` (K-State unifying addresses). Updated the
3-phage script and added `END,FAIL` notifications to the new long-read script (it had none).

---

## ⚠️ Open item: SRR5131136 is dead data (not a pipeline bug)

The 280,540-contig assembly is host-DNA carryover, now quantified: **longest contig only
10,818 bp** for a ~40 kb T7 phage — there is *no* phage-sized contig at all. A length filter
can't recover a phage that isn't there.

- A `FILTER_CONTIGS` step (seqkit, `min_contig_length` param) before annotation is still
  worth adding — it helps CheckV/VIBRANT/DIAMOND performance and result sanity — but it is
  **not** what unblocked Pharokka (that was the container).
- Consider dropping SRR5131136 from the test set, or revisiting its input reads. Both
  SRR5131134/136 `proteins.faa` are ~27–34 MB, so the whole 3-phage set is contaminated,
  not clean isolates.

---

## 🚀 Next Steps
- [ ] Confirm short-read run reaches `PHINDER_SUMMARY` (resubmit once for the MULTIQC fix; `-resume` carries everything)
- [ ] Confirm long-read smoke test now exits 0 clean (MULTIQC fix)
- [ ] **Merge `main` → `feature/long-reads`** so the long-read branch inherits the CheckV `maxForks` + Pharokka `_0` fixes
- [ ] Full long-read run (`SKIP_DOWNSTREAM=false`) — validate CheckV/Pharokka/VIBRANT/Phanotate on the Flye/Medaka assembly
- [ ] Add the `FILTER_CONTIGS` step; decide test-set data quality (drop/replace SRR5131136)
- [ ] Phase 2 — hybrid mode (`modules/polypolish.nf`, `input_mode='hybrid'`)

---

## 🔑 Lessons reinforced
- **A successful `pull` ≠ a working image** — bit us twice today (NanoPlot 1.42, Pharokka `_2`).
  Runtime smoke tests (`--version`, *and a real run* for dep checks) catch what manifest checks can't.
- **Read the actual error before fixing.** The Pharokka failure *looked* like the 280 k-contig
  data problem; it was a dependency-check crash that had nothing to do with the assembly.
- **The newest container build is not always the right one.** Pharokka `_0` (original) beats `_2`
  (rebuild) because the rebuild drifted dependency versions.
- **DIAMOND sizes memory from total node RAM, not its cgroup share** — serialize DIAMOND-backed
  steps under the local executor.

**Last Updated:** 2026-06-10
**Maintained By:** Tyler Doerksen
