# PHINDER Development Session - June 8, 2026

**Branch:** `session-notes` (fixes landed on `main`)
**Status:** Container + database blockers cleared; full 3-phage run resubmitted (job 9425907)

---

## 🎯 Session Goals

1. Diagnose why the latest 3-phage run (job 9362601) failed
2. Stop the recurring "fix one container, hit the next" cycle by auditing everything at once
3. Confirm whether the Unicycler→SPAdes switch actually resolved the assembly crashes
4. Get the pipeline running end-to-end through PHINDER_SUMMARY

---

## 🐛 Major Issues Found & Fixed

### Issue 1: PHAROKKA container removed from Docker Hub (the failure that killed job 9362601)
**Error:**
```
Error executing process > 'PHINDER_PIPELINE:PHAROKKA (1)'
  Failed to pull singularity image
    docker://gbouras13/pharokka:1.7.0
    FATAL: ... failed to get checksum ... reading manifest 1.7.0 in
    docker.io/gbouras13/pharokka: requested access to the resource is denied
```
**Root cause:** `gbouras13/pharokka` no longer resolves on Docker Hub (pull → 401; repo not even in search). The `access denied` on the manifest is the "image is gone/private" signal, not an auth problem.

**Fix:** Switched to BioContainers — same pattern as the earlier VIBRANT (`4796c1e`) and DIAMOND (`b11dbca`) fixes.

---

### Issue 2: MULTIQC container never existed (latent — caught before it bit)
Found by auditing **all** container references against their registries, not just the one that failed.

**Finding:** `ewels/multiqc:1.25.1` → manifest 404. `ewels/multiqc` was abandoned at `v1.19` (last pushed 2023); the project moved to the `multiqc/multiqc` namespace. So `1.25.1` never existed under `ewels/`. MULTIQC runs late in the DAG, so this would have failed the *next* run the moment PHAROKKA was fixed — a second wasted submit/wait/fail cycle, avoided.

**Fix:** `docker://quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0` (same pinned version, verified pullable).

**Container audit result (all refs):** the two above were broken; everything else — `staphb/{quast,fastqc,fastp}` and all 8 `quay.io/biocontainers/*` — verified OK. (Earlier worry that staphb images were unreliable was unfounded.)

---

### Issue 3: PHAROKKA database not installed + misleading config/docs
**Finding:** `pharokka_db = null` with a comment claiming "Auto-downloaded if not specified." Pharokka does **not** download its PHROGs DB at runtime — `install_databases.py` is a separate one-time step, and the BioContainers image does not bundle the DB. So the annotation run would have failed "database not found" right after the container was fixed.

**Fix:**
- Defaulted `pharokka_db = /fastscratch/tylerdoe/databases/pharokka_db`
- Rewrote `docs/DATABASE_SETUP.md` Pharokka section + storage table to show the real one-time install step (runnable via the pipeline's own container)
- Installed the DB on Beocat and verified all components present (PHROGs `all_phrogs.h3m`, VFDB, CARD, INPHARED Mash `1Aug2023_genomes.fa.msh`, `phrog_annot_v4.tsv`) — 689 MB total

---

### Issue 4: BioContainers pharokka `1.7.0--pyhdfd78af_0` is itself broken (numpy/pandas ABI)
Surfaced when the DB install crashed:
```
ValueError: numpy.dtype size changed, may indicate binary incompatibility.
Expected 96 from C header, got 88 from PyObject
```
**Root cause:** That specific build ships an incompatible numpy/pandas pairing (the classic numpy 1.x vs 2.x ABI break). Critically, `pharokka.py` imports the same `post_processing`→`pandas` chain, so the PHAROKKA *step* would crash identically — not just the DB installer.

**Smoke-tested candidate builds on Beocat** (`import pandas, numpy` + `pharokka.py --version`):

| Build | pandas / numpy | Result |
|---|---|---|
| `1.7.0--pyhdfd78af_0` | — | ❌ ValueError |
| `1.7.5--pyhdfd78af_2` | 2.3.2 / 2.3.2 | ✅ **chosen** |
| `1.8.2--pyhdfd78af_0` | 2.3.3 / 2.3.4 | ✅ |
| `1.9.1--pyhdfd78af_1` | 3.0.0 / 2.4.2 | ✅ (pandas 3.x too new) |

**Fix:** Pinned PHAROKKA to `1.7.5--pyhdfd78af_2` — closest working build to the original 1.7.0 pin, least behavioral drift.

---

## ✅ Assembly Question Resolved (Unicycler vs SPAdes)

The recurring "main issue" was historically **assembly crashes (Unicycler)**. Confirmed this is **fixed and is not the current blocker**:

- Unicycler is fully removed — no module, zero references; SPAdes is the sole assembler (since `028d43c`)
- SPAdes config is robust: 8 CPU, 32 GB × attempt (auto-escalates to 64 GB on retry), 12 h × attempt, `maxRetries=2`, `--only-assembler` (skips the memory-hungry error-correction stage that was the likely crash source)
- Job 9362601's log showed SPAdes `1 of 3, cached ✔` then died at the PHAROKKA *pull* — assembly succeeded; the container was the blocker

**The failure mode migrated:** assembly (Unicycler) → fixed via SPAdes → downstream container availability → now cleared.

---

## ⚠️ Open Item: assembly *quality* (not a crash)

`SRR5131136_assembly.fasta` is **82 MB** — a purified phage genome should be ~30–200 kb. That much sequence strongly suggests host-DNA carryover / contaminated input. Not a pipeline failure, but worth investigating:
```bash
grep -c '^>' .../assemblies/SRR5131136_assembly.fasta   # contig count
```
If this is thousands of contigs, consider host-read filtering before assembly. Revisit after a clean end-to-end run.

---

## 📦 Commits (on `main`)

- `d854793` — Fix broken container images (pharokka, multiqc) + correct pharokka DB docs
- `c7a2cad` — Pin Pharokka to 1.7.5--pyhdfd78af_2 (fix numpy/pandas ABI break)

---

## 🚀 Next Steps

### Immediate
- [ ] Watch job **9425907** (3-phage) clear the first-ever real `PHAROKKA` execution
- [ ] Confirm `VIBRANT → DIAMOND_PROPHAGE → PHANOTATE → MULTIQC → PHINDER_SUMMARY` complete
- [ ] Check `phinder_report.html` renders with annotation data

### Follow-up
- [ ] Investigate the 82 MB assembly (contig count, host contamination, read-filtering)
- [ ] Runtime-smoke-test VIBRANT/DIAMOND pandas imports (passed pull check, not import check)
- [ ] Merge `main` into `session-notes` so this branch reflects current code state

---

## 🔑 Lessons / Patterns

- **Audit the whole class of problem at once.** Container references were failing one-at-a-time across multiple days; probing every registry in one pass caught the MULTIQC landmine before it cost another run.
- **A successful `pull` ≠ a working image.** The pharokka 1.7.0 biocontainer pulled fine but crashed on `import pandas`. Runtime smoke tests (`import`, `--version`) catch what manifest checks can't.
- **Trust but verify vendor docs.** The "Pharokka auto-downloads its DB" claim (in our own docs) was wrong and would have caused the next failure.

---

**Session End Status:** All known blockers cleared (assembly, 2 containers, pharokka ABI, pharokka DB; CheckV + prophage DBs verified present). 3-phage run resubmitted as job 9425907 — awaiting first clean end-to-end completion.

**Last Updated:** 2026-06-08
**Maintained By:** Tyler Doerksen
