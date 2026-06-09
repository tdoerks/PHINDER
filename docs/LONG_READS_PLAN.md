# PHINDER Long-Read & Hybrid Assembly — Implementation Plan

**Branch:** `feature/long-reads` (refreshed onto `main` @ ba109da — 2026-06-09)
**Status:** Planning → Phase 1
**Scope:** Add Nanopore/PacBio (pure long-read) and Illumina+long-read (hybrid)
assembly **alongside** the existing short-read SPAdes path. Short-read isolate
mode on `main` is unchanged.

---

## 🎯 Decision: assembler stack = Flye → Medaka → Polypolish (NOT Unicycler)

The original plan (now superseded) built hybrid on **Unicycler hybrid mode**.
We are deliberately **not** doing that.

**Why Unicycler is rejected:** PHINDER started on Unicycler and abandoned it
(commit `028d43c9`, 2026-06-04) after **persistent SIGSEGV crashes (signal -11)
in Unicycler's internal SPAdes step**, during read mapping on Lambda phage data
on Beocat. Two days of memory/resource tuning (16→48 GB, CPU/time changes) did
not resolve it; the project switched to bare SPAdes for stability. Building
hybrid on Unicycler walks straight back into a proven-unstable tool on this
exact data and cluster. (This was a segfault in the assembly core — not a
container/version issue like the CheckV/Pharokka ones we fixed.)

**Chosen stack** (also what Unicycler's own author now recommends over Unicycler):

| Mode | Chain |
|---|---|
| Pure long-read (`input_mode=long_reads`) | **Flye** → **Medaka** (Nanopore polish) |
| Hybrid (`input_mode=hybrid`)             | **Flye** → **Medaka** → **Polypolish** (Illumina short-read polish) |

**Why this stack:**
- Flye is rock-solid and purpose-built for the 5–200 kb phage genome range;
  avoids the SIGSEGV entirely.
- Current best practice for phage/bacterial long-read & hybrid assembly.
- Modular — each step is a clean Nextflow process. Pure-long-read mode is just
  hybrid minus the final Polypolish step, so we build incrementally.
- SPAdes stays untouched on `main` for short-read isolates; these are separate
  input modes selected by `--input_mode`.

---

## 🧪 Test data strategy (two tiers)

**Tier 1 — development (build the plumbing, ground-truth, no data hunting):**
Simulate long reads from the Lambda/T4/T7 references we already use, pair with
the real Illumina reads already in the repo.
- Nanopore sim: **Badread** (`quay.io/biocontainers/badread`)
- PacBio sim: **pbsim3** (optional)
- Gives controlled hybrid + long-read inputs with a known answer to assemble
  against — ideal for wiring/validating each module before real data.

**Tier 2 — validation (Phase 3, real data):**
- **PRJEB56639** (ENA, *Briefings in Bioinformatics* 2024 viromics benchmark):
  has Illumina MiSeq (3 runs), Nanopore MinION (7), PacBio Sequel (3) of a mock
  phage community. NOTE: labeled "viral metagenome" — a mock *community*, not
  pure isolates, so it's a benchmark fit, not a clean isolate test.
- **Better:** SRA search for a phage *isolate* deposited with BOTH Illumina and
  Nanopore (common in recent phage genome papers). TODO: find accessions.

---

## 🏗️ Architecture

### Input modes (via `--input_mode`)
- `reads`   — Illumina short reads → SPAdes (EXISTING, unchanged)
- `assembly`— pre-assembled FASTA (EXISTING)
- `sra`     — download from SRA (EXISTING)
- `long_reads` — Nanopore/PacBio only → Flye + Medaka (NEW)
- `hybrid`  — Illumina + long reads → Flye + Medaka + Polypolish (NEW)

### Samplesheet formats
Long-read only:
```csv
sample,long_reads,platform
Lambda,lambda_ont.fastq.gz,nanopore
```
Hybrid:
```csv
sample,read1,read2,long_reads,platform
Lambda,lambda_R1.fq.gz,lambda_R2.fq.gz,lambda_ont.fq.gz,nanopore
```
`platform` ∈ {nanopore, pacbio} — drives Flye preset (`--nano-hq` / `--pacbio-hifi`)
and whether Medaka runs (Nanopore only).

---

## 🔧 New modules

| Module | Container | Role |
|---|---|---|
| `nanoplot.nf`  | `quay.io/biocontainers/nanoplot`  | long-read QC (feeds MultiQC) |
| `flye.nf`      | `quay.io/biocontainers/flye`      | long-read assembly (preset by platform) |
| `medaka.nf`    | `ontresearch/medaka` or biocontainer | Nanopore consensus polish |
| `polypolish.nf`| `quay.io/biocontainers/polypolish` | Illumina short-read polish (hybrid) |
| `badread.nf`   | `quay.io/biocontainers/badread`   | (dev only) simulate long reads for testing |

All pinned to explicit version tags (lesson from the Pharokka/CheckV container
saga — never use `:latest`, audit the registry tag exists before wiring).

---

## 📋 Phases

### Phase 1 — Pure long-read mode  ← START HERE
1. `nanoplot.nf` (long-read QC) + `flye.nf` + `medaka.nf`
2. Samplesheet parser: accept `long_reads` + `platform`; add `input_mode=long_reads`
3. Skip fastp/FASTQC (short-read QC) on this path; run NanoPlot instead
4. Feed polished assembly into the EXISTING downstream (CheckV, Pharokka,
   VIBRANT, DIAMOND, Phanotate, Summary) — they take a FASTA, so they're reused
5. Test with Badread-simulated reads from Lambda/T4/T7

### Phase 2 — Hybrid mode
1. `polypolish.nf` (align Illumina with bwa/minimap2, then polypolish)
2. `input_mode=hybrid`: Flye(long) → Medaka → Polypolish(short)
3. Samplesheet with both short + long columns
4. Test with simulated long + real Illumina

### Phase 3 — Validation & benchmarking
1. Real data (PRJEB56639 and/or a real phage isolate w/ paired platforms)
2. Run all modes; compare QUAST metrics, CheckV completeness, Pharokka gene calls
3. Document recommended mode per data type in README

---

## 🔑 Constraints / lessons carried in
- **Pin every container tag**; verify it resolves before wiring (Pharokka 404 /
  numpy-ABI saga, CheckV DIAMOND build skew).
- **Reproduce assembly crashes in isolation** before tuning resources — Flye is
  the choice precisely to avoid Unicycler's SIGSEGV.
- Downstream phage modules already accept a FASTA, so long-read/hybrid only
  needs to produce a clean assembly and hand it off — minimal downstream change.

---

**Last Updated:** 2026-06-09
**Maintainer:** Tyler Doerksen
**Next:** implement Phase 1 — `flye.nf` + `nanoplot.nf` + `medaka.nf`, simulated test data
