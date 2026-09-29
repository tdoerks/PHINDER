# PHINDER Test Phages

This document lists the test phages used for validating the PHINDER pipeline.

## Quick Test (3 phages)

**File:** `test_phages_sra.txt`
**Runtime:** ~2-4 hours
**Script:** `bin/run_phinder_full_test_beocat.sh`

| SRA Accession | Phage | Host | Genome Size | Lifestyle | Notes |
|--------------|-------|------|-------------|-----------|-------|
| SRR17327631 | Lambda | *E. coli* | 48.5 kb | Temperate | WGS, HiSeq 4000, ~15 Mb (~300x) |
| SRR19649190 | T4 | *E. coli* | 169 kb | Lytic | WGS, MiSeq, ~33 Mb (~200x) |
| ERR10819273 | T7 | *E. coli* | 40 kb | Lytic | WGS, NextSeq 500, ~10 Mb (~250x) |

> ⚠️ **Corrected 2026-09-29.** Until then this test used SRR5131134/5/6, which are **not phages** —
> they are RNA-Seq of scallop (*Azumapecten farreri*) foot tissue (~3 Gb each, study SRP018107).
> Every "lambda/T4/T7" result from those runs (incl. the `results_3phages_*` runs and the dsDNA
> rows of the SPAdes-mode comparison) came from a scallop transcriptome. The replacements above
> were verified by organism (NCBI taxid 10710 / 10665 / 10760) in ENA; they come from
> experimental-evolution/mutant studies, so expect a few SNPs vs RefSeq.

## Single Phage Test

**File:** `test_lambda_sra.txt`
**Runtime:** ~1 hour
**Script:** `bin/run_phinder_sra_beocat.sh`

| SRA Accession | Phage | Host | Genome Size | Lifestyle | Notes |
|--------------|-------|------|-------------|-----------|-------|
| SRR17327631 | Lambda | *E. coli* | 48.5 kb | Temperate | Quick validation (verified WGS) |

## Extended Test (20 phages)

> ⚠️ The 17 "diverse collection" accessions below are sequential run IDs with placeholder notes and have
> **not been verified** to be phage isolates. Prefer the stress test (below) for capability testing.

**File:** `test_phages_20_sra.txt`
**Runtime:** ~12-24 hours
**Script:** `bin/run_phinder_20phages_beocat.sh`

### Classic Model Phages (3)
| SRA Accession | Phage | Host | Genome Size | Lifestyle | Notes |
|--------------|-------|------|-------------|-----------|-------|
| SRR17327631 | Lambda | *E. coli* | 48.5 kb | Temperate | Integration/excision |
| SRR19649190 | T4 | *E. coli* | 169 kb | Lytic | Modified DNA bases |
| ERR10819273 | T7 | *E. coli* | 40 kb | Lytic | RNA polymerase |

### Diverse Collection (17)
| SRA Accession | Expected Type | Notes |
|--------------|---------------|-------|
| SRR8437269 | Phage isolate | Various hosts |
| SRR8437270 | Phage isolate | Genome diversity |
| SRR8437271 | Phage isolate | Size variation |
| SRR8437272 | Phage isolate | Lifestyle testing |
| SRR8437273 | Phage isolate | Quality checks |
| SRR8437274 | Phage isolate | Assembly testing |
| SRR8437275 | Phage isolate | Annotation testing |
| SRR11537895 | Recent isolate | Modern dataset |
| SRR11537896 | Recent isolate | 2020 sequences |
| SRR11537897 | Recent isolate | Updated tools |
| SRR11537898 | Recent isolate | Benchmark data |
| SRR13145901 | Latest isolate | 2021+ sequences |
| SRR13145902 | Latest isolate | Current standards |
| SRR13145903 | Latest isolate | Latest protocols |
| SRR13145904 | Latest isolate | Quality metrics |
| SRR13145905 | Latest isolate | Complete genomes |
| SRR13145906 | Latest isolate | Validation set |

## Stress Test (capability matrix, known answers)

**Branch:** `stress-test` · **Script:** `bin/run_phinder_stress_beocat.sh`

Simulated 2x150 reads from 19 verified RefSeq phage genomes plus E. coli K-12, and real reads already on
Beocat, run through all three input modes, then scored against the known answer (assembly recovery,
CheckV/geNomad detection, BacPhlip/VIBRANT lifestyle, per-module status, dashboard presence).

| Axis | Samples |
|------|---------|
| Genome size / architecture | phiX174, M13, MS2, Phi6 (3 segments), PM2, PRD1, T7, T5, T4, phiKZ (280 kb), phage G (498 kb) |
| Lifestyle | lambda, P22, Mu, P1 (temperate) vs the lytic phages above |
| Host / GC | phage K (Staph, low GC), L5 + D29 (Mycobacterium, high GC), crAss001 (Bacteroides) |
| Coverage | T7 at 5x, 20x, 100x, 1,000x, 10,000x |
| Host contamination | T7 + 10/50/90% E. coli reads; lambda + 50% |
| Mixed isolates | T7 + lambda, T4 + T7 |
| Negative controls | E. coli K-12 only, random sequence |
| Input modes | reads (all above + real reads), assembly (reference FASTAs), sra (verified lambda SRR17327631, T4 SRR19649190, T7 ERR10819273) |

```bash
python3 bin/stress_fetch_refs.py --outdir stress_data/refs          # downloads + verifies each accession
python3 bin/stress_simulate.py --refs stress_data/refs --outdir stress_data \
    --real-reads /path/to/samplesheet_spades_compare.csv          # optional real reads
sbatch bin/run_phinder_stress_beocat.sh                             # all 3 modes, then scoring
```

### Breadth tier — 1000 diverse phages (which *kinds* of phage break PHINDER)

`assets/stress_breadth_panel.tsv` (committed) lists 1000 complete RefSeq genomes with bacterial or archaeal
hosts, chosen for diversity rather than abundance: all 12 virus classes (incl. ssDNA, ssRNA, dsRNA, tailless
dsDNA), all 129 families, 230 host genera across 22 phyla, 131 archaeal viruses, 65 jumbo phages (>200 kb).
No host genus contributes more than ~30. Each genome is simulated clean at 100x, and the scorer adds
geNomad-vs-NCBI taxonomy agreement and a "where it breaks" breakdown by class, family, host and size.

```bash
python3 bin/stress_fetch_refs.py --panel assets/stress_breadth_panel.tsv --outdir stress_data/breadth_refs   # ~1 min, 70 MB
python3 bin/stress_simulate.py --breadth-panel assets/stress_breadth_panel.tsv --refs stress_data/breadth_refs --outdir stress_data   # ~9 min, ~3.5 GB
MODES=breadth sbatch bin/run_phinder_stress_beocat.sh     # ~12k SLURM tasks, 2-3 days; resubmit to resume
```

Rebuild the panel (new RefSeq release) from NCBI Datasets summaries, one query per virus class:
`datasets summary virus genome taxon <Class> --refseq --complete-only --as-json-lines > cat_<Class>.jsonl`, then
`python3 bin/stress_breadth_panel.py cat_*.jsonl --n 1000 --out assets/stress_breadth_panel.tsv`.
Segmented genomes (e.g. Cystoviridae) are one RefSeq record per segment, so a breadth sample is one segment.

Output: `stress_runs/stress_scorecard.html` + `.tsv` + `_breakdown.tsv`. Failed tasks are retried twice then recorded
(`conf/stress.config`), so one crash does not stop the run. PhageTerm calls on *simulated* reads are
not meaningful (read ends are uniform) — judge packaging only on the real-read samples.


### Quick Test (recommended first)
```bash
cd /fastscratch/tylerdoe/PHINDER
sbatch bin/run_phinder_sra_beocat.sh
```

### Full 3-Phage Test
```bash
cd /fastscratch/tylerdoe/PHINDER
sbatch bin/run_phinder_full_test_beocat.sh
```

### Extended 20-Phage Test
```bash
cd /fastscratch/tylerdoe/PHINDER
sbatch bin/run_phinder_20phages_beocat.sh
```

## Expected Results

Each phage will generate:
- **FastQC reports** - Read quality metrics
- **fastp reports** - Quality trimming results
- **Assembly** - Unicycler contigs (FASTA)
- **QUAST** - Assembly quality metrics
- **CheckV** - Genome completeness, contamination
- **Pharokka** - Gene annotations, functional predictions
- **VIBRANT** - Lifestyle prediction (lytic vs. temperate)
- **PHANOTATE** - Alternative gene calling
- **DIAMOND** - Prophage database hits

All results are summarized in:
- `multiqc/multiqc_report.html` - Aggregated QC metrics
- `summary/phinder_summary.html` - Per-sample phage analysis

## Download Results

```bash
# Download summary report
scp tylerdoe@beocat.cis.ksu.edu:/fastscratch/tylerdoe/PHINDER/results_*/summary/phinder_summary.html .

# Download MultiQC
scp tylerdoe@beocat.cis.ksu.edu:/fastscratch/tylerdoe/PHINDER/results_*/multiqc/multiqc_report.html .
```

## Troubleshooting

### SRA Download Failures
Some SRA accessions may fail to download due to:
- Network connectivity issues
- SRA database availability
- SLURM controller connectivity

**Solution:** Use the 20-phage test - even if some fail, others should succeed.

### SLURM Controller Errors
If you see `"Unable to contact slurm controller"`:
- This is a temporary Beocat infrastructure issue
- Wait 10-15 minutes and resubmit
- Check Beocat status: https://beocat.ksu.edu/status

### Memory Errors
If processes fail with OOM (out of memory):
- Increase `--mem` in SBATCH script
- Adjust process memory in `nextflow.config`

## Notes

- All phages are **publicly available** from NCBI SRA
- Accessions were selected for diversity in size, host, and lifestyle
- Some accessions may become deprecated - check SRA for alternatives
- Download times vary by file size and network conditions
