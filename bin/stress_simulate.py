#!/usr/bin/env python3
"""
Simulate Illumina paired-end reads for the PHINDER stress test and write the samplesheets
plus a ground-truth table the scorer compares against.

Pure Python (no ART/InSilicoSeq dependency): 2x150 bp, insert ~N(350, 50), position-dependent
substitution errors (~0.1% at read start rising to ~1% at the end), fixed seed per sample.
Circular references wrap around the origin; multi-record references (Phi6 segments,
E. coli) are sampled proportionally to record length.

Usage:
    python3 bin/stress_simulate.py --refs stress_data/refs --outdir stress_data \
        [--real-reads samplesheets/samplesheet_spades_compare.csv]
    # breadth tier (1000-genome panel, clean 100x each):
    python3 bin/stress_simulate.py --breadth-panel assets/stress_breadth_panel.tsv \
        --refs stress_data/breadth_refs --outdir stress_data

Outputs (in --outdir):
    reads/<sample>_R1.fastq.gz, reads/<sample>_R2.fastq.gz
    samplesheet_stress_reads.csv       sample,read1,read2           (--input_mode reads)
    samplesheet_stress_assembly.csv    sample,assembly              (--input_mode assembly)
    stress_sra.txt                     one SRR per line             (--input_mode sra)
    stress_truth.tsv                   expected results per sample
  breadth tier instead writes samplesheet_stress_breadth.csv + stress_truth_breadth.tsv
"""
import argparse
import csv
import gzip
import math
import random
import sys
from pathlib import Path

READ_LEN = 150
INSERT_MEAN, INSERT_SD = 350, 50
ERR_START, ERR_END = 0.001, 0.01
HOST = "Ecoli_K12"

# (sample_id, category, [(ref_id, coverage)], host_fraction, notes)
# host_fraction = fraction of READ PAIRS drawn from the E. coli host on top of the phage reads.
SIM_DESIGN = [
    # Genome size / architecture — everything at 100x
    ("sim_phiX174", "architecture", [("phiX174", 100)], 0, "5.4 kb ssDNA, circular"),
    ("sim_M13",     "architecture", [("M13", 100)], 0, "filamentous ssDNA (Inoviridae)"),
    ("sim_MS2",     "architecture", [("MS2", 100)], 0, "3.6 kb ssRNA (as cDNA)"),
    ("sim_Phi6",    "architecture", [("Phi6", 100)], 0, "segmented dsRNA — 3 contigs expected"),
    ("sim_PM2",     "architecture", [("PM2", 100)], 0, "tailless, lipid, circular dsDNA"),
    ("sim_PRD1",    "architecture", [("PRD1", 100)], 0, "tailless, lipid, linear dsDNA"),
    ("sim_T7",      "architecture", [("T7", 100)], 0, "40 kb reference point"),
    ("sim_T5",      "architecture", [("T5", 100)], 0, "122 kb"),
    ("sim_T4",      "architecture", [("T4", 100)], 0, "169 kb"),
    ("sim_phiKZ",   "architecture", [("phiKZ", 100)], 0, "280 kb jumbo"),
    ("sim_phageG",  "architecture", [("phageG", 100)], 0, "498 kb jumbo — largest known phage class"),
    # Lifestyle
    ("sim_lambda",  "lifestyle", [("lambda", 100)], 0, "temperate, cos"),
    ("sim_P22",     "lifestyle", [("P22", 100)], 0, "temperate, pac"),
    ("sim_Mu",      "lifestyle", [("Mu", 100)], 0, "temperate, transposable"),
    ("sim_P1",      "lifestyle", [("P1", 100)], 0, "temperate, plasmid-like prophage"),
    # Host range / GC
    ("sim_phageK",  "host_gc", [("phageK", 100)], 0, "Staphylococcus, ~30% GC"),
    ("sim_L5",      "host_gc", [("L5", 100)], 0, "Mycobacterium, ~63% GC, temperate"),
    ("sim_D29",     "host_gc", [("D29", 100)], 0, "Mycobacterium, ~63% GC, lytic"),
    ("sim_crAss001", "host_gc", [("crAss001", 100)], 0, "Bacteroides crAss-like"),
    # Coverage (T7 100x above is the reference point)
    ("sim_T7_cov5",     "coverage", [("T7", 5)], 0, "very low coverage"),
    ("sim_T7_cov20",    "coverage", [("T7", 20)], 0, "low coverage"),
    ("sim_T7_cov1000",  "coverage", [("T7", 1000)], 0, "high coverage"),
    ("sim_T7_cov10000", "coverage", [("T7", 10000)], 0, "extreme coverage (phiX-calibration-like)"),
    # Host contamination (the SA_4 case)
    ("sim_T7_host10",     "contamination", [("T7", 100)], 0.10, "10% of pairs from E. coli"),
    ("sim_T7_host50",     "contamination", [("T7", 100)], 0.50, "50% of pairs from E. coli"),
    ("sim_T7_host90",     "contamination", [("T7", 100)], 0.90, "90% of pairs from E. coli"),
    ("sim_lambda_host50", "contamination", [("lambda", 100)], 0.50, "temperate phage + host"),
    # Mixed isolates
    ("sim_mix_T7_lambda", "mixed", [("T7", 100), ("lambda", 100)], 0, "two phages, similar size"),
    ("sim_mix_T4_T7",     "mixed", [("T4", 100), ("T7", 100)], 0, "two phages, 4x size difference"),
    # Negative controls
    ("sim_Ecoli_only", "negative", [], 1.0, "E. coli K-12 only at 30x — must NOT be called a phage"),
    ("sim_random",     "negative", [], 0, "random sequence — must NOT be called a phage"),
]
ECOLI_ONLY_COV = 30
RANDOM_PAIRS = 20000

# Assembly-mode run: the reference FASTAs themselves (tests the post-assembly modules in isolation)
ASSEMBLY_DESIGN = [
    ("asm_T4", "T4"), ("asm_phiKZ", "phiKZ"), ("asm_phageG", "phageG"), ("asm_Phi6", "Phi6"),
    ("asm_MS2", "MS2"), ("asm_M13", "M13"), ("asm_lambda", "lambda"), ("asm_Ecoli_K12", HOST),
]

# SRA-mode run: real phage data via DOWNLOAD_SRA (identity checked by the scorer's k-mer recovery).
# Verified by ENA organism (2026-09-29), Illumina paired WGS:
SRA_DESIGN = [
    ("SRR17327631", "lambda"),   # taxid 10710, HiSeq 4000, ~15 Mb (~300x)
    ("SRR19649190", "T4"),       # taxid 10665, MiSeq, ~33 Mb (~200x)
    ("ERR10819273", "T7"),       # taxid 10760, NextSeq 500, ~10 Mb (~250x)
]

# Real reads already on Beocat (SPAdes-compare samplesheet); sample -> claimed reference.
# SRR5131134/5/6 were long documented as lambda/T4/T7 but are scallop (Azumapecten farreri)
# RNA-Seq (ENA, checked 2026-09-29) — not phages, so they are excluded from the stress test.
EXCLUDE_REAL = {"SRR5131134", "SRR5131135", "SRR5131136"}
REAL_CLAIMS = {
    "phiX174": ("phiX174", "real reads SRR001665"),
    "Ecoli_K12": (None, "real reads SRR001666 — negative control"),
    "MS2": ("MS2", "real reads SRR31435157"),
    "Phi6": ("Phi6", "real reads SRR30985651"),
}

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    recs, seq = [], []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if seq:
                    recs.append("".join(seq))
                seq = []
            else:
                seq.append(line.strip().upper())
    if seq:
        recs.append("".join(seq))
    return recs


def load_manifest(refs_dir):
    manifest = {}
    with open(Path(refs_dir) / "refs_manifest.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            row["records"] = read_fasta(Path(refs_dir) / f"{row['ref_id']}.fasta")
            manifest[row["ref_id"]] = row
    return manifest


class Simulator:
    def __init__(self, seed):
        self.rng = random.Random(seed)
        # Per-position error probability and matching Phred quality string
        self.err = [ERR_START + (ERR_END - ERR_START) * i / (READ_LEN - 1) for i in range(READ_LEN)]
        self.qual = "".join(chr(33 + min(41, int(-10 * math.log10(p)))) for p in self.err)
        self.mean_errs = sum(self.err)

    def _poisson(self, lam):
        L, k, p = math.exp(-lam), 0, 1.0
        while True:
            p *= self.rng.random()
            if p <= L:
                return k
            k += 1

    def _mutate(self, read):
        k = self._poisson(self.mean_errs)
        if not k:
            return read
        r = list(read)
        for pos in self.rng.choices(range(READ_LEN), weights=self.err, k=k):
            r[pos] = self.rng.choice([b for b in "ACGT" if b != r[pos]])
        return "".join(r)

    def pairs_from(self, records, circular, n_pairs):
        rng = self.rng
        lengths = [len(s) for s in records]
        for _ in range(n_pairs):
            s = records[rng.choices(range(len(records)), weights=lengths)[0]] if len(records) > 1 else records[0]
            L = len(s)
            ins = max(READ_LEN, min(L, int(rng.gauss(INSERT_MEAN, INSERT_SD))))
            if circular:
                start = rng.randrange(L)
                frag = s[start:start + ins] if start + ins <= L else s[start:] + s[:start + ins - L]
            else:
                start = rng.randrange(L - ins + 1)
                frag = s[start:start + ins]
            if rng.random() < 0.5:
                frag = revcomp(frag)
            r1 = frag[:READ_LEN].ljust(READ_LEN, "N")
            r2 = revcomp(frag)[:READ_LEN].ljust(READ_LEN, "N")
            yield self._mutate(r1), self._mutate(r2)

    def random_pairs(self, n_pairs):
        for _ in range(n_pairs):
            yield ("".join(self.rng.choices("ACGT", k=READ_LEN)),
                   "".join(self.rng.choices("ACGT", k=READ_LEN)))


def n_pairs_for(total_len, cov):
    return max(1, round(cov * total_len / (2 * READ_LEN)))


def write_sample(sample, sources, sim, reads_dir):
    """sources: list of generators yielding (r1, r2), streamed to disk (read order is irrelevant
    to fastp/SPAdes, and 10,000x coverage would not fit in memory as a list)."""
    r1p, r2p = reads_dir / f"{sample}_R1.fastq.gz", reads_dir / f"{sample}_R2.fastq.gz"
    i = 0
    with gzip.open(r1p, "wt", compresslevel=3) as f1, gzip.open(r2p, "wt", compresslevel=3) as f2:
        for gen in sources:
            for a, b in gen:
                i += 1
                f1.write(f"@{sample}.{i}/1\n{a}\n+\n{sim.qual}\n")
                f2.write(f"@{sample}.{i}/2\n{b}\n+\n{sim.qual}\n")
    return r1p, r2p, i


TRUTH_COLS = ["sample", "mode", "category", "refs", "host_fraction", "expect_phage", "expect_lifestyle",
              "genome_type", "notes", "class", "family", "host_genus", "host_domain", "size_bin", "length"]


def simulate_breadth(args, refs, out, reads_dir):
    """One clean sample per panel genome at --breadth-cov; truth carries NCBI taxonomy + host."""
    with open(args.breadth_panel) as f:
        panel = list(csv.DictReader(f, delimiter="\t"))
    truth, rows, skipped = [], [], []
    for i, p in enumerate(panel):
        s = p["sample"]
        if s not in refs:
            skipped.append(s)
            continue
        r1p, r2p = reads_dir / f"{s}_R1.fastq.gz", reads_dir / f"{s}_R2.fastq.gz"
        if not (r1p.exists() and r2p.exists()):   # resumable: 1000 samples take a while
            ref = refs[s]
            sim = Simulator(args.seed + 100_000 + i)
            gen = sim.pairs_from(ref["records"], ref["topology"] == "circular",
                                 n_pairs_for(int(ref["length"]), args.breadth_cov))
            write_sample(s, [gen], sim, reads_dir)
        rows.append((s, r1p, r2p))
        truth.append(dict(sample=s, mode="breadth", category="breadth", refs=f"{s}:{args.breadth_cov}",
                          host_fraction=0, expect_phage="yes", expect_lifestyle="NA",
                          genome_type=f"{p['class']}/{p['family']}",
                          notes=f"{p['virus_name']} ({p['accession']}); selected for {p['selected_for']}",
                          **{k: p[k] for k in ("class", "family", "host_genus", "host_domain", "size_bin", "length")}))
        if (i + 1) % 100 == 0:
            print(f"  {i + 1}/{len(panel)}")
    with open(out / "samplesheet_stress_breadth.csv", "w") as f:
        f.write("sample,read1,read2\n")
        for s, r1, r2 in rows:
            f.write(f"{s},{r1},{r2}\n")
    with open(out / "stress_truth_breadth.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=TRUTH_COLS, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(truth)
    print(f"\n{len(rows)} breadth samples -> {out}/samplesheet_stress_breadth.csv")
    if skipped:
        print(f"SKIPPED {len(skipped)} (genome not downloaded — rerun stress_fetch_refs.py --panel): "
              f"{', '.join(skipped[:10])}{' ...' if len(skipped) > 10 else ''}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refs", default="stress_data/refs")
    ap.add_argument("--outdir", default="stress_data")
    ap.add_argument("--real-reads", help="CSV with sample,fastq_1,fastq_2 (e.g. the SPAdes-compare samplesheet)")
    ap.add_argument("--only", nargs="*", help="Simulate only these sample_ids")
    ap.add_argument("--seed", type=int, default=20260928)
    ap.add_argument("--truth-only", action="store_true",
                    help="Rewrite samplesheets + truth table but not the reads (existing reads are "
                         "reused, so an in-progress -resume run is not invalidated)")
    ap.add_argument("--breadth-panel", help="Panel TSV from stress_breadth_panel.py — simulate the breadth tier instead")
    ap.add_argument("--breadth-cov", type=int, default=100)
    args = ap.parse_args()

    out = Path(args.outdir).resolve()
    reads_dir = out / "reads"
    reads_dir.mkdir(parents=True, exist_ok=True)
    refs = load_manifest(args.refs)
    refs_dir = Path(args.refs).resolve()
    if args.breadth_panel:
        return simulate_breadth(args, refs, out, reads_dir)

    truth, reads_rows, skipped = [], [], []

    def truth_row(sample, mode, category, comps, host_frac, notes, expect_phage=None):
        lifes = {refs[r]["lifestyle"] for r, _ in comps if r in refs}
        life = lifes.pop() if len(lifes) == 1 else ("mixed" if lifes else "NA")
        gtype = "; ".join(refs[r]["genome_type"] for r, _ in comps if r in refs) or "none"
        if expect_phage is None:
            expect_phage = "yes" if comps else "no"
        truth.append(dict(sample=sample, mode=mode, category=category,
                          refs=";".join(f"{r}:{c}" for r, c in comps), host_fraction=host_frac,
                          expect_phage=expect_phage, expect_lifestyle=life, genome_type=gtype, notes=notes))

    # ── Simulated reads ───────────────────────────────────────────────────────
    for i, (sample, category, comps, host_frac, notes) in enumerate(SIM_DESIGN):
        if args.only and sample not in args.only:
            continue
        needed = [r for r, _ in comps] + ([HOST] if host_frac else [])
        missing = [r for r in needed if r not in refs]
        if missing:
            skipped.append(f"{sample} (missing ref: {', '.join(missing)})")
            continue
        sim = Simulator(args.seed + i)
        sources, phage_pairs = [], 0
        for ref_id, cov in comps:
            ref = refs[ref_id]
            n = n_pairs_for(int(ref["length"]), cov)
            phage_pairs += n
            sources.append(sim.pairs_from(ref["records"], ref["topology"] == "circular", n))
        host = refs.get(HOST)
        if host_frac and comps:
            n_host = round(phage_pairs * host_frac / (1 - host_frac))
            sources.append(sim.pairs_from(host["records"], True, n_host))
        elif host_frac == 1.0:
            sources.append(sim.pairs_from(host["records"], True, n_pairs_for(int(host["length"]), ECOLI_ONLY_COV)))
        if sample == "sim_random":
            sources.append(sim.random_pairs(RANDOM_PAIRS))
        if args.truth_only and (reads_dir / f"{sample}_R1.fastq.gz").exists():
            r1, r2 = reads_dir / f"{sample}_R1.fastq.gz", reads_dir / f"{sample}_R2.fastq.gz"
            print(f"  {sample:<20} (reads kept)")
        else:
            r1, r2, n = write_sample(sample, sources, sim, reads_dir)
            print(f"  {sample:<20} {n:>9,} pairs  ({category})")
        reads_rows.append((sample, r1, r2))
        truth_row(sample, "reads", category, comps, host_frac, notes)

    # ── Real reads already on disk ────────────────────────────────────────────
    if args.real_reads:
        with open(args.real_reads) as f:
            for row in csv.DictReader(f):
                s = row["sample"]
                if s in EXCLUDE_REAL:
                    skipped.append(f"real_{s} (excluded: scallop RNA-Seq, not a phage)")
                    continue
                r1, r2 = Path(row["fastq_1"]), Path(row["fastq_2"])
                if not (r1.exists() and r2.exists()):
                    skipped.append(f"real_{s} (reads not found: {r1})")
                    continue
                claimed, note = REAL_CLAIMS.get(s, (None, "real reads; no reference claim"))
                comps = [(claimed, "real")] if claimed else []
                reads_rows.append((f"real_{s}", r1, r2))
                truth_row(f"real_{s}", "reads", "real_reads", comps, 0, note,
                          expect_phage="yes" if claimed else "no")

    with open(out / "samplesheet_stress_reads.csv", "w") as f:
        f.write("sample,read1,read2\n")
        for s, r1, r2 in reads_rows:
            f.write(f"{s},{r1},{r2}\n")

    # ── Assembly mode ─────────────────────────────────────────────────────────
    with open(out / "samplesheet_stress_assembly.csv", "w") as f:
        f.write("sample,assembly\n")
        for sample, ref_id in ASSEMBLY_DESIGN:
            if ref_id not in refs:
                skipped.append(f"{sample} (missing ref: {ref_id})")
                continue
            f.write(f"{sample},{refs_dir / (ref_id + '.fasta')}\n")
            comps = [] if ref_id == HOST else [(ref_id, "ref")]
            truth_row(sample, "assembly", "assembly_mode", comps, 0, "reference FASTA as input")

    # ── SRA mode ──────────────────────────────────────────────────────────────
    with open(out / "stress_sra.txt", "w") as f:
        for srr, claimed in SRA_DESIGN:
            f.write(srr + "\n")
            truth_row(srr, "sra", "sra_mode", [(claimed, "real")], 0, f"DOWNLOAD_SRA path; verified {claimed} WGS (ENA organism check)")

    with open(out / "stress_truth.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=TRUTH_COLS, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(truth)

    print(f"\n{len(reads_rows)} reads-mode samples, {len(truth)} truth rows -> {out}")
    if skipped:
        print("SKIPPED:\n  " + "\n  ".join(skipped))


if __name__ == "__main__":
    main()
