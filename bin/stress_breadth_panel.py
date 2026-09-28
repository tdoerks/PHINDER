#!/usr/bin/env python3
"""
Build the PHINDER breadth-test panel: N complete RefSeq phage/archaeal-virus genomes chosen
to cover as much diversity as possible, not just the most-sequenced hosts.

Input is NCBI Datasets virus summaries (JSON Lines), e.g.:
    datasets summary virus genome taxon Caudoviricetes --refseq --complete-only --as-json-lines > cat.jsonl

Selection (deterministic for a given seed), in priority order:
  1. every non-Caudoviricetes class (non-tailed phages: ssDNA, ssRNA, dsRNA, tailless dsDNA)
  2. every family
  3. every host genus
  4. jumbo phages (>200 kb) and archaeal viruses, up to a quota each
  5. fill: round-robin across host genera, preferring unseen family x size-bin combinations
The output TSV records why each genome was picked, and is meant to be committed so the
test is reproducible without re-querying NCBI.

Usage:
    python3 bin/stress_breadth_panel.py cat_*.jsonl --n 1000 --out assets/stress_breadth_panel.tsv
"""
import argparse
import json
import random
from collections import defaultdict

SIZE_BINS = [(0, 10_000, "<10kb"), (10_000, 30_000, "10-30kb"), (30_000, 60_000, "30-60kb"),
             (60_000, 100_000, "60-100kb"), (100_000, 200_000, "100-200kb"), (200_000, 10**9, ">200kb")]
JUMBO_QUOTA = 60
ARCHAEA_QUOTA = 50
MAX_LEN = 700_000  # SPAdes/geNomad sanity ceiling; nothing real is excluded by this


def size_bin(n):
    return next(label for lo, hi, label in SIZE_BINS if lo <= n < hi)


def host_genus(name):
    words = [w for w in name.replace("[", "").replace("]", "").split()
             if w not in ("Candidatus", "uncultured", "unclassified")]
    return words[0] if words else "?"


def load(paths):
    recs = {}
    for p in paths:
        with open(p) as f:
            for line in f:
                try:
                    d = json.loads(line)
                except json.JSONDecodeError:
                    continue
                acc = d.get("accession")
                host_lin = d.get("host", {}).get("lineage", [])
                host_ids = {x["tax_id"] for x in host_lin}
                domain = "Bacteria" if 2 in host_ids else "Archaea" if 2157 in host_ids else None
                length = d.get("length") or 0
                if not acc or not domain or not (1_000 <= length <= MAX_LEN):
                    continue
                vlin = [x["name"] for x in d.get("virus", {}).get("lineage", [])]
                host_names = [x["name"] for x in host_lin]
                recs[acc] = {
                    "accession": acc,
                    "virus_name": d.get("virus", {}).get("organism_name", ""),
                    "isolate": d.get("isolate", {}).get("name", ""),
                    "realm": vlin[1] if len(vlin) > 1 else "",
                    "class": next((x for x in vlin if x.endswith("viricetes")), "unclassified"),
                    "family": next((x for x in vlin if x.endswith("viridae")), "unclassified"),
                    "genus": next((x for x in vlin if x.endswith("virus") and " " not in x), ""),
                    "host_domain": domain,
                    "host_phylum": next((x for x in host_names if x.endswith("ota") and x not in ("Eukaryota",)), ""),
                    "host_genus": host_genus(d.get("host", {}).get("organism_name") or "?"),
                    "length": length,
                    "size_bin": size_bin(length),
                }
    return recs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("catalogs", nargs="+", help="datasets virus summary JSON Lines file(s)")
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20260928)
    ap.add_argument("--out", default="assets/stress_breadth_panel.tsv")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    recs = load(args.catalogs)
    pool = sorted(recs.values(), key=lambda r: r["accession"])
    rng.shuffle(pool)
    print(f"{len(pool)} eligible genomes (bacterial/archaeal host, complete, RefSeq)")

    chosen, why = {}, {}

    def take(r, reason):
        if r["accession"] not in chosen and len(chosen) < args.n:
            chosen[r["accession"]] = r
            why[r["accession"]] = reason
            return True
        return False

    def cover(key, reason, filt=lambda r: True):
        seen = {r[key] for r in chosen.values()}
        for r in pool:
            if filt(r) and r[key] not in seen and take(r, reason):
                seen.add(r[key])

    # 1. every non-tailed class, up to 10 each (they are rare and the most likely to break things)
    per_class = defaultdict(int)
    for r in pool:
        if r["class"] != "Caudoviricetes" and per_class[r["class"]] < 10 and take(r, f"class:{r['class']}"):
            per_class[r["class"]] += 1
    # 2-3. every family, every host genus
    cover("family", "family")
    cover("host_genus", "host_genus")
    # 4. jumbo + archaeal quotas
    for label, quota, filt in (("jumbo", JUMBO_QUOTA, lambda r: r["length"] >= 200_000),
                               ("archaea", ARCHAEA_QUOTA, lambda r: r["host_domain"] == "Archaea")):
        have = sum(filt(r) for r in chosen.values())
        for r in pool:
            if have >= quota:
                break
            if filt(r) and take(r, label):
                have += 1
    # 5. fill: round-robin over host genera, preferring unseen family x size-bin combos
    by_genus = defaultdict(list)
    for r in pool:
        if r["accession"] not in chosen:
            by_genus[r["host_genus"]].append(r)
    combos = {(r["family"], r["size_bin"]) for r in chosen.values()}
    for g in by_genus:
        by_genus[g].sort(key=lambda r: (r["family"], r["size_bin"]) in combos)
    genera = sorted(by_genus)
    rng.shuffle(genera)
    while len(chosen) < args.n and any(by_genus.values()):
        for g in genera:
            if by_genus[g] and len(chosen) < args.n:
                r = by_genus[g].pop(0)
                take(r, "fill")
                combos.add((r["family"], r["size_bin"]))

    cols = ["sample", "accession", "virus_name", "isolate", "realm", "class", "family", "genus",
            "host_domain", "host_phylum", "host_genus", "length", "size_bin", "selected_for"]
    rows = sorted(chosen.values(), key=lambda r: (r["class"], r["family"], r["length"]))
    with open(args.out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            r = dict(r, sample="bx_" + r["accession"].replace(".", "_"), selected_for=why[r["accession"]])
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    def n_unique(k):
        return len({r[k] for r in chosen.values()})
    print(f"Selected {len(chosen)}: {n_unique('class')} classes, {n_unique('family')} families, "
          f"{n_unique('host_genus')} host genera, {n_unique('host_phylum')} host phyla, "
          f"{sum(r['host_domain'] == 'Archaea' for r in chosen.values())} archaeal, "
          f"{sum(r['length'] >= 200_000 for r in chosen.values())} jumbo")
    for _, _, label in SIZE_BINS:
        print(f"  {label:>10}: {sum(r['size_bin'] == label for r in chosen.values())}")
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
