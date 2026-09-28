#!/usr/bin/env python3
"""
Download + verify reference genomes for the PHINDER stress test.

Each reference is fetched from NCBI nuccore (E-utilities) and VERIFIED before use:
the FASTA title must contain an expected keyword and the total length must be within
tolerance of the expected size. Candidate accessions are tried in order; the first that
verifies wins. Anything that fails verification is reported and left out — fix the
accession here rather than trusting an unverified genome.

Usage:
    python3 bin/stress_fetch_refs.py --outdir stress_data/refs
    python3 bin/stress_fetch_refs.py --panel assets/stress_breadth_panel.tsv --outdir stress_data/breadth_refs

--panel mode downloads every accession in a breadth panel (batched, 100 per request) and verifies
each record's accession and length against the panel before writing <sample>.fasta.

Outputs:
    <outdir>/<ref_id>.fasta
    <outdir>/refs_manifest.tsv   (ref_id, accessions, title, length, topology, lifestyle, genome_type)
"""
import argparse
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id="

# ref_id: candidates (each candidate = list of accessions forming the genome, e.g. Phi6 segments),
#         title keywords (any match; NCBI titles may use ICTV binomials), expected total bp, topology for read simulation, lifestyle truth, genome type
REFS = {
    "phiX174":   dict(cands=[["NC_001422"]], kw=("phix174", "phi x174", "sinsheimervirus"), size=5386,   topo="circular", life="lytic",     gtype="ssDNA"),
    "M13":       dict(cands=[["NC_003287"]], kw=("phage m13", "m13"),     size=6407,   topo="circular", life="NA",        gtype="ssDNA (filamentous, chronic)"),
    "MS2":       dict(cands=[["NC_001417"]], kw=("ms2", "emesvirus zinderi"),     size=3569,   topo="linear",   life="lytic",     gtype="ssRNA"),
    "Phi6":      dict(cands=[["NC_003714", "NC_003715", "NC_003716"]], kw=("phi6", "phi-6", "cystovirus"), size=13385, topo="linear", life="lytic", gtype="dsRNA (3 segments)"),
    "PM2":       dict(cands=[["NC_000867"]], kw=("pm2", "corticovirus"),     size=10079,  topo="circular", life="lytic",     gtype="dsDNA (tailless, lipid)"),
    "PRD1":      dict(cands=[["NC_001421"]], kw=("prd1", "alphatectivirus"),    size=14927,  topo="linear",   life="lytic",     gtype="dsDNA (tailless, lipid)"),
    "T7":        dict(cands=[["NC_001604"]], kw=("phage t7", "virus t7", "teseptimavirus"),      size=39937,  topo="linear",   life="lytic",     gtype="dsDNA"),
    "lambda":    dict(cands=[["NC_001416"]], kw=("lambda", "lambdavirus"),  size=48502,  topo="linear",   life="temperate", gtype="dsDNA"),
    "P22":       dict(cands=[["NC_002371"]], kw=("phage p22", "virus p22", "lederbergvirus"),     size=41724,  topo="linear",   life="temperate", gtype="dsDNA"),
    "Mu":        dict(cands=[["NC_000929"]], kw=("phage mu", "virus mu", "muvirus"),      size=36717,  topo="linear",   life="temperate", gtype="dsDNA (transposable)"),
    "P1":        dict(cands=[["NC_005856"]], kw=("phage p1", "virus p1", "punavirus"),      size=94800,  topo="linear",   life="temperate", gtype="dsDNA (plasmid prophage)"),
    "L5":        dict(cands=[["NC_001335"]], kw=("phage l5", "fromanvirus"),      size=52297,  topo="linear",   life="temperate", gtype="dsDNA (high GC)"),
    "D29":       dict(cands=[["NC_001900"], ["AF022214"]], kw=("phage d29", "virus d29", "fromanvirus"), size=49136, topo="linear", life="lytic", gtype="dsDNA (high GC)"),
    "crAss001":  dict(cands=[["MH675552"]],  kw=("crass001",), size=102679, topo="linear",  life="NA",        gtype="dsDNA (crAss-like)"),
    "T5":        dict(cands=[["NC_005859"]], kw=("phage t5", "virus t5", "tequintavirus"),      size=121750, topo="linear",   life="lytic",     gtype="dsDNA"),
    "phageK":    dict(cands=[["NC_005880"], ["KF766114"]], kw=("phage k ", "phage k,", "kayvirus"), size=140000, topo="linear", life="lytic", gtype="dsDNA (low GC)"),
    "T4":        dict(cands=[["NC_000866"]], kw=("phage t4", "virus t4", "tequatrovirus"),      size=168903, topo="linear",   life="lytic",     gtype="dsDNA"),
    "phiKZ":     dict(cands=[["NC_004629"]], kw=("phikz", "phi kz", "phikzvirus"),   size=280334, topo="linear",   life="lytic",     gtype="dsDNA (jumbo)"),
    "phageG":    dict(cands=[["NC_023719"], ["JN638751"]], kw=("phage g ", "phage g,", "donellivirus"), size=497513, topo="linear", life="lytic", gtype="dsDNA (jumbo)"),
    # Host / negative control (not a phage)
    "Ecoli_K12": dict(cands=[["NC_000913"]], kw=("k-12",),    size=4641652, topo="circular", life="NA",       gtype="bacterial"),
}
SIZE_TOL = 0.20  # wrong organism is usually off by far more than 20%


def efetch(acc, retries=3):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(EUTILS + acc, timeout=120) as r:
                txt = r.read().decode()
            if txt.startswith(">"):
                return txt
        except (urllib.error.URLError, TimeoutError) as e:
            print(f"    {acc}: {e} (attempt {attempt + 1})", file=sys.stderr)
        time.sleep(2 * (attempt + 1))
    return None


def parse_fasta(txt):
    recs, hdr, seq = [], None, []
    for line in txt.splitlines():
        if line.startswith(">"):
            if hdr is not None:
                recs.append((hdr, "".join(seq)))
            hdr, seq = line[1:].strip(), []
        elif line.strip():
            seq.append(line.strip().upper())
    if hdr is not None:
        recs.append((hdr, "".join(seq)))
    return recs


# Classes whose genomes are circular (affects only where simulated fragments may span the origin)
CIRCULAR_CLASSES = {"Microviricetes", "Faserviricetes", "Laserviricetes"}


def fetch_panel(panel_tsv, out, batch=100):
    import csv
    with open(panel_tsv) as f:
        panel = list(csv.DictReader(f, delimiter="\t"))
    todo = [r for r in panel if not (out / f"{r['sample']}.fasta").exists()]
    print(f"{len(panel)} panel genomes, {len(panel) - len(todo)} already downloaded, fetching {len(todo)}")
    by_acc = {r["accession"].split(".")[0]: r for r in panel}
    failed = []
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        data = ("db=nuccore&rettype=fasta&retmode=text&id=" + ",".join(r["accession"] for r in chunk)).encode()
        txt = None
        for attempt in range(4):
            try:
                with urllib.request.urlopen(EUTILS.split("?")[0], data=data, timeout=300) as resp:
                    txt = resp.read().decode()
                if txt.startswith(">"):
                    break
            except (urllib.error.URLError, TimeoutError) as e:
                print(f"    batch {i // batch + 1}: {e} (attempt {attempt + 1})", file=sys.stderr)
            txt = None
            time.sleep(5 * (attempt + 1))
        got = set()
        for hdr, seq in parse_fasta(txt or ""):
            acc = hdr.split()[0].split(".")[0]
            r = by_acc.get(acc)
            if r is None:
                continue
            if abs(len(seq) - int(r["length"])) > 0.01 * int(r["length"]):
                print(f"    {r['accession']}: length {len(seq):,} != panel {int(r['length']):,} — skipped")
                continue
            with open(out / f"{r['sample']}.fasta", "w") as f:
                f.write(f">{hdr.split()[0]}\n")
                for j in range(0, len(seq), 80):
                    f.write(seq[j:j + 80] + "\n")
            got.add(r["sample"])
        failed += [r["sample"] for r in chunk if r["sample"] not in got]
        print(f"  batch {i // batch + 1}/{(len(todo) + batch - 1) // batch}: {len(got)}/{len(chunk)} ok")
        time.sleep(0.5)

    # Manifest in the same format as the curated refs, so stress_simulate.py can load either
    n = 0
    with open(out / "refs_manifest.tsv", "w") as f:
        f.write("ref_id\taccessions\ttitle\tlength\ttopology\tlifestyle\tgenome_type\n")
        for r in panel:
            if (out / f"{r['sample']}.fasta").exists():
                topo = "circular" if r["class"] in CIRCULAR_CLASSES else "linear"
                f.write(f"{r['sample']}\t{r['accession']}\t{r['virus_name']}\t{r['length']}\t{topo}\tNA\t"
                        f"{r['class']}/{r['family']}\n")
                n += 1
    print(f"\n{n}/{len(panel)} panel genomes ready -> {out}/refs_manifest.tsv")
    if failed:
        print(f"FAILED ({len(failed)}) — rerun the same command to retry: {', '.join(failed[:20])}"
              f"{' ...' if len(failed) > 20 else ''}")
        sys.exit(1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="stress_data/refs")
    ap.add_argument("--only", nargs="*", help="Fetch only these ref_ids")
    ap.add_argument("--panel", help="Breadth panel TSV (bin/stress_breadth_panel.py) — fetch those instead")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    if args.panel:
        return fetch_panel(args.panel, out)
    manifest, failed = [], []

    for ref_id, spec in REFS.items():
        if args.only and ref_id not in args.only:
            continue
        fasta = out / f"{ref_id}.fasta"
        print(f"[{ref_id}]")
        ok = None
        for cand in spec["cands"]:
            recs = []
            for acc in cand:
                txt = efetch(acc)
                time.sleep(0.4)  # stay under NCBI's 3 req/s without an API key
                if txt is None:
                    recs = None
                    break
                recs += parse_fasta(txt)
            if not recs:
                print(f"    {'+'.join(cand)}: download failed")
                continue
            total = sum(len(s) for _, s in recs)
            titles = " | ".join(h for h, _ in recs)
            kw_ok = all(any(k in h.lower() + " " for k in spec["kw"]) for h, _ in recs)
            size_ok = abs(total - spec["size"]) <= SIZE_TOL * spec["size"]
            status = "OK" if kw_ok and size_ok else "MISMATCH"
            print(f"    {'+'.join(cand)}: {total:,} bp  [{status}]  {titles[:110]}")
            if kw_ok and size_ok:
                ok = (cand, recs, total, titles)
                break

        if ok is None:
            failed.append(ref_id)
            continue
        cand, recs, total, titles = ok
        with open(fasta, "w") as f:
            for h, s in recs:
                f.write(f">{h.split()[0]}\n")
                for i in range(0, len(s), 80):
                    f.write(s[i:i + 80] + "\n")
        manifest.append([ref_id, "+".join(cand), titles, str(total), spec["topo"], spec["life"], spec["gtype"]])

    with open(out / "refs_manifest.tsv", "w") as f:
        f.write("ref_id\taccessions\ttitle\tlength\ttopology\tlifestyle\tgenome_type\n")
        for row in manifest:
            f.write("\t".join(row) + "\n")

    print(f"\n{len(manifest)} references verified -> {out}/refs_manifest.tsv")
    if failed:
        print(f"FAILED VERIFICATION ({len(failed)}): {', '.join(failed)}")
        print("Fix the accession(s) in REFS before simulating; samples using these refs will be skipped.")
        sys.exit(1)


if __name__ == "__main__":
    main()
