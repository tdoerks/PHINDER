#!/usr/bin/env python3
"""
Score a PHINDER stress-test run against the ground truth from bin/stress_simulate.py.

Per sample it reports:
  * module status   — from the Nextflow trace (retried-then-ignored failures show as FAILED)
  * assembly        — contigs, largest contig, k-mer recovery of each expected reference,
                      purity (share of assembled bp from contigs that match an expected phage)
  * detection       — CheckV best quality/completeness, geNomad virus count + viral fraction
  * lifestyle       — BacPhlip and VIBRANT calls (first row, as the dashboard shows) vs truth
  * dashboard       — whether the sample made it into summary/phinder_summary.tsv
  * taxonomy        — breadth tier: geNomad's class/family vs the NCBI lineage of the input genome
  * breakdown       — pass/failure rates grouped by class, family, host genus/domain and size bin

Module outputs are read by EXACT path: a missing file is reported as missing, never
substituted with another sample's result.

Usage:
    python3 bin/score_stress_test.py --truth stress_data/stress_truth.tsv --refs stress_data/refs \
        --run reads=stress_runs/reads/results --run assembly=stress_runs/assembly/results \
        --run sra=stress_runs/sra/results --out stress_runs/stress_scorecard
    (--truth and --refs are repeatable; runs whose directory does not exist are skipped)
"""
import argparse
import csv
import html
from collections import defaultdict
from pathlib import Path

K = 21
RECOVERY_PASS = 90.0      # % of reference k-mers recovered
NEG_VIRAL_FRAC_MAX = 20.0  # negatives: geNomad-viral share of assembly must stay below this
SAMPLE_MODULES = ["DOWNLOAD_SRA", "FASTQC", "FASTP", "SPADES", "QUAST", "CHECKV", "PHAROKKA",
                  "VIBRANT", "DIAMOND_PROPHAGE", "PHANOTATE", "BACPHLIP", "AMRFINDERPLUS",
                  "GENOMAD", "PHAGETERM"]
QUALITY_RANK = {"Complete": 5, "High-quality": 4, "Medium-quality": 3, "Low-quality": 2, "Not-determined": 1}
COMP = str.maketrans("ACGTN", "TGCAN")


# ─── Sequence helpers ─────────────────────────────────────────────────────────

def read_fasta(path):
    recs, name, seq = [], None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(seq)))
                name, seq = line[1:].split()[0] if line[1:].strip() else "", []
            else:
                seq.append(line.strip().upper())
    if name is not None:
        recs.append((name, "".join(seq)))
    return recs


class RefKmers:
    """Forward + reverse-complement k-mer sets of one reference (cached)."""
    _cache = {}

    def __new__(cls, fasta):
        if fasta not in cls._cache:
            obj = super().__new__(cls)
            obj.fwd, obj.rc = set(), set()
            for _, s in read_fasta(fasta):
                for i in range(len(s) - K + 1):
                    obj.fwd.add(s[i:i + K])
            obj.rc = {k.translate(COMP)[::-1] for k in obj.fwd}
            cls._cache[fasta] = obj
        return cls._cache[fasta]


def assembly_vs_refs(contigs, ref_fastas):
    """Return ({ref: recovery%}, purity%) — purity = assembled bp in contigs with >=50% k-mers
    from any expected reference, over all assembled bp."""
    refs = {r: RefKmers(p) for r, p in ref_fastas.items()}
    found = {r: set() for r in refs}
    total_bp, matched_bp = 0, 0
    for _, s in contigs:
        total_bp += len(s)
        n = hits = 0
        for i in range(len(s) - K + 1):
            km = s[i:i + K]
            n += 1
            for r, rk in refs.items():
                if km in rk.fwd:
                    found[r].add(km); hits += 1; break
                if km in rk.rc:
                    found[r].add(km.translate(COMP)[::-1]); hits += 1; break
        if n and hits / n >= 0.5:
            matched_bp += len(s)
    recovery = {r: 100.0 * len(found[r]) / max(1, len(refs[r].fwd)) for r in refs}
    purity = 100.0 * matched_bp / total_bp if total_bp else None
    return recovery, purity


# ─── Output parsers (strict paths) ────────────────────────────────────────────

def parse_trace(run_dir):
    """{tag: {process: (status, exit, realtime, peak_rss)}}, plus run-level (untagged) tasks."""
    trace = Path(run_dir) / "pipeline_trace.txt"
    per_sample, run_level = defaultdict(dict), {}
    if not trace.exists():
        return None, None
    with open(trace) as f:
        rows = sorted(csv.DictReader(f, delimiter="\t"), key=lambda r: int(r.get("task_id") or 0))
    for r in rows:
        name = r["name"]
        proc = name.split(" (")[0].split(":")[-1]
        tag = name[name.find("(") + 1:name.rfind(")")] if "(" in name else None
        rec = (r["status"], r.get("exit", ""), r.get("realtime", ""), r.get("peak_rss", ""))
        target = per_sample[tag] if tag and proc in SAMPLE_MODULES else run_level
        prev = target.get(proc)
        # a later attempt overrides, but never let a failure mask an earlier success
        if prev is None or prev[0] not in ("COMPLETED", "CACHED"):
            target[proc] = rec
    return per_sample, run_level


def first_row(path, delim="\t"):
    if not path or not Path(path).exists():
        return None
    with open(path) as f:
        return next(csv.DictReader(f, delimiter=delim), None)


def parse_checkv(run_dir, s):
    p = Path(run_dir) / "checkv" / f"{s}_checkv" / "quality_summary.tsv"
    if not p.exists():
        return None
    with open(p) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if not rows:
        return {"best_quality": "no contigs", "completeness": None, "n_good": 0}

    def key(r):
        try:
            c = float(r.get("completeness") or 0)
        except ValueError:
            c = 0
        return (QUALITY_RANK.get(r.get("checkv_quality", ""), 0), c)
    best = max(rows, key=key)
    return {"best_quality": best.get("checkv_quality", "?"), "completeness": best.get("completeness", ""),
            "n_good": sum(QUALITY_RANK.get(r.get("checkv_quality", ""), 0) >= 3 for r in rows)}


def parse_genomad(run_dir, s, asm_bp):
    hits = list((Path(run_dir) / "genomad").glob(f"{s}/**/{s}_virus_summary.tsv"))
    if not hits:
        return None
    with open(hits[0]) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    viral_bp = sum(int(r.get("length") or 0) for r in rows)
    best = max(rows, key=lambda r: float(r.get("virus_score") or 0), default=None)
    tax = ""
    if best:
        parts = [p for p in (best.get("taxonomy") or "").split(";") if p]
        tax = parts[-1] if parts else "unclassified"
    return {"n_virus": len(rows), "best_score": best.get("virus_score") if best else None, "taxon": tax,
            "lineage": (best.get("taxonomy") or "") if best else "",
            "viral_frac": 100.0 * viral_bp / asm_bp if asm_bp else None}


def parse_bacphlip(run_dir, s):
    row = first_row(Path(run_dir) / "bacphlip" / f"{s}.bacphlip")
    if not row:
        return None
    try:
        v = float(row.get("Virulent", 0))
    except ValueError:
        return None
    return {"call": "lytic" if v >= 0.5 else "temperate", "p_virulent": v}


def parse_vibrant(run_dir, s):
    hits = list((Path(run_dir) / "vibrant" / f"{s}_vibrant").glob("**/VIBRANT_genome_quality_*.tsv"))
    if not hits:
        return None
    with open(hits[0]) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if not rows:
        return {"call": "no virus", "n": 0}
    t = rows[0].get("type", "")
    return {"call": "temperate" if t in ("lysogenic", "temperate") else t, "n": len(rows)}


def parse_phageterm(run_dir, s):
    p = Path(run_dir) / "phageterm" / s / f"{s}_phageterm_strategy.txt"
    if not p.exists():
        return None
    lines = p.read_text().strip().splitlines()
    return lines[0].strip() if lines else "empty"


def dashboard_samples(run_dir):
    row_file = Path(run_dir) / "summary" / "phinder_summary.tsv"
    if not row_file.exists():
        return None
    with open(row_file) as f:
        return {r["sample_id"] for r in csv.DictReader(f, delimiter="\t")}


# ─── Scoring ──────────────────────────────────────────────────────────────────

def find_ref(refs_dirs, ref_id):
    for d in refs_dirs:
        p = Path(d) / f"{ref_id}.fasta"
        if p.exists():
            return str(p)
    return None


def taxonomy_match(genomad, t):
    """'family' / 'class' / 'none' — how deep geNomad's lineage agrees with NCBI's."""
    if not genomad or not genomad.get("lineage"):
        return "none"
    lin = genomad["lineage"].lower()
    if t.get("family") and t["family"] != "unclassified" and t["family"].lower() in lin:
        return "family"
    if t.get("class") and t["class"] != "unclassified" and t["class"].lower() in lin:
        return "class"
    return "none"


def score_sample(t, run_dir, refs_dirs, trace, dash):
    s, mode = t["sample"], t["mode"]
    res = {"sample": s, "mode": mode, "category": t["category"], "notes": t["notes"],
           "expect_phage": t["expect_phage"], "expect_lifestyle": t["expect_lifestyle"]}
    checks = {}

    # Module status
    mods = trace.get(s, {}) if trace is not None else {}
    res["modules"] = mods
    failed = [m for m, r in mods.items() if r[0] not in ("COMPLETED", "CACHED")]
    res["failed_modules"] = failed
    checks["modules"] = None if trace is None else (not failed and bool(mods))

    # Assembly
    if mode == "assembly":
        asm = Path(find_ref(refs_dirs, s[4:] if s.startswith("asm_") else s) or "/nonexistent")
    else:
        asm = Path(run_dir) / "assemblies" / f"{s}_assembly.fasta"
    contigs = read_fasta(asm) if asm.exists() else None
    comps = [c.split(":")[0] for c in t["refs"].split(";") if c]
    asm_bp = sum(len(q) for _, q in contigs) if contigs else 0
    if contigs is None:
        res.update(n_contigs=None, largest=None, asm_bp=None, recovery={}, purity=None)
        checks["assembly"] = False
    else:
        res.update(n_contigs=len(contigs), asm_bp=asm_bp,
                   largest=max((len(q) for _, q in contigs), default=0))
        ref_fastas = {r: find_ref(refs_dirs, r) for r in comps if find_ref(refs_dirs, r)}
        rec, pur = assembly_vs_refs(contigs, ref_fastas) if ref_fastas else ({}, None)
        res.update(recovery=rec, purity=pur)
        checks["assembly"] = all(v >= RECOVERY_PASS for v in rec.values()) if rec else (None if comps else True)

    # Detection
    res["checkv"] = parse_checkv(run_dir, s)
    res["genomad"] = parse_genomad(run_dir, s, asm_bp)
    g, c = res["genomad"], res["checkv"]
    if t["expect_phage"] == "yes":
        checks["checkv"] = bool(c) and QUALITY_RANK.get(c["best_quality"], 0) >= 3
        checks["genomad"] = bool(g) and g["n_virus"] > 0
    else:
        checks["genomad"] = bool(g) and (g["viral_frac"] or 0) < NEG_VIRAL_FRAC_MAX
        checks["checkv"] = None

    # Lifestyle
    res["bacphlip"] = parse_bacphlip(run_dir, s)
    res["vibrant"] = parse_vibrant(run_dir, s)
    exp = t["expect_lifestyle"]
    for tool in ("bacphlip", "vibrant"):
        r = res[tool]
        checks[tool] = None if exp not in ("lytic", "temperate") else (bool(r) and r["call"] == exp)

    res["meta"] = {k: t.get(k, "") for k in ("class", "family", "host_genus", "host_domain", "size_bin", "length")}
    if t.get("class"):
        res["tax_match"] = taxonomy_match(res["genomad"], t)
        # class is the fair bar: geNomad's database can lag newly created ICTV families
        checks["taxonomy"] = res["tax_match"] in ("family", "class") if t["class"] != "unclassified" else None
    else:
        res["tax_match"] = ""

    res["phageterm"] = parse_phageterm(run_dir, s)
    res["on_dashboard"] = None if dash is None else s in dash
    checks["dashboard"] = res["on_dashboard"]

    decided = [v for v in checks.values() if v is not None]
    res["checks"] = checks
    res["verdict"] = ("PASS" if decided and all(decided) else
                      "FAIL" if not checks.get("modules") and not checks.get("assembly") else "PARTIAL")
    return res


def fmt(v, nd=1):
    if v is None:
        return "–"
    return f"{v:.{nd}f}" if isinstance(v, float) else str(v)


def tsv_row(r):
    c, g, b, v = r["checkv"] or {}, r["genomad"] or {}, r["bacphlip"] or {}, r["vibrant"] or {}
    spades = r["modules"].get("SPADES", ("", "", "", ""))
    return {
        "sample": r["sample"], "mode": r["mode"], "category": r["category"], "verdict": r["verdict"],
        "failed_modules": ",".join(f"{m}({r['modules'][m][1]})" for m in r["failed_modules"]),
        "n_contigs": fmt(r["n_contigs"]), "largest_contig": fmt(r["largest"]), "assembly_bp": fmt(r["asm_bp"]),
        "recovery": ";".join(f"{k}:{v:.1f}" for k, v in r["recovery"].items()),
        "purity_pct": fmt(r["purity"]),
        "checkv_best": c.get("best_quality", "missing" if not r["checkv"] else ""),
        "checkv_completeness": c.get("completeness", ""), "checkv_n_med_plus": c.get("n_good", ""),
        "genomad_n_virus": g.get("n_virus", "missing" if not r["genomad"] else ""),
        "genomad_viral_frac_pct": fmt(g.get("viral_frac")), "genomad_taxon": g.get("taxon", ""),
        "expect_lifestyle": r["expect_lifestyle"], "bacphlip": b.get("call", "missing"),
        "bacphlip_p_virulent": fmt(b.get("p_virulent"), 3), "vibrant": v.get("call", "missing"),
        "phageterm": r["phageterm"] or "–", "on_dashboard": fmt(r["on_dashboard"]),
        "spades_realtime": spades[2], "spades_peak_rss": spades[3],
        "ncbi_class": r["meta"]["class"], "ncbi_family": r["meta"]["family"], "genomad_tax_match": r["tax_match"],
        "host_genus": r["meta"]["host_genus"], "host_domain": r["meta"]["host_domain"],
        "size_bin": r["meta"]["size_bin"], "genome_length": r["meta"]["length"], "notes": r["notes"],
    }


# ─── HTML ─────────────────────────────────────────────────────────────────────

CSS = """
:root{--bg:#fafaf8;--fg:#1d1d1b;--mut:#6b6b66;--line:#e2e1dc;--ok:#1f7a4d;--okbg:#e3f3ea;--bad:#b3261e;
--badbg:#fbe6e4;--mid:#8a6100;--midbg:#fbf1d9;--na:#9a9a94;--nabg:#f0efeb}
@media (prefers-color-scheme:dark){:root{--bg:#171716;--fg:#ecebe6;--mut:#a3a29c;--line:#34332f;--ok:#7fd3a4;
--okbg:#173626;--bad:#ff9b92;--badbg:#3d1c19;--mid:#f1c86b;--midbg:#3a2e12;--na:#7c7b75;--nabg:#242320}}
body{background:var(--bg);color:var(--fg);font:14px/1.45 system-ui,sans-serif;margin:0;padding:24px 16px}
h1{font-size:22px;margin:0 0 4px}h2{font-size:16px;margin:28px 0 8px}p.sub{color:var(--mut);margin:0 0 16px}
.wrap{overflow-x:auto;border:1px solid var(--line);border-radius:8px}
table{border-collapse:collapse;width:100%;font-size:12.5px}
th,td{padding:5px 7px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}
th{position:sticky;top:0;background:var(--bg);font-weight:600;color:var(--mut)}
th.rot{writing-mode:vertical-rl;transform:rotate(180deg);height:110px;padding:6px 3px}
td.m{text-align:center;padding:5px 3px;font-weight:600}
.ok{background:var(--okbg);color:var(--ok)}.bad{background:var(--badbg);color:var(--bad)}
.mid{background:var(--midbg);color:var(--mid)}.na{background:var(--nabg);color:var(--na)}
tr.cat td{background:var(--nabg);font-weight:600;color:var(--mut);text-transform:uppercase;font-size:11px;letter-spacing:.04em}
.pill{display:inline-block;padding:1px 8px;border-radius:10px;font-weight:600;font-size:11.5px}
.stats{display:flex;gap:12px;flex-wrap:wrap;margin:0 0 8px}.stat{border:1px solid var(--line);border-radius:8px;padding:8px 14px}
.stat b{display:block;font-size:20px}
"""


def cell(ok, text):
    cls = {True: "ok", False: "bad", None: "na"}[ok]
    return f'<td class="{cls}">{html.escape(text)}</td>'


def breakdown(results):
    """Rows of per-group rates for every grouping dimension present in the truth metadata."""
    out = []
    dims = [("class", "Virus class"), ("family", "Virus family"), ("host_domain", "Host domain"),
            ("host_genus", "Host genus"), ("size_bin", "Genome size"), ("category", "Stress category")]
    for key, label in dims:
        groups = defaultdict(list)
        for r in results:
            g = r["category"] if key == "category" else r["meta"].get(key)
            if g:
                groups[g].append(r)
        for g, rs in groups.items():
            n = len(rs)
            fails = defaultdict(int)
            for r in rs:
                for m in r["failed_modules"]:
                    fails[m] += 1

            def pct(pred):
                return round(100.0 * sum(1 for r in rs if pred(r)) / n, 1)
            out.append({
                "dimension": label, "group": g, "n": n,
                "pass_pct": pct(lambda r: r["verdict"] == "PASS"),
                "any_module_failed_pct": pct(lambda r: bool(r["failed_modules"])),
                "recovered_pct": pct(lambda r: r["checks"].get("assembly") is True),
                "checkv_med_plus_pct": pct(lambda r: r["checks"].get("checkv") is True),
                "genomad_detected_pct": pct(lambda r: bool(r["genomad"]) and r["genomad"]["n_virus"] > 0),
                "tax_class_match_pct": pct(lambda r: r.get("tax_match") in ("family", "class")),
                "top_failing_module": max(fails, key=fails.get) + f" ({fails[max(fails, key=fails.get)]})" if fails else "",
            })
    return out


def build_html(results, run_levels, out_path, groups=None):
    counts = defaultdict(int)
    for r in results:
        counts[r["verdict"]] += 1
    mods_present = [m for m in SAMPLE_MODULES if any(m in r["modules"] for r in results)]
    head = "".join(f'<th class="rot">{m}</th>' for m in mods_present)
    rows, last_cat = [], None
    for r in results:
        cat = f'{r["mode"]} · {r["category"]}'
        if cat != last_cat:
            rows.append(f'<tr class="cat"><td colspan="{len(mods_present) + 12}">{html.escape(cat)}</td></tr>')
            last_cat = cat
        vcls = {"PASS": "ok", "PARTIAL": "mid", "FAIL": "bad"}[r["verdict"]]
        tds = [f'<td title="{html.escape(r["notes"])}">{html.escape(r["sample"])}</td>',
               f'<td><span class="pill {vcls}">{r["verdict"]}</span></td>']
        for m in mods_present:
            st = r["modules"].get(m)
            if st is None:
                tds.append('<td class="m na">·</td>')
            elif st[0] in ("COMPLETED", "CACHED"):
                tds.append(f'<td class="m ok" title="{st[2]} {st[3]}">✓</td>')
            else:
                tds.append(f'<td class="m bad" title="{st[0]} exit {st[1]}">{html.escape(st[1] or "✗")}</td>')
        ch = r["checks"]
        rec = " ".join(f"{k} {v:.0f}%" for k, v in r["recovery"].items()) or ("no asm" if r["n_contigs"] is None else "–")
        tds.append(cell(ch.get("assembly"), rec))
        tds.append(f'<td>{fmt(r["n_contigs"])} / {fmt(r["largest"])}</td>')
        tds.append(f'<td>{fmt(r["purity"], 0)}</td>')
        c = r["checkv"]
        tds.append(cell(ch.get("checkv"), f'{c["best_quality"]} {c["completeness"] or ""}' if c else "missing"))
        g = r["genomad"]
        gtxt = f'{g["n_virus"]} · {fmt(g["viral_frac"], 0)}% · {g["taxon"]}' if g else "missing"
        tds.append(cell(ch.get("genomad"), gtxt))
        tds.append(f'<td>{html.escape(r["expect_lifestyle"])}</td>')
        b, v = r["bacphlip"], r["vibrant"]
        tds.append(cell(ch.get("bacphlip"), f'{b["call"]} ({b["p_virulent"]:.2f})' if b else "missing"))
        tds.append(cell(ch.get("vibrant"), v["call"] if v else "missing"))
        tds.append(f'<td>{html.escape(r["phageterm"] or "–")}</td>')
        tds.append(cell(r["on_dashboard"], {True: "yes", False: "MISSING", None: "no summary"}[r["on_dashboard"]]))
        rows.append("<tr>" + "".join(tds) + "</tr>")

    run_html = []
    for mode, lv in run_levels.items():
        if lv is None:
            run_html.append(f"<li><b>{mode}</b>: no trace file (run did not start?)</li>")
            continue
        items = ", ".join(f'{p} {"✓" if s[0] in ("COMPLETED", "CACHED") else "✗ " + s[1]}' for p, s in lv.items())
        run_html.append(f"<li><b>{mode}</b>: {html.escape(items) or 'no run-level tasks'}</li>")

    doc = f"""<title>PHINDER Stress Scorecard</title><style>{CSS}</style>
<h1>PHINDER stress test scorecard</h1>
<p class="sub">Simulated and real phage data with known answers. Hover a sample for notes, a ✓ for runtime/memory,
a red cell for the exit code. Recovery = % of reference {K}-mers present in the assembly (pass ≥ {RECOVERY_PASS:.0f}%).
Purity = % of assembled bp from contigs matching an expected phage.</p>
<div class="stats"><div class="stat"><b>{len(results)}</b>samples</div><div class="stat"><b>{counts['PASS']}</b>pass</div>
<div class="stat"><b>{counts['PARTIAL']}</b>partial</div><div class="stat"><b>{counts['FAIL']}</b>fail</div></div>
<div class="wrap"><table><thead><tr><th>Sample</th><th>Verdict</th>{head}<th>Recovery</th><th>Contigs / largest</th>
<th>Purity %</th><th>CheckV best</th><th>geNomad n · viral% · taxon</th><th>Expected</th><th>BacPhlip</th><th>VIBRANT</th>
<th>PhageTerm</th><th>Dashboard</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div>
<h2>Run-level tasks</h2><ul>{''.join(run_html)}</ul>
{breakdown_html(groups or [])}
"""
    Path(out_path).write_text(doc)


def breakdown_html(groups):
    if not groups:
        return ""
    parts = ['<h2>Where it breaks — by group</h2><p class="sub">Sorted worst first within each dimension '
             '(groups with n ≥ 3; full table in the _breakdown.tsv).</p>']
    by_dim = defaultdict(list)
    for g in groups:
        by_dim[g["dimension"]].append(g)
    for dim, gs in by_dim.items():
        gs = sorted((g for g in gs if g["n"] >= 3), key=lambda g: (g["pass_pct"], -g["n"]))
        if not gs:
            continue
        rows = []
        for g in gs[:40]:
            def c(v, good=90):
                cls = "ok" if v >= good else "mid" if v >= 50 else "bad"
                return f'<td class="{cls}">{v:.0f}%</td>'
            rows.append(f'<tr><td>{html.escape(str(g["group"]))}</td><td>{g["n"]}</td>{c(g["pass_pct"])}'
                        f'{c(100 - g["any_module_failed_pct"])}{c(g["recovered_pct"])}{c(g["checkv_med_plus_pct"])}'
                        f'{c(g["genomad_detected_pct"])}{c(g["tax_class_match_pct"])}'
                        f'<td>{html.escape(g["top_failing_module"])}</td></tr>')
        parts.append(f'<h2>{html.escape(dim)}</h2><div class="wrap"><table><thead><tr><th>Group</th><th>n</th>'
                     '<th>Pass</th><th>No module failed</th><th>Recovered ≥90%</th><th>CheckV ≥ medium</th>'
                     '<th>geNomad detected</th><th>Taxonomy class ok</th><th>Top failing module</th></tr></thead>'
                     f'<tbody>{"".join(rows)}</tbody></table></div>')
    return "".join(parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", action="append", required=True, help="truth TSV (repeatable)")
    ap.add_argument("--refs", action="append", required=True, help="reference FASTA dir (repeatable)")
    ap.add_argument("--run", action="append", required=True, help="mode=results_dir (repeatable)")
    ap.add_argument("--out", default="stress_scorecard", help="Output prefix (.tsv and .html)")
    args = ap.parse_args()

    runs = {m: d for m, d in (r.split("=", 1) for r in args.run) if Path(d).is_dir()}
    print(f"Scoring runs: {', '.join(runs) or 'none found'}")
    traces, run_levels, dashes = {}, {}, {}
    for mode, d in runs.items():
        traces[mode], run_levels[mode] = parse_trace(d)
        dashes[mode] = dashboard_samples(d)

    truth = []
    for tp in args.truth:
        if Path(tp).exists():
            with open(tp) as f:
                truth += [t for t in csv.DictReader(f, delimiter="\t") if t["mode"] in runs]

    results = []
    for t in truth:
        r = score_sample(t, runs[t["mode"]], args.refs, traces[t["mode"]], dashes[t["mode"]])
        results.append(r)
        print(f'{r["verdict"]:<8} {r["sample"]:<24} failed={",".join(r["failed_modules"]) or "-"} '
              f'recovery={";".join(f"{k}:{v:.0f}" for k, v in r["recovery"].items()) or "-"}')

    rows = [tsv_row(r) for r in results]
    with open(f"{args.out}.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["sample"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    groups = breakdown(results)
    if groups:
        with open(f"{args.out}_breakdown.tsv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(groups[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(groups)
    build_html(results, run_levels, f"{args.out}.html", groups)
    print(f"\n-> {args.out}.tsv\n-> {args.out}_breakdown.tsv\n-> {args.out}.html")


if __name__ == "__main__":
    main()
