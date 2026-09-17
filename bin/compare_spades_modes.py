#!/usr/bin/env python3
"""
compare_spades_modes.py
Aggregate CheckV + QUAST results across all SPAdes mode × sample combinations
and produce a comparison TSV and HTML report.

Directory naming convention (set by spades_mode_compare.nf):
  checkv_results/<sample>__<mode>_checkv/quality_summary.tsv
  quast_results/<sample>__<mode>_quast/report.tsv

The samplesheet CSV (sample,fastq_1,fastq_2,sample_type) maps sample → type.
"""

import argparse
import csv
import os
import re
import sys
from pathlib import Path

MODES = ['isolate', 'meta', 'metaviral', 'rnaviral', 'careful', 'standard']

CHECKV_TIERS = ['Complete', 'High-quality', 'Medium-quality', 'Low-quality', 'Not-determined']

# ── helpers ───────────────────────────────────────────────────────────────────

def parse_samplesheet(path):
    """Return {sample_id: sample_type}."""
    types = {}
    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            types[row['sample'].strip()] = row.get('sample_type', 'unknown').strip()
    return types


def parse_checkv(checkv_dir: Path, effective_id: str):
    """
    Parse quality_summary.tsv from a CheckV output directory.
    Returns a dict with summary counts per quality tier and genome types found.
    effective_id: 'sample__mode'
    """
    tsv = checkv_dir / 'quality_summary.tsv'
    result = {
        'n_contigs': 0,
        'n_complete': 0,
        'n_high_quality': 0,
        'n_medium_quality': 0,
        'n_low_quality': 0,
        'n_not_determined': 0,
        'checkv_genome_types': '',
        'assembly_status': 'failed',
        'total_phage_bp': 0,
    }
    if not tsv.exists():
        return result

    tiers = []
    genome_types = set()
    total_bp = 0
    n = 0
    with open(tsv) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            n += 1
            q = row.get('checkv_quality', 'Not-determined')
            tiers.append(q)
            gt = row.get('genome_copies', '') or row.get('host_genes', '')
            # genome type lives in 'gene_count' metadata — look for dtR / atr flags
            # Real field is typically absent; use 'viral_genes' as proxy for 'is viral'
            viral_genes = int(row.get('viral_genes', 0) or 0)
            if viral_genes > 0:
                genome_types.add('viral')
            try:
                bp = int(row.get('contig_length', 0) or 0)
                if q in ('Complete', 'High-quality', 'Medium-quality'):
                    total_bp += bp
            except ValueError:
                pass

    result['n_contigs'] = n
    result['n_complete'] = tiers.count('Complete')
    result['n_high_quality'] = tiers.count('High-quality')
    result['n_medium_quality'] = tiers.count('Medium-quality')
    result['n_low_quality'] = tiers.count('Low-quality')
    result['n_not_determined'] = tiers.count('Not-determined')
    result['checkv_genome_types'] = ','.join(sorted(genome_types)) or 'none'
    result['assembly_status'] = 'success'
    result['total_phage_bp'] = total_bp
    return result


def parse_quast(quast_dir: Path):
    """Parse QUAST report.tsv. Returns key assembly metrics dict."""
    result = {
        'quast_n_contigs': '',
        'quast_largest_contig': '',
        'quast_total_length': '',
        'quast_N50': '',
        'quast_N90': '',
        'quast_GC_pct': '',
    }
    tsv = quast_dir / 'report.tsv'
    if not tsv.exists():
        return result

    field_map = {
        '# contigs (>= 0 bp)': 'quast_n_contigs',
        '# contigs': 'quast_n_contigs',
        'Largest contig': 'quast_largest_contig',
        'Total length': 'quast_total_length',
        'Total length (>= 0 bp)': 'quast_total_length',
        'N50': 'quast_N50',
        'N90': 'quast_N90',
        'GC (%)': 'quast_GC_pct',
    }
    with open(tsv) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                key = parts[0].strip()
                if key in field_map:
                    result[field_map[key]] = parts[1].strip()
    return result


# ── main ──────────────────────────────────────────────────────────────────────

def build_table(checkv_root, quast_root, sample_types):
    rows = []

    # Discover all effective IDs from checkv dirs: <sample>__<mode>_checkv
    for d in sorted(Path(checkv_root).iterdir()):
        if not d.is_dir():
            continue
        # Strip trailing _checkv
        name = d.name
        if name.endswith('_checkv'):
            effective_id = name[:-len('_checkv')]
        else:
            effective_id = name

        # Parse sample and mode from double-underscore separator
        parts = effective_id.split('__')
        if len(parts) < 2:
            continue
        mode = parts[-1]
        sample_id = '__'.join(parts[:-1])

        if mode not in MODES:
            continue

        sample_type = sample_types.get(sample_id, 'unknown')

        checkv_stats = parse_checkv(d, effective_id)

        # Find matching QUAST dir
        quast_dir = Path(quast_root) / f'{effective_id}_quast'
        quast_stats = parse_quast(quast_dir)

        # Did assembly actually succeed? (stub FASTA has ">ASSEMBLY_FAILED" header)
        assembly_ok = checkv_stats['assembly_status'] == 'success' and checkv_stats['n_contigs'] > 0

        rows.append({
            'sample': sample_id,
            'sample_type': sample_type,
            'mode': mode,
            'assembly_status': 'success' if assembly_ok else 'failed',
            **checkv_stats,
            **quast_stats,
        })

    return rows


def write_tsv(rows, out_path):
    if not rows:
        print("WARNING: no rows to write", file=sys.stderr)
        Path(out_path).write_text('')
        return
    fieldnames = list(rows[0].keys())
    with open(out_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows → {out_path}")


def _tier_color(n_complete, n_hq, n_mq, n_lq):
    """Return CSS background colour representing best tier achieved."""
    if n_complete > 0:
        return '#15803d'   # green  — complete
    if n_hq > 0:
        return '#2563eb'   # blue   — high quality
    if n_mq > 0:
        return '#d97706'   # amber  — medium quality
    if n_lq > 0:
        return '#6b7280'   # grey   — low quality
    return '#991b1b'       # red    — failed / no phage recovered


def write_html(rows, out_path, samplesheet_path):
    if not rows:
        Path(out_path).write_text('<p>No data.</p>')
        return

    # Organise data: {sample: {mode: row}}
    samples_ordered = []
    by_sample = {}
    for r in rows:
        s = r['sample']
        if s not in by_sample:
            by_sample[s] = {}
            samples_ordered.append(s)
        by_sample[s][r['mode']] = r

    # Heatmap table (sample × mode, coloured by best CheckV tier)
    heatmap_rows = []
    for s in samples_ordered:
        stype = next(iter(by_sample[s].values()))['sample_type']
        cells = [f'<td><b>{s}</b><br><small style="color:#9ca3af">{stype}</small></td>']
        for m in MODES:
            r = by_sample[s].get(m)
            if r:
                bg = _tier_color(r['n_complete'], r['n_high_quality'],
                                 r['n_medium_quality'], r['n_low_quality'])
                label_parts = []
                if r['n_complete']:
                    label_parts.append(f"{r['n_complete']}✓")
                if r['n_high_quality']:
                    label_parts.append(f"{r['n_high_quality']}HQ")
                if r['n_medium_quality']:
                    label_parts.append(f"{r['n_medium_quality']}MQ")
                label = ' '.join(label_parts) or ('ok' if r['assembly_status'] == 'success' else 'fail')
                n50 = r.get('quast_N50', '')
                title = f"N50={n50} | contigs={r['quast_n_contigs']}"
                cells.append(
                    f'<td style="background:{bg};color:#fff;text-align:center;padding:6px 10px" title="{title}">'
                    f'{label}</td>'
                )
            else:
                cells.append('<td style="background:#1f2937;color:#6b7280;text-align:center">—</td>')
        heatmap_rows.append('<tr>' + ''.join(cells) + '</tr>')

    # Detailed table
    detail_rows = []
    for r in rows:
        bg = _tier_color(r['n_complete'], r['n_high_quality'],
                         r['n_medium_quality'], r['n_low_quality'])
        status_badge = (f'<span style="background:{bg};color:#fff;padding:2px 8px;border-radius:4px">'
                        f'{r["assembly_status"]}</span>')
        tier_summary = (f"C:{r['n_complete']} HQ:{r['n_high_quality']} "
                        f"MQ:{r['n_medium_quality']} LQ:{r['n_low_quality']}")
        detail_rows.append(
            f"<tr>"
            f"<td>{r['sample']}</td>"
            f"<td>{r['sample_type']}</td>"
            f"<td><code>{r['mode']}</code></td>"
            f"<td>{status_badge}</td>"
            f"<td>{r['quast_n_contigs']}</td>"
            f"<td>{r['quast_N50']}</td>"
            f"<td>{r['quast_largest_contig']}</td>"
            f"<td>{r['quast_total_length']}</td>"
            f"<td>{tier_summary}</td>"
            f"<td>{r['n_complete']}</td>"
            f"<td>{r['total_phage_bp']}</td>"
            f"</tr>"
        )

    mode_headers = ''.join(f'<th>{m}</th>' for m in MODES)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>SPAdes Mode Comparison — PHINDER</title>
<style>
  body {{ background:#111827; color:#e5e7eb; font-family:'Segoe UI',sans-serif; margin:0; padding:24px; }}
  h1 {{ color:#f9fafb; font-size:1.5em; margin-bottom:4px; }}
  h2 {{ color:#9ca3af; font-size:1em; font-weight:400; margin-bottom:24px; }}
  h3 {{ color:#d1d5db; font-size:1.1em; margin-top:32px; margin-bottom:8px; }}
  table {{ border-collapse:collapse; width:100%; margin-bottom:32px; }}
  th {{ background:#1f2937; color:#9ca3af; font-weight:600; padding:8px 12px; text-align:left; font-size:0.82em; text-transform:uppercase; letter-spacing:.05em; }}
  td {{ border-bottom:1px solid #1f2937; padding:7px 12px; font-size:0.88em; vertical-align:middle; }}
  tr:hover td {{ background:#1f2937; }}
  .legend {{ display:flex; gap:16px; margin-bottom:16px; flex-wrap:wrap; }}
  .leg {{ display:flex; align-items:center; gap:6px; font-size:0.82em; }}
  .dot {{ width:14px; height:14px; border-radius:3px; }}
</style>
</head>
<body>
<h1>SPAdes Mode Comparison</h1>
<h2>PHINDER — phage recovery across 6 SPAdes assembly modes</h2>

<div class="legend">
  <div class="leg"><div class="dot" style="background:#15803d"></div> Complete genome</div>
  <div class="leg"><div class="dot" style="background:#2563eb"></div> High-quality (&gt;90%)</div>
  <div class="leg"><div class="dot" style="background:#d97706"></div> Medium-quality (50–90%)</div>
  <div class="leg"><div class="dot" style="background:#6b7280"></div> Low-quality (&lt;50%)</div>
  <div class="leg"><div class="dot" style="background:#991b1b"></div> Failed / no phage</div>
</div>

<h3>Heatmap — best CheckV tier per sample × mode</h3>
<p style="color:#9ca3af;font-size:0.82em">Hover cells for N50 and contig count. Numbers: ✓=complete, HQ=high-quality, MQ=medium-quality.</p>
<table>
<thead><tr><th>Sample (type)</th>{mode_headers}</tr></thead>
<tbody>
{''.join(heatmap_rows)}
</tbody>
</table>

<h3>Full detail table</h3>
<table>
<thead>
<tr>
  <th>Sample</th><th>Type</th><th>Mode</th><th>Status</th>
  <th>Contigs</th><th>N50</th><th>Largest</th><th>Total bp</th>
  <th>CheckV tiers</th><th>Complete</th><th>Phage bp</th>
</tr>
</thead>
<tbody>
{''.join(detail_rows)}
</tbody>
</table>

<p style="color:#4b5563;font-size:0.75em;margin-top:40px">
  Generated by PHINDER compare_spades_modes.py — samplesheet: {samplesheet_path}
</p>
</body>
</html>
"""
    Path(out_path).write_text(html)
    print(f"Wrote HTML report → {out_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--checkv-dir',   required=True, help='Dir containing per-sample CheckV output dirs')
    ap.add_argument('--quast-dir',    required=True, help='Dir containing per-sample QUAST output dirs')
    ap.add_argument('--samplesheet',  required=True, help='CSV with sample,fastq_1,fastq_2,sample_type')
    ap.add_argument('--outdir',       default='.',   help='Output directory')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    sample_types = parse_samplesheet(args.samplesheet)

    rows = build_table(args.checkv_dir, args.quast_dir, sample_types)

    tsv_out  = os.path.join(args.outdir, 'spades_mode_comparison.tsv')
    html_out = os.path.join(args.outdir, 'spades_mode_comparison.html')

    write_tsv(rows, tsv_out)
    write_html(rows, html_out, args.samplesheet)


if __name__ == '__main__':
    main()
