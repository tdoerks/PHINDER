#!/usr/bin/env python3
"""
PHINDER Summary Report Generator
Generates tabbed HTML dashboard and TSV reports from PHINDER pipeline outputs
"""

import os
import sys
import glob
import argparse
from pathlib import Path
from datetime import datetime
import csv


# ─── Parsers ──────────────────────────────────────────────────────────────────

def parse_checkv_quality(checkv_dir):
    quality_file = Path(checkv_dir) / "quality_summary.tsv"
    if not quality_file.exists():
        return {}
    with open(quality_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            return {
                'completeness': row.get('completeness', 'N/A'),
                'completeness_method': row.get('completeness_method', 'N/A'),
                'contamination': row.get('contamination', '0'),
                'checkv_quality': row.get('checkv_quality', 'Not determined'),
                'miuvig_quality': row.get('miuvig_quality', 'N/A'),
                'warnings': row.get('warnings', 'None')
            }
    return {}


def parse_quast_report(quast_dir):
    report_file = Path(quast_dir) / "report.tsv"
    if not report_file.exists():
        return {}
    metrics = {}
    with open(report_file) as f:
        for line in f:
            if '\t' in line:
                key, value = line.strip().split('\t', 1)
                metrics[key] = value
    return {
        'total_length': metrics.get('Total length', 'N/A'),
        'num_contigs': metrics.get('# contigs', 'N/A'),
        'largest_contig': metrics.get('Largest contig', 'N/A'),
        'n50': metrics.get('N50', 'N/A'),
        'gc_percent': metrics.get('GC (%)', 'N/A'),
        'n_per_100kb': metrics.get("# N's per 100 kbp", '0')
    }


def parse_pharokka_results(pharokka_dir):
    """Parse Pharokka annotation results from _cds_functions.tsv category summary"""
    functions_files = list(Path(pharokka_dir).glob("*_cds_functions.tsv"))
    if not functions_files:
        return {'total_cds': 0, 'annotated': 0, 'unknown': 0, 'annotation_rate': '0%',
                'trnas': 0, 'crisprs': 0, 'amr_genes': 0, 'virulence_factors': 0,
                'categories': {}}

    total = 0; unknown = 0; trnas = 0; crisprs = 0; amr_genes = 0; virulence_factors = 0
    categories = {}

    seen_header = False
    with open(functions_files[0]) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            desc, count_str = parts[0], parts[1]
            if desc == 'Description':
                if seen_header:
                    break  # file is sometimes duplicated — stop at second header
                seen_header = True
                continue
            try:
                count = int(count_str)
            except ValueError:
                continue
            if desc == 'CDS':
                total = count
            elif desc == 'unknown function':
                unknown = count
            elif desc == 'tRNAs':
                trnas = count
            elif desc == 'CRISPRs':
                crisprs = count
            elif desc == 'CARD_AMR_Genes':
                amr_genes = count
            elif desc == 'VFDB_Virulence_Factors':
                virulence_factors = count
            elif desc not in ('tmRNAs',):
                categories[desc] = count

    annotated = total - unknown
    return {
        'total_cds': total,
        'annotated': annotated,
        'unknown': unknown,
        'annotation_rate': f"{(annotated/total*100):.1f}%" if total > 0 else "0%",
        'trnas': trnas,
        'crisprs': crisprs,
        'amr_genes': amr_genes,
        'virulence_factors': virulence_factors,
        'categories': categories
    }


def parse_amrfinderplus_results(amrfinder_dir, sample_id):
    """Parse AMRFinder Plus output TSV (--plus mode includes AMR, VIRULENCE, STRESS)"""
    report_file = Path(amrfinder_dir) / f"{sample_id}_amrfinder.tsv"
    if not report_file.exists():
        files = list(Path(amrfinder_dir).glob(f"*_amrfinder.tsv"))
        if not files:
            return {'total': 0, 'amr': 0, 'virulence': 0, 'stress': 0, 'genes': []}
        report_file = files[0]

    total = 0; amr = 0; virulence = 0; stress = 0; genes = []
    with open(report_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            total += 1
            etype = row.get('Element type', '').upper()
            symbol = row.get('Gene symbol', '')
            if etype == 'AMR':
                amr += 1
            elif etype == 'VIRULENCE':
                virulence += 1
            elif etype == 'STRESS':
                stress += 1
            if symbol:
                genes.append(f"{symbol} ({etype})")
    return {'total': total, 'amr': amr, 'virulence': virulence,
            'stress': stress, 'genes': genes}


def parse_genomad_results(genomad_dir, sample_id):
    """Parse geNomad virus_summary.tsv for taxonomy, virus score, topology"""
    # publishDir = genomad/{sample_id}/, geNomad output dir = {sample_id}_genomad/
    summary_file = Path(genomad_dir) / sample_id / f"{sample_id}_genomad" / f"{sample_id}_summary" / f"{sample_id}_virus_summary.tsv"
    if not summary_file.exists():
        # Fallback: glob in case directory structure differs
        hits = list(Path(genomad_dir).glob(f"**/{sample_id}_virus_summary.tsv"))
        if not hits:
            return {'virus_score': None, 'taxonomy': '', 'topology': 'N/A',
                    'hallmarks': 0, 'family': None, 'genus': None}
        summary_file = hits[0]

    best = None
    with open(summary_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            try:
                score = float(row.get('virus_score', 0) or 0)
            except ValueError:
                score = 0
            if best is None or score > float(best.get('virus_score', 0) or 0):
                best = row

    if not best:
        return {'virus_score': None, 'taxonomy': '', 'topology': 'N/A',
                'hallmarks': 0, 'family': None, 'genus': None}

    taxonomy = best.get('taxonomy', '') or ''
    parts = [p.strip() for p in taxonomy.split(';')]
    # geNomad taxonomy: Viruses;kingdom;phylum;class;order;family;genus;species
    family = parts[6] if len(parts) > 6 and parts[6] else None
    genus  = parts[7] if len(parts) > 7 and parts[7] else None

    return {
        'virus_score': best.get('virus_score', 'N/A'),
        'taxonomy': taxonomy,
        'topology': best.get('topology', 'N/A'),
        'hallmarks': int(best.get('n_hallmarks', 0) or 0),
        'family': family,
        'genus': genus
    }


def parse_phageterm_results(phageterm_dir, sample_id):
    """Parse PhageTerm packaging strategy output (reads mode only)"""
    strategy_file = Path(phageterm_dir) / sample_id / f"{sample_id}_phageterm_strategy.txt"
    if not strategy_file.exists():
        alt = list(Path(phageterm_dir).glob(f"{sample_id}/*_phageterm_strategy.txt"))
        if not alt:
            return {'strategy': None}  # None = module not run (assembly mode)
        strategy_file = alt[0]
    with open(strategy_file) as f:
        content = f.read().strip()
    strategy = content.split('\n')[0].strip() if content else 'Unknown'
    return {'strategy': strategy if strategy else 'Unknown'}


def parse_bacphlip_results(bacphlip_dir, sample_id):
    pred_file = Path(bacphlip_dir) / f"{sample_id}.bacphlip"
    if not pred_file.exists():
        pred_files = list(Path(bacphlip_dir).glob("*.bacphlip"))
        if not pred_files:
            return {'lifestyle': 'Not determined', 'virulent_prob': 'N/A', 'temperate_prob': 'N/A'}
        pred_file = pred_files[0]
    with open(pred_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            virulent = float(row.get('Virulent', 0))
            temperate = float(row.get('Temperate', 0))
            lifestyle = 'Virulent (lytic)' if virulent >= 0.5 else 'Temperate (lysogenic)'
            return {
                'lifestyle': lifestyle,
                'virulent_prob': f"{virulent:.3f}",
                'temperate_prob': f"{temperate:.3f}",
                'virulent_pct': virulent * 100
            }
    return {'lifestyle': 'Not determined', 'virulent_prob': 'N/A', 'temperate_prob': 'N/A', 'virulent_pct': 0}


def parse_vibrant_results(vibrant_dir):
    quality_files = list(Path(vibrant_dir).glob("**/VIBRANT_genome_quality_*.tsv"))
    if not quality_files:
        return {'lifestyle': 'Not determined', 'confidence': 'N/A'}
    with open(quality_files[0]) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            lifestyle = row.get('type', 'Not determined')
            if lifestyle == 'lytic':
                lifestyle = 'Lytic'
            elif lifestyle in ('lysogenic', 'temperate'):
                lifestyle = 'Temperate (lysogenic)'
            return {'lifestyle': lifestyle, 'confidence': row.get('score', 'N/A')}
    return {'lifestyle': 'Not determined', 'confidence': 'N/A'}


# ─── Data collection ──────────────────────────────────────────────────────────

def collect_sample_data(outdir):
    samples = {}
    quast_base = Path(outdir) / "quast"
    sample_ids = set()
    if quast_base.exists():
        for quast_dir in quast_base.glob("*_quast"):
            sample_ids.add(quast_dir.name.replace('_quast', ''))
    if not sample_ids:
        assembly_dir = Path(outdir) / "assemblies"
        if assembly_dir.exists():
            for f in assembly_dir.glob("*_assembly.fasta"):
                sample_ids.add(f.stem.replace('_assembly', ''))

    for sample_id in sorted(sample_ids):
        samples[sample_id] = {
            'sample_id': sample_id,
            'checkv': {}, 'quast': {}, 'pharokka': {}, 'vibrant': {}, 'bacphlip': {},
            'amrfinderplus': {}, 'genomad': {}, 'phageterm': {}
        }
        checkv_dir = Path(outdir) / "checkv" / f"{sample_id}_checkv"
        if checkv_dir.exists():
            samples[sample_id]['checkv'] = parse_checkv_quality(checkv_dir)
        quast_dir = Path(outdir) / "quast" / f"{sample_id}_quast"
        if quast_dir.exists():
            samples[sample_id]['quast'] = parse_quast_report(quast_dir)
        pharokka_dir = Path(outdir) / "pharokka" / f"{sample_id}_pharokka"
        if pharokka_dir.exists():
            samples[sample_id]['pharokka'] = parse_pharokka_results(pharokka_dir)
        vibrant_dir = Path(outdir) / "vibrant" / f"{sample_id}_vibrant"
        if vibrant_dir.exists():
            samples[sample_id]['vibrant'] = parse_vibrant_results(vibrant_dir)
        bacphlip_dir = Path(outdir) / "bacphlip"
        if bacphlip_dir.exists():
            samples[sample_id]['bacphlip'] = parse_bacphlip_results(bacphlip_dir, sample_id)
        amrfinder_dir = Path(outdir) / "amrfinderplus"
        if amrfinder_dir.exists():
            samples[sample_id]['amrfinderplus'] = parse_amrfinderplus_results(amrfinder_dir, sample_id)
        genomad_dir = Path(outdir) / "genomad"
        if genomad_dir.exists():
            samples[sample_id]['genomad'] = parse_genomad_results(genomad_dir, sample_id)
        phageterm_dir = Path(outdir) / "phageterm"
        if phageterm_dir.exists():
            samples[sample_id]['phageterm'] = parse_phageterm_results(phageterm_dir, sample_id)

    return samples


# ─── HTML row builders ────────────────────────────────────────────────────────

def _qbadge(quality):
    q = quality.lower().replace('-', '').replace(' ', '')
    cls_map = {'complete': 'complete', 'highquality': 'high', 'mediumquality': 'medium',
               'lowquality': 'low', 'notdetermined': 'nd', 'fragment': 'fragment'}
    cls = cls_map.get(q, 'nd')
    return f'<span class="badge badge-{cls}">{quality}</span>'


def _lbadge(lifestyle):
    ls = lifestyle.lower()
    if 'lytic' in ls or 'virulent' in ls:
        return f'<span class="badge badge-lytic">{lifestyle}</span>'
    elif 'temperate' in ls or 'lysogenic' in ls:
        return f'<span class="badge badge-temperate">{lifestyle}</span>'
    else:
        return f'<span class="badge badge-nd">{lifestyle}</span>'


def _probbar(virulent_pct):
    lytic_w = max(0, min(100, virulent_pct))
    temp_w = 100 - lytic_w
    return (f'<div class="probbar">'
            f'<div class="lytic-fill" style="width:{lytic_w:.1f}%"></div>'
            f'<div class="temp-fill" style="width:{temp_w:.1f}%"></div>'
            f'</div>')


def _catpills(categories):
    pills = []
    for cat, cnt in sorted(categories.items()):
        if cnt > 0 and cat not in ('CDS', 'unknown function'):
            label = cat.replace('_', ' ').title()
            pills.append(f'<span class="catpill">{label} <span>{cnt}</span></span>')
    return ''.join(pills) if pills else '<span style="color:var(--muted)">—</span>'


def _fmt_size(val):
    try:
        n = int(val)
        if n >= 1_000_000:
            return f"{n/1_000_000:.2f} Mbp"
        elif n >= 1_000:
            return f"{n/1_000:.1f} kbp"
        return str(n)
    except (ValueError, TypeError):
        return str(val) if val else 'N/A'


def build_overview_rows(samples):
    rows = []
    for sid, d in samples.items():
        size = _fmt_size(d['quast'].get('total_length', ''))
        gc = d['quast'].get('gc_percent', 'N/A')
        gc_str = f"{gc}%" if gc != 'N/A' else 'N/A'
        quality = d['checkv'].get('checkv_quality', 'Not determined')
        comp = d['checkv'].get('completeness', 'N/A')
        comp_str = f"{comp}%" if comp not in ('N/A', '', 'NA') else 'N/A'
        genes = d['pharokka'].get('total_cds', 0)
        ann_rate = d['pharokka'].get('annotation_rate', 'N/A')
        vib = d['vibrant'].get('lifestyle', 'Not determined')
        bp = d['bacphlip'].get('lifestyle', 'Not determined')
        amr = d['pharokka'].get('amr_genes', 0)
        amr_badge = (f'<span style="color:var(--bad);font-weight:700">{amr}</span>'
                     if amr > 0 else f'<span style="color:var(--good)">{amr}</span>')
        pkg = d['phageterm'].get('strategy', None)
        pkg_cell = (f'<span style="color:var(--muted)">—</span>' if pkg is None
                    else f'<span class="badge badge-nd">{pkg}</span>')
        gd = d.get('genomad', {})
        family = gd.get('family')
        genus  = gd.get('genus')
        if family or genus:
            tax_cell = (f'<span style="color:var(--ink)">{family or "?"}</span>'
                        + (f'<br><small style="color:var(--muted)">{genus}</small>' if genus else ''))
        elif gd.get('virus_score') is not None:
            tax_cell = '<span style="color:var(--muted)">Unclassified</span>'
        else:
            tax_cell = '<span style="color:var(--muted)">—</span>'
        rows.append(f"""<tr>
          <td><strong>{sid}</strong></td>
          <td>{size}</td>
          <td>{gc_str}</td>
          <td>{_qbadge(quality)}</td>
          <td>{comp_str}</td>
          <td>{genes}</td>
          <td>{ann_rate}</td>
          <td>{_lbadge(vib)}</td>
          <td>{_lbadge(bp)}</td>
          <td>{pkg_cell}</td>
          <td>{tax_cell}</td>
          <td>{amr_badge}</td>
        </tr>""")
    return '\n'.join(rows)


def build_lifestyle_rows(samples):
    rows = []
    for sid, d in samples.items():
        vib = d['vibrant'].get('lifestyle', 'Not determined')
        bp = d['bacphlip'].get('lifestyle', 'Not determined')
        vp = d['bacphlip'].get('virulent_prob', 'N/A')
        tp = d['bacphlip'].get('temperate_prob', 'N/A')
        pct = d['bacphlip'].get('virulent_pct', 0)

        # Normalize for comparison
        vib_lytic = 'lytic' in vib.lower() or 'virulent' in vib.lower()
        bp_lytic = 'lytic' in bp.lower() or 'virulent' in bp.lower()
        agree = vib_lytic == bp_lytic
        agree_cell = ('<span class="agree">✓ Agree</span>' if agree
                      else '<span class="disagree">⚠ Disagree</span>')

        bar = _probbar(pct) if pct != 0 or vp != 'N/A' else '—'
        rows.append(f"""<tr>
          <td><strong>{sid}</strong></td>
          <td>{_lbadge(vib)}</td>
          <td>{_lbadge(bp)}</td>
          <td>{bar} {vp}</td>
          <td>{tp}</td>
          <td>{agree_cell}</td>
          <td><span style="color:var(--muted)">—</span></td>
        </tr>""")
    return '\n'.join(rows)


def build_annotation_rows(samples):
    rows = []
    for sid, d in samples.items():
        p = d['pharokka']
        cats = _catpills(p.get('categories', {}))
        rows.append(f"""<tr>
          <td><strong>{sid}</strong></td>
          <td>{p.get('total_cds', 0)}</td>
          <td>{p.get('annotated', 0)}</td>
          <td>{p.get('unknown', 0)}</td>
          <td>{p.get('annotation_rate', 'N/A')}</td>
          <td>{p.get('trnas', 0)}</td>
          <td>{p.get('crisprs', 0)}</td>
          <td>{cats}</td>
        </tr>""")
    return '\n'.join(rows)


def build_quality_rows(samples):
    rows = []
    for sid, d in samples.items():
        c = d['checkv']
        quality = c.get('checkv_quality', 'Not determined')
        comp_raw = c.get('completeness', 'N/A')
        try:
            comp_f = float(comp_raw)
            bar_color = 'var(--good)' if comp_f >= 90 else ('var(--warn)' if comp_f >= 50 else 'var(--bad)')
            comp_cell = (f'{comp_f:.1f}%'
                         f'<div class="compbar"><div class="fill" style="width:{min(comp_f,100):.1f}%;background:{bar_color}"></div></div>')
        except (ValueError, TypeError):
            comp_cell = comp_raw if comp_raw else 'N/A'
        cont = c.get('contamination', '0')
        cont_str = f'<span style="color:var(--bad)">{cont}%</span>' if cont not in ('0', '0.0', 'N/A', '') and float(cont or 0) > 0 else (cont or '0') + '%'
        warn = c.get('warnings', 'None') or 'None'
        warn_str = f'<span style="color:var(--warn)">{warn}</span>' if warn and warn != 'None' else '<span style="color:var(--muted)">None</span>'
        rows.append(f"""<tr>
          <td><strong>{sid}</strong></td>
          <td>{_qbadge(quality)}</td>
          <td>{comp_cell}</td>
          <td>{c.get('completeness_method', 'N/A')}</td>
          <td>{cont_str}</td>
          <td>{c.get('miuvig_quality', 'N/A')}</td>
          <td>{warn_str}</td>
        </tr>""")
    return '\n'.join(rows)


def build_taxonomy_rows(samples):
    rows = []
    for sid, d in samples.items():
        gd = d.get('genomad', {})
        vs = gd.get('virus_score')
        if vs is None:
            rows.append(f"""<tr>
              <td><strong>{sid}</strong></td>
              <td colspan="6" style="color:var(--muted);text-align:center">geNomad not run</td>
            </tr>""")
            continue
        try:
            vs_f = float(vs)
            score_color = 'var(--good)' if vs_f >= 0.9 else ('var(--warn)' if vs_f >= 0.7 else 'var(--bad)')
            score_cell = f'<span style="color:{score_color};font-weight:700">{vs_f:.3f}</span>'
        except (ValueError, TypeError):
            score_cell = str(vs)

        taxonomy = gd.get('taxonomy', '')
        parts = [p.strip() for p in taxonomy.split(';')] if taxonomy else []

        def _tpart(idx, fallback='—'):
            return parts[idx] if len(parts) > idx and parts[idx] else f'<span style="color:var(--muted)">{fallback}</span>'

        rows.append(f"""<tr>
          <td><strong>{sid}</strong></td>
          <td>{score_cell}</td>
          <td>{gd.get('topology', 'N/A')}</td>
          <td>{gd.get('hallmarks', 0)}</td>
          <td>{_tpart(4)}</td>
          <td>{_tpart(5)}</td>
          <td>{_tpart(6)}</td>
          <td>{_tpart(7)}</td>
        </tr>""")
    return '\n'.join(rows)


def build_safety_rows(samples):
    rows = []
    for sid, d in samples.items():
        p = d['pharokka']
        af = d.get('amrfinderplus', {})
        amr = p.get('amr_genes', 0)
        vf = p.get('virulence_factors', 0)
        integ = p.get('categories', {}).get('integration/excision', 0)
        af_total = af.get('total', None)  # None = not run
        af_amr = af.get('amr', 0)
        af_vir = af.get('virulence', 0)
        af_stress = af.get('stress', 0)
        af_genes = af.get('genes', [])

        safe = amr == 0 and vf == 0 and (af_total is None or af_amr == 0)
        amr_cell = (f'<span style="color:var(--bad);font-weight:700">{amr}</span>'
                    if amr > 0 else f'<span style="color:var(--good)">{amr}</span>')
        vf_cell = (f'<span style="color:var(--bad);font-weight:700">{vf}</span>'
                   if vf > 0 else f'<span style="color:var(--good)">{vf}</span>')
        integ_cell = (f'<span style="color:var(--warn)">{integ}</span>'
                      if integ > 0 else f'<span style="color:var(--muted)">{integ}</span>')
        if af_total is None:
            af_cell = '<span style="color:var(--muted)">—</span>'
        elif af_total == 0:
            af_cell = '<span style="color:var(--good)">0</span>'
        else:
            gene_list = ', '.join(af_genes[:5]) + ('…' if len(af_genes) > 5 else '')
            af_cell = (f'<span style="color:var(--bad);font-weight:700">'
                       f'AMR:{af_amr} VIR:{af_vir} STRESS:{af_stress}</span>'
                       f'<br><small style="color:var(--muted)">{gene_list}</small>')
        safety_badge = ('<span class="badge badge-safe">✓ Safe</span>' if safe
                        else '<span class="badge badge-warn">⚠ Review Required</span>')
        rows.append(f"""<tr>
          <td><strong>{sid}</strong></td>
          <td>{amr_cell}</td>
          <td>{vf_cell}</td>
          <td>{af_cell}</td>
          <td>{integ_cell}</td>
          <td>{safety_badge}</td>
        </tr>""")
    return '\n'.join(rows)


# ─── Report generators ────────────────────────────────────────────────────────

def generate_html_report(samples, output_file):
    template_path = Path(__file__).parent / 'phinder_dashboard_template.html'
    if not template_path.exists():
        # Fallback: generate simple report without template
        _generate_simple_html(samples, output_file)
        return

    with open(template_path) as f:
        html = f.read()

    # Summary stats
    n_phages = len(samples)
    n_hq = sum(1 for s in samples.values()
                if s['checkv'].get('checkv_quality') in ('Complete', 'High-quality'))
    n_genes = sum(s['pharokka'].get('total_cds', 0) for s in samples.values())
    n_lytic = sum(1 for s in samples.values()
                  if ('lytic' in s['vibrant'].get('lifestyle', '').lower()
                      or 'virulent' in s['bacphlip'].get('lifestyle', '').lower()))
    n_amr = sum(s['pharokka'].get('amr_genes', 0) for s in samples.values())

    html = (html
        .replace('__TIMESTAMP__', datetime.now().strftime('%Y-%m-%d %H:%M'))
        .replace('__N_PHAGES__', str(n_phages))
        .replace('__N_HQ__', str(n_hq))
        .replace('__N_GENES__', str(n_genes))
        .replace('__N_LYTIC__', str(n_lytic))
        .replace('__N_AMR__', str(n_amr))
        .replace('__OVERVIEW_ROWS__', build_overview_rows(samples))
        .replace('__LIFESTYLE_ROWS__', build_lifestyle_rows(samples))
        .replace('__ANNOTATION_ROWS__', build_annotation_rows(samples))
        .replace('__QUALITY_ROWS__', build_quality_rows(samples))
        .replace('__TAXONOMY_ROWS__', build_taxonomy_rows(samples))
        .replace('__SAFETY_ROWS__', build_safety_rows(samples))
    )

    with open(output_file, 'w') as f:
        f.write(html)


def _generate_simple_html(samples, output_file):
    """Minimal fallback if template not found"""
    lines = ['<!DOCTYPE html><html><head><meta charset="UTF-8"><title>PHINDER</title></head><body>',
             f'<h1>PHINDER Summary — {datetime.now().strftime("%Y-%m-%d %H:%M")}</h1>',
             f'<p>{len(samples)} phages analyzed</p><ul>']
    for sid, d in samples.items():
        lines.append(f'<li><strong>{sid}</strong>: {d["checkv"].get("checkv_quality","?")} | '
                     f'{d["pharokka"].get("total_cds",0)} genes | '
                     f'{d["vibrant"].get("lifestyle","?")}</li>')
    lines += ['</ul></body></html>']
    with open(output_file, 'w') as f:
        f.write('\n'.join(lines))


def generate_tsv_report(samples, output_file):
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            'sample_id', 'genome_size', 'num_contigs', 'n50', 'gc_percent',
            'checkv_quality', 'completeness', 'contamination', 'miuvig_quality',
            'total_genes', 'annotated_genes', 'unknown_genes', 'annotation_rate',
            'trnas', 'crisprs', 'amr_genes', 'virulence_factors',
            'vibrant_lifestyle', 'bacphlip_lifestyle', 'bacphlip_virulent_prob',
            'phageterm_packaging_strategy'
        ])
        for sample_id, data in samples.items():
            writer.writerow([
                sample_id,
                data['quast'].get('total_length', 'N/A'),
                data['quast'].get('num_contigs', 'N/A'),
                data['quast'].get('n50', 'N/A'),
                data['quast'].get('gc_percent', 'N/A'),
                data['checkv'].get('checkv_quality', 'N/A'),
                data['checkv'].get('completeness', 'N/A'),
                data['checkv'].get('contamination', 'N/A'),
                data['checkv'].get('miuvig_quality', 'N/A'),
                data['pharokka'].get('total_cds', 0),
                data['pharokka'].get('annotated', 0),
                data['pharokka'].get('unknown', 0),
                data['pharokka'].get('annotation_rate', '0%'),
                data['pharokka'].get('trnas', 0),
                data['pharokka'].get('crisprs', 0),
                data['pharokka'].get('amr_genes', 0),
                data['pharokka'].get('virulence_factors', 0),
                data['vibrant'].get('lifestyle', 'N/A'),
                data['bacphlip'].get('lifestyle', 'N/A'),
                data['bacphlip'].get('virulent_prob', 'N/A'),
                data['phageterm'].get('strategy', 'N/A')
            ])


# ─── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Generate PHINDER summary report')
    parser.add_argument('--outdir', required=True, help='PHINDER output directory')
    parser.add_argument('--output-html', default='phinder_summary.html')
    parser.add_argument('--output-tsv', default='phinder_summary.tsv')
    args = parser.parse_args()

    print("=" * 60)
    print("PHINDER Summary Report Generator")
    print("=" * 60)

    samples = collect_sample_data(args.outdir)
    print(f"Found {len(samples)} samples")

    generate_html_report(samples, args.output_html)
    print(f"✓ {args.output_html}")

    generate_tsv_report(samples, args.output_tsv)
    print(f"✓ {args.output_tsv}")

    print("=" * 60)
    print("✅ PHINDER summary complete!")

if __name__ == '__main__':
    main()
