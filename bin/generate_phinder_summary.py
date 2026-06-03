#!/usr/bin/env python3
"""
PHINDER Summary Report Generator
Generates comprehensive HTML and TSV reports from PHINDER pipeline outputs
"""

import os
import sys
import json
import glob
import argparse
from pathlib import Path
from datetime import datetime
import csv

def parse_checkv_quality(checkv_dir):
    """Parse CheckV quality summary"""
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
    """Parse QUAST assembly metrics"""
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
    """Parse Pharokka annotation results"""
    functions_file = Path(pharokka_dir) / f"{pharokka_dir.name.replace('_pharokka', '')}_cds_functions.tsv"

    if not functions_file.exists():
        # Try alternative naming
        functions_files = list(Path(pharokka_dir).glob("*_cds_functions.tsv"))
        if functions_files:
            functions_file = functions_files[0]
        else:
            return {'total_cds': 0, 'annotated': 0, 'hypothetical': 0}

    total = 0
    annotated = 0
    hypothetical = 0

    with open(functions_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            total += 1
            annotation = row.get('annot', '').lower()
            if 'hypothetical' in annotation or not annotation:
                hypothetical += 1
            else:
                annotated += 1

    return {
        'total_cds': total,
        'annotated': annotated,
        'hypothetical': hypothetical,
        'annotation_rate': f"{(annotated/total*100):.1f}%" if total > 0 else "0%"
    }

def parse_vibrant_results(vibrant_dir):
    """Parse VIBRANT lifestyle prediction"""
    quality_file = list(Path(vibrant_dir).glob("**/VIBRANT_genome_quality*.tsv"))

    if not quality_file:
        return {'lifestyle': 'Not determined', 'confidence': 'N/A'}

    with open(quality_file[0]) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            return {
                'lifestyle': row.get('type', 'Not determined'),
                'confidence': row.get('score', 'N/A')
            }

    return {'lifestyle': 'Not determined', 'confidence': 'N/A'}

def collect_sample_data(outdir):
    """Collect data for all samples"""
    samples = {}

    # Find all assemblies
    assembly_dir = Path(outdir) / "assemblies"
    if assembly_dir.exists():
        for assembly_file in assembly_dir.glob("*_assembly.fasta"):
            sample_id = assembly_file.stem.replace('_assembly', '')

            samples[sample_id] = {
                'sample_id': sample_id,
                'assembly_file': str(assembly_file),
                'checkv': {},
                'quast': {},
                'pharokka': {},
                'vibrant': {}
            }

            # CheckV
            checkv_dir = Path(outdir) / "checkv" / f"{sample_id}_checkv"
            if checkv_dir.exists():
                samples[sample_id]['checkv'] = parse_checkv_quality(checkv_dir)

            # QUAST
            quast_dir = Path(outdir) / "quast" / f"{sample_id}_quast"
            if quast_dir.exists():
                samples[sample_id]['quast'] = parse_quast_report(quast_dir)

            # Pharokka
            pharokka_dir = Path(outdir) / "pharokka" / f"{sample_id}_pharokka"
            if pharokka_dir.exists():
                samples[sample_id]['pharokka'] = parse_pharokka_results(pharokka_dir)

            # VIBRANT
            vibrant_dir = Path(outdir) / "vibrant" / f"{sample_id}_vibrant"
            if vibrant_dir.exists():
                samples[sample_id]['vibrant'] = parse_vibrant_results(vibrant_dir)

    return samples

def generate_html_report(samples, output_file):
    """Generate interactive HTML report"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PHINDER Summary Report</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
            color: #333;
        }}

        .container {{
            max-width: 1400px;
            margin: 0 auto;
        }}

        .header {{
            background: white;
            padding: 30px;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.2);
            margin-bottom: 30px;
            text-align: center;
        }}

        .header h1 {{
            font-size: 2.5em;
            color: #667eea;
            margin-bottom: 10px;
        }}

        .header .subtitle {{
            color: #666;
            font-size: 1.1em;
        }}

        .header .timestamp {{
            margin-top: 15px;
            color: #999;
            font-size: 0.9em;
        }}

        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}

        .stat-card {{
            background: white;
            padding: 25px;
            border-radius: 10px;
            box-shadow: 0 5px 15px rgba(0,0,0,0.1);
            text-align: center;
            transition: transform 0.2s;
        }}

        .stat-card:hover {{
            transform: translateY(-5px);
        }}

        .stat-card .value {{
            font-size: 2.5em;
            font-weight: bold;
            color: #667eea;
            margin-bottom: 10px;
        }}

        .stat-card .label {{
            color: #666;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}

        .sample-card {{
            background: white;
            padding: 25px;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
            margin-bottom: 20px;
        }}

        .sample-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 20px;
            padding-bottom: 15px;
            border-bottom: 2px solid #f0f0f0;
        }}

        .sample-title {{
            font-size: 1.5em;
            color: #667eea;
            font-weight: bold;
        }}

        .quality-badge {{
            padding: 8px 16px;
            border-radius: 20px;
            font-weight: bold;
            font-size: 0.9em;
        }}

        .quality-complete {{
            background: #10b981;
            color: white;
        }}

        .quality-high {{
            background: #3b82f6;
            color: white;
        }}

        .quality-medium {{
            background: #f59e0b;
            color: white;
        }}

        .quality-low {{
            background: #ef4444;
            color: white;
        }}

        .quality-undetermined {{
            background: #6b7280;
            color: white;
        }}

        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 20px;
        }}

        .metric-box {{
            background: #f9fafb;
            padding: 15px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}

        .metric-label {{
            color: #666;
            font-size: 0.85em;
            margin-bottom: 5px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}

        .metric-value {{
            font-size: 1.3em;
            font-weight: bold;
            color: #333;
        }}

        .lifestyle-lytic {{
            color: #ef4444;
        }}

        .lifestyle-temperate {{
            color: #f59e0b;
        }}

        .lifestyle-lysogenic {{
            color: #8b5cf6;
        }}

        footer {{
            background: white;
            padding: 20px;
            border-radius: 10px;
            text-align: center;
            margin-top: 30px;
            color: #666;
        }}

        .section-title {{
            font-size: 1.2em;
            color: #667eea;
            margin-bottom: 15px;
            font-weight: bold;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 PHINDER Summary Report</h1>
            <p class="subtitle">PHage Isolate characterizatioN, Discovery & Evaluation Resource</p>
            <p class="timestamp">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>

        <div class="stats-grid">
            <div class="stat-card">
                <div class="value">{len(samples)}</div>
                <div class="label">Phages Analyzed</div>
            </div>
            <div class="stat-card">
                <div class="value">{sum(1 for s in samples.values() if s['checkv'].get('checkv_quality') == 'Complete')}</div>
                <div class="label">Complete Genomes</div>
            </div>
            <div class="stat-card">
                <div class="value">{sum(s['pharokka'].get('total_cds', 0) for s in samples.values())}</div>
                <div class="label">Total Genes</div>
            </div>
            <div class="stat-card">
                <div class="value">{sum(1 for s in samples.values() if 'lytic' in s['vibrant'].get('lifestyle', '').lower())}</div>
                <div class="label">Lytic Phages</div>
            </div>
        </div>
"""

    # Add sample cards
    for sample_id, data in samples.items():
        quality = data['checkv'].get('checkv_quality', 'Not determined')
        quality_class = f"quality-{quality.lower().replace(' ', '-')}"

        lifestyle = data['vibrant'].get('lifestyle', 'Not determined')
        lifestyle_class = ''
        if 'lytic' in lifestyle.lower():
            lifestyle_class = 'lifestyle-lytic'
        elif 'temperate' in lifestyle.lower() or 'lysogenic' in lifestyle.lower():
            lifestyle_class = 'lifestyle-lysogenic'

        html += f"""
        <div class="sample-card">
            <div class="sample-header">
                <div class="sample-title">📦 {sample_id}</div>
                <div class="quality-badge {quality_class}">{quality}</div>
            </div>

            <div class="section-title">Assembly Metrics</div>
            <div class="metrics-grid">
                <div class="metric-box">
                    <div class="metric-label">Genome Size</div>
                    <div class="metric-value">{data['quast'].get('total_length', 'N/A')}</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">Contigs</div>
                    <div class="metric-value">{data['quast'].get('num_contigs', 'N/A')}</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">N50</div>
                    <div class="metric-value">{data['quast'].get('n50', 'N/A')}</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">GC Content</div>
                    <div class="metric-value">{data['quast'].get('gc_percent', 'N/A')}%</div>
                </div>
            </div>

            <div class="section-title">Quality Assessment (CheckV)</div>
            <div class="metrics-grid">
                <div class="metric-box">
                    <div class="metric-label">Completeness</div>
                    <div class="metric-value">{data['checkv'].get('completeness', 'N/A')}</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">Contamination</div>
                    <div class="metric-value">{data['checkv'].get('contamination', 'N/A')}</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">MIUViG Quality</div>
                    <div class="metric-value">{data['checkv'].get('miuvig_quality', 'N/A')}</div>
                </div>
            </div>

            <div class="section-title">Annotation (Pharokka)</div>
            <div class="metrics-grid">
                <div class="metric-box">
                    <div class="metric-label">Total Genes</div>
                    <div class="metric-value">{data['pharokka'].get('total_cds', 0)}</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">Annotated</div>
                    <div class="metric-value">{data['pharokka'].get('annotated', 0)}</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">Hypothetical</div>
                    <div class="metric-value">{data['pharokka'].get('hypothetical', 0)}</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">Annotation Rate</div>
                    <div class="metric-value">{data['pharokka'].get('annotation_rate', '0%')}</div>
                </div>
            </div>

            <div class="section-title">Lifestyle Prediction (VIBRANT)</div>
            <div class="metrics-grid">
                <div class="metric-box">
                    <div class="metric-label">Predicted Lifestyle</div>
                    <div class="metric-value {lifestyle_class}">{lifestyle}</div>
                </div>
            </div>
        </div>
"""

    html += """
        <footer>
            <p><strong>PHINDER</strong> - PHage Isolate characterizatioN, Discovery & Evaluation Resource</p>
            <p>Generated by PHINDER pipeline | <a href="https://github.com/tdoerks/PHINDER">github.com/tdoerks/PHINDER</a></p>
        </footer>
    </div>
</body>
</html>
"""

    with open(output_file, 'w') as f:
        f.write(html)

def generate_tsv_report(samples, output_file):
    """Generate TSV summary"""

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')

        # Header
        writer.writerow([
            'sample_id',
            'genome_size',
            'num_contigs',
            'n50',
            'gc_percent',
            'checkv_quality',
            'completeness',
            'contamination',
            'total_genes',
            'annotated_genes',
            'annotation_rate',
            'lifestyle'
        ])

        # Data rows
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
                data['pharokka'].get('total_cds', 0),
                data['pharokka'].get('annotated', 0),
                data['pharokka'].get('annotation_rate', '0%'),
                data['vibrant'].get('lifestyle', 'Not determined')
            ])

def main():
    parser = argparse.ArgumentParser(description='Generate PHINDER summary report')
    parser.add_argument('--outdir', required=True, help='PHINDER output directory')
    parser.add_argument('--output-html', default='phinder_summary.html', help='Output HTML file')
    parser.add_argument('--output-tsv', default='phinder_summary.tsv', help='Output TSV file')

    args = parser.parse_args()

    print("=" * 60)
    print("PHINDER Summary Report Generator")
    print("=" * 60)
    print()

    print("Collecting sample data...")
    samples = collect_sample_data(args.outdir)
    print(f"  Found {len(samples)} samples")
    print()

    print("Generating HTML report...")
    generate_html_report(samples, args.output_html)
    print(f"  ✓ {args.output_html}")
    print()

    print("Generating TSV report...")
    generate_tsv_report(samples, args.output_tsv)
    print(f"  ✓ {args.output_tsv}")
    print()

    print("=" * 60)
    print("✅ PHINDER summary complete!")
    print("=" * 60)

if __name__ == '__main__':
    main()
