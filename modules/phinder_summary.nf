process PHINDER_SUMMARY {
    publishDir "${params.outdir}/summary", mode: 'copy'
    container = 'quay.io/biocontainers/pandas:1.5.2'

    input:
    val(ready)  // Signal that all analyses are complete

    output:
    path "phinder_summary.html", emit: html
    path "phinder_summary.tsv", emit: tsv
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3

    import sys
    sys.path.insert(0, '${projectDir}/bin')

    from generate_phinder_summary import collect_sample_data, generate_html_report, generate_tsv_report

    print("=" * 60)
    print("PHINDER Summary Report Generation")
    print("=" * 60)
    print()

    print("Collecting sample data from: ${params.outdir}")
    samples = collect_sample_data('${params.outdir}')
    print(f"  Found {len(samples)} samples")
    print()

    print("Generating HTML report...")
    generate_html_report(samples, 'phinder_summary.html')
    print("  ✓ phinder_summary.html")
    print()

    print("Generating TSV report...")
    generate_tsv_report(samples, 'phinder_summary.tsv')
    print("  ✓ phinder_summary.tsv")
    print()

    print("=" * 60)
    print("✅ PHINDER summary complete!")
    print("=" * 60)

    # Versions
    with open('versions.yml', 'w') as f:
        f.write('"PHINDER_SUMMARY":\\n')
        f.write('    phinder_summary: "1.0.0"\\n')
    """
}
