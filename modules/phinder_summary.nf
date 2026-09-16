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

    from generate_phinder_summary import collect_sample_data, collect_fastani_data, generate_html_report, generate_tsv_report

    print("=" * 60)
    print("PHINDER Summary Report Generation")
    print("=" * 60)
    print()

    outdir = '${launchDir}/${params.outdir}'
    print(f"Collecting sample data from: {outdir}")
    samples = collect_sample_data(outdir)
    print(f"  Found {len(samples)} samples")
    fastani_data = collect_fastani_data(outdir)
    print()

    run_meta = {
        'pipeline_version': '${workflow.manifest.version}',
        'nextflow_version': '${workflow.nextflow.version}',
        'run_name': '${workflow.runName}',
        'input': '${params.input}',
        'input_mode': '${params.input_mode}',
        'outdir': '${params.outdir}',
        'parameters': {
            'assembler': '${params.assembler}',
            'skip_assembly': '${params.skip_assembly}',
            'skip_fastqc': '${params.skip_fastqc}',
            'skip_fastp': '${params.skip_fastp}',
            'skip_pharokka': '${params.skip_pharokka}',
            'skip_vibrant': '${params.skip_vibrant}',
            'skip_phanotate': '${params.skip_phanotate}',
            'skip_diamond': '${params.skip_diamond}',
            'skip_bacphlip': '${params.skip_bacphlip}',
            'skip_checkv': '${params.skip_checkv}',
            'skip_amrfinderplus': '${params.skip_amrfinderplus}',
            'skip_genomad': '${params.skip_genomad}',
            'skip_fastani': '${params.skip_fastani}',
            'skip_phageterm': '${params.skip_phageterm}',
            'checkv_db': '${params.checkv_db}',
            'pharokka_db': '${params.pharokka_db}',
            'prophage_db': '${params.prophage_db}',
            'amrfinder_db': '${params.amrfinder_db}',
            'genomad_db': '${params.genomad_db}',
        },
    }

    print("Generating HTML report...")
    generate_html_report(samples, fastani_data, 'phinder_summary.html',
                         run_meta=run_meta, outdir=outdir)
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
