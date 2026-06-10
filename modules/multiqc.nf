process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    container = 'docker://quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0'

    input:
    path('*')

    output:
    // MultiQC prefixes outputs with the sanitized --title (e.g.
    // "PHINDER-Report_multiqc_report.html" / "..._multiqc_report_data"), so
    // match with globs rather than fixed names — otherwise Nextflow reports a
    // missing output and fails the process even though multiqc exits 0.
    path "*multiqc_report.html", emit: report
    path "*_data", emit: data
    path "versions.yml", emit: versions

    script:
    def config_opt = params.multiqc_config ? "--config ${params.multiqc_config}" : ""
    def title_opt = params.multiqc_title ? "--title \"${params.multiqc_title}\"" : ""
    """
    multiqc . \\
        ${config_opt} \\
        ${title_opt} \\
        --force \\
        --interactive

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //g')
    END_VERSIONS
    """
}
