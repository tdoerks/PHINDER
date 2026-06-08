process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    container = 'docker://quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0'

    input:
    path('*')

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data
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
