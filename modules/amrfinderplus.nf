process AMRFINDERPLUS {
    tag "$sample_id"
    publishDir "${params.outdir}/amrfinderplus", mode: 'copy'
    container = 'staphb/ncbi-amrfinderplus:3.12.8'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_amrfinder.tsv"), emit: results
    path "${sample_id}_amrfinder.tsv", emit: report
    path "versions.yml", emit: versions

    script:
    def db_arg = params.amrfinder_db ? "--database ${params.amrfinder_db}" : ""
    """
    amrfinder \\
        --nucleotide ${assembly} \\
        --plus \\
        --name ${sample_id} \\
        --output ${sample_id}_amrfinder.tsv \\
        --threads ${task.cpus} \\
        ${db_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version 2>&1 | head -1)
    END_VERSIONS
    """
}
