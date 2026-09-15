process IPHOP {
    tag "$sample_id"
    publishDir "${params.outdir}/iphop/${sample_id}", mode: 'copy'
    // iPHoP: Integrated Phage-Host Prediction
    // Container: quay.io/biocontainers/iphop:1.3.3--pyhdfd78af_0
    // Database: ~400GB download required — see docs/iphop_setup.md
    // Setup: iphop download --db_dir /path/to/iphop_db --no_prompt
    container = 'quay.io/biocontainers/iphop:1.3.3--pyhdfd78af_0'

    input:
    tuple val(sample_id), path(assembly)

    output:
    path "${sample_id}_iphop/", emit: results
    path "${sample_id}_iphop/Host_prediction_to_genome_m*.csv", emit: predictions, optional: true
    path "versions.yml", emit: versions

    script:
    if (!params.iphop_db) error "iPHoP requires params.iphop_db — set in nextflow.config or via --iphop_db"
    """
    iphop predict \\
        --fa_file ${assembly} \\
        --db_dir ${params.iphop_db} \\
        --out_dir ${sample_id}_iphop \\
        --min_score ${params.iphop_min_score ?: 75} \\
        --num_threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iphop: \$(iphop --version 2>&1 | head -1)
    END_VERSIONS
    """
}
