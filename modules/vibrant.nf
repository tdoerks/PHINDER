process VIBRANT {
    tag "$sample_id"
    publishDir "${params.outdir}/vibrant", mode: 'copy'
    container = 'docker://staphb/vibrant:4.0.0'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_vibrant"), emit: results
    path "${sample_id}_vibrant/VIBRANT_${sample_id}/VIBRANT_results_${sample_id}/VIBRANT_genome_quality_${sample_id}.tsv", emit: quality, optional: true
    path "${sample_id}_vibrant/VIBRANT_${sample_id}/VIBRANT_phages_${sample_id}/${sample_id}.phages_combined.fna", emit: phages, optional: true
    path "versions.yml", emit: versions

    script:
    """
    # Create output directory
    mkdir -p ${sample_id}_vibrant

    # Run VIBRANT for phage identification and lifestyle prediction
    VIBRANT_run.py \\
        -i ${assembly} \\
        -t ${task.cpus} \\
        -folder ${sample_id}_vibrant

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vibrant: "4.0.0"
    END_VERSIONS
    """
}
