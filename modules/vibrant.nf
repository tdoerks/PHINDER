process VIBRANT {
    tag "$sample_id"
    publishDir "${params.outdir}/vibrant", mode: 'copy'
    // staphb/vibrant has KEGG + VOG databases bundled — no external db setup needed
    container = 'staphb/vibrant'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_vibrant"), emit: results
    path "${sample_id}_vibrant/VIBRANT_${assembly.baseName}/VIBRANT_results_${assembly.baseName}/VIBRANT_genome_quality_${assembly.baseName}.tsv", emit: quality, optional: true
    path "${sample_id}_vibrant/VIBRANT_${assembly.baseName}/VIBRANT_phages_${assembly.baseName}/${assembly.baseName}.phages_combined.fna", emit: phages, optional: true
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p ${sample_id}_vibrant

    VIBRANT_run.py \\
        -i ${assembly} \\
        -t ${task.cpus} \\
        -folder ${sample_id}_vibrant

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vibrant: \$(VIBRANT_run.py --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """
}
