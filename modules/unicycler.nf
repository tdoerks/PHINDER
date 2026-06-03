process UNICYCLER {
    tag "$sample_id"
    publishDir "${params.outdir}/assemblies", mode: 'copy'
    container = 'docker://staphb/unicycler:0.5.0'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_assembly.fasta"), emit: assembly
    path "${sample_id}_unicycler.log", emit: log
    path "versions.yml", emit: versions

    script:
    """
    unicycler \\
        -1 ${read1} \\
        -2 ${read2} \\
        -o unicycler_output \\
        -t ${task.cpus} \\
        --verbosity 2

    # Rename output assembly
    mv unicycler_output/assembly.fasta ${sample_id}_assembly.fasta
    mv unicycler_output/unicycler.log ${sample_id}_unicycler.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        unicycler: \$(unicycler --version | sed 's/Unicycler v//g')
    END_VERSIONS
    """
}
