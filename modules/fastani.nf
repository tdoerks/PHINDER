process FASTANI {
    publishDir "${params.outdir}/fastani", mode: 'copy'
    container = 'staphb/fastani:1.34'

    input:
    val(sample_ids)
    path(assemblies)

    output:
    path "fastani_all_vs_all.tsv", emit: matrix
    path "versions.yml", emit: versions

    script:
    def ids = sample_ids instanceof List ? sample_ids : [sample_ids]
    def fns = (assemblies instanceof List ? assemblies : [assemblies]).collect { it.name }
    // Single-line cp commands joined with && to avoid multi-line interpolation
    // (multi-line interpolation breaks Nextflow's stripIndent and heredoc terminators)
    def renames = [ids, fns].transpose().collect { id, fn -> "cp '${fn}' '${id}.fasta'" }.join(' && ')
    """
    ${renames}

    ls *.fasta > genome_list.txt

    fastANI \\
        --ql genome_list.txt \\
        --rl genome_list.txt \\
        -o fastani_all_vs_all.tsv \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastani: \$(fastANI --version 2>&1 | head -1)
    END_VERSIONS
    """
}
