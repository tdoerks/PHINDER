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
    // Build a sample_id → filename map so output uses readable names
    def ids  = sample_ids instanceof List ? sample_ids : [sample_ids]
    def fns  = assemblies  instanceof List ? assemblies.collect { it.name } : [assemblies.name]
    def pairs = [ids, fns].transpose().collect { id, fn -> "${id}\t${fn}" }.join('\n')
    """
    cat > sample_map.tsv << 'MAPEOF'
    ${pairs}
    MAPEOF

    python3 - << 'PYEOF'
    import shutil, pathlib
    with open('sample_map.tsv') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            sid, fn = line.split('\\t', 1)
            src = pathlib.Path(fn.strip())
            if src.exists():
                shutil.copy(src, f'{sid.strip()}.fasta')
    PYEOF

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
