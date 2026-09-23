process BACPHLIP {
    tag "$sample_id"
    publishDir "${params.outdir}/bacphlip", mode: 'copy'
    container = 'quay.io/biocontainers/mulled-v2-e16bfb0f667f2f3c236b32087aaf8c76a0cd2864:c64689d7d5c51670ff5841ec4af982edbe7aa406-0'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}.bacphlip"), emit: results
    path "${sample_id}.bacphlip", emit: predictions
    path "versions.yml", emit: versions

    script:
    """
    # Filter contigs < 200bp to avoid empty 6-frame translations crashing BacPhlip
    awk '/^>/{header=\$0; seq=""} !/^>/{seq=seq\$0} /^>/ && seq!="" && length(seq)>=200{print header"\\n"seq} END{if(seq!="" && length(seq)>=200) print header"\\n"seq}' ${assembly} > filtered_assembly.fasta
    # Fall back to full assembly if filtering removed everything
    [ -s filtered_assembly.fasta ] || cp ${assembly} filtered_assembly.fasta

    contig_count=\$(grep -c '^>' filtered_assembly.fasta)
    multi_flag=\$([ "\$contig_count" -gt 1 ] && echo "--multi_fasta" || echo "")

    bacphlip -i filtered_assembly.fasta \$multi_flag

    # bacphlip names output after the input file — rename to sample_id
    mv filtered_assembly.fasta.bacphlip ${sample_id}.bacphlip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bacphlip: \$(bacphlip --version 2>&1 | sed 's/bacphlip //g')
    END_VERSIONS
    """
}
