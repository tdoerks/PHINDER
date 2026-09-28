process BACPHLIP {
    tag "$sample_id"
    publishDir "${params.outdir}/bacphlip", mode: 'copy'
    container = 'quay.io/biocontainers/mulled-v2-e16bfb0f667f2f3c236b32087aaf8c76a0cd2864:c64689d7d5c51670ff5841ec4af982edbe7aa406-0'
    errorStrategy = 'ignore'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}.bacphlip"), emit: results, optional: true
    path "${sample_id}.bacphlip", emit: predictions, optional: true
    path "versions.yml", emit: versions

    script:
    """
    # Filter contigs < 500bp — short contigs produce empty AA translations that crash hmmsearch
    awk '/^>/{if(seq && length(seq)>=500) print h"\\n"seq; h=\$0; seq=""; next} {seq=seq\$0} END{if(seq && length(seq)>=500) print h"\\n"seq}' ${assembly} > filtered_assembly.fasta
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
