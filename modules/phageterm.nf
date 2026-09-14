process PHAGETERM {
    tag "$sample_id"
    publishDir "${params.outdir}/phageterm/${sample_id}", mode: 'copy'
    container = 'docker://eluengor/phageterm:1.0.12-tfm'

    input:
    tuple val(sample_id), path(reads), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_phageterm_strategy.txt"), emit: strategy
    path "${sample_id}_phageterm_strategy.txt", emit: report
    path "versions.yml", emit: versions

    script:
    def read1 = reads[0]
    def read2 = reads[1]
    """
    # PhageTerm outputs to current directory using the -n prefix
    SCRIPT=\$(command -v PhageTerm.py 2>/dev/null || command -v PhageTerm 2>/dev/null || echo "PhageTerm.py")

    "\$SCRIPT" \\
        -f ${read1} \\
        -p ${read2} \\
        -r ${assembly} \\
        -n ${sample_id} \\
        -c ${task.cpus} 2>&1 || true

    # Extract packaging strategy from output files
    if [ -f "${sample_id}_packaging_strategy.txt" ]; then
        cp "${sample_id}_packaging_strategy.txt" "${sample_id}_phageterm_strategy.txt"
    elif [ -f "${sample_id}_results.txt" ]; then
        grep -iE "packaging|strategy|termini" "${sample_id}_results.txt" | head -3 > "${sample_id}_phageterm_strategy.txt" || echo "Unknown" > "${sample_id}_phageterm_strategy.txt"
    else
        echo "No packaging signal detected" > "${sample_id}_phageterm_strategy.txt"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phageterm: "1.0.12"
    END_VERSIONS
    """
}
