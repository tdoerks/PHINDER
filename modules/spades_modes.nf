/*
========================================================================================
    SPADES_MODE — SPAdes assembly with explicit mode selection
    Used by the spades-mode-test comparison workflow.
    Tries all 6 modes per sample; graceful failure writes a stub FASTA.
========================================================================================
*/

process SPADES_MODE {
    tag "${sample_id}:${mode}"
    publishDir "${params.outdir}/spades_modes/${sample_id}__${mode}", mode: 'copy'
    container = 'quay.io/biocontainers/spades:3.15.5--h95f258a_1'
    errorStrategy = 'ignore'

    input:
    tuple val(sample_id), path(read1), path(read2), val(mode)

    output:
    tuple val(sample_id), val(mode), path("${sample_id}__${mode}_assembly.fasta"), emit: assembly
    path "${sample_id}__${mode}_spades.log",                                        emit: log
    path "versions.yml",                                                            emit: versions

    script:
    // --careful only works with standard mode; the others use --only-assembler
    def flags = mode == 'isolate'   ? '--isolate --only-assembler'   :
                mode == 'meta'      ? '--meta --only-assembler'      :
                mode == 'metaviral' ? '--metaviral --only-assembler' :
                mode == 'rnaviral'  ? '--rnaviral --only-assembler'  :
                mode == 'careful'   ? '--careful'                    :
                                      '--only-assembler'  // 'standard'
    """
    set +e
    spades.py \\
        -1 ${read1} \\
        -2 ${read2} \\
        -o spades_out \\
        --threads ${task.cpus} \\
        --memory ${task.memory.toGiga()} \\
        ${flags} \\
        2>&1 | tee ${sample_id}__${mode}_spades.log
    SPADES_EXIT=\$?

    # Pick the best output file for this mode:
    #   metaviral → metaviral.contigs.fasta (viral-specific output)
    #   rnaviral  → rna_viral.contigs.fasta
    #   others    → scaffolds.fasta → contigs.fasta (fallback)
    ASSEMBLY=""
    if [ "${mode}" = "metaviral" ] && [ -f spades_out/metaviral.contigs.fasta ] && [ -s spades_out/metaviral.contigs.fasta ]; then
        ASSEMBLY=spades_out/metaviral.contigs.fasta
    elif [ "${mode}" = "rnaviral" ] && [ -f spades_out/rna_viral.contigs.fasta ] && [ -s spades_out/rna_viral.contigs.fasta ]; then
        ASSEMBLY=spades_out/rna_viral.contigs.fasta
    elif [ -f spades_out/scaffolds.fasta ] && [ -s spades_out/scaffolds.fasta ]; then
        ASSEMBLY=spades_out/scaffolds.fasta
    elif [ -f spades_out/contigs.fasta ] && [ -s spades_out/contigs.fasta ]; then
        ASSEMBLY=spades_out/contigs.fasta
    fi

    if [ -n "\$ASSEMBLY" ]; then
        cp "\$ASSEMBLY" ${sample_id}__${mode}_assembly.fasta
        echo "MODE_STATUS: success (exit \$SPADES_EXIT)" >> ${sample_id}__${mode}_spades.log
    else
        # Write a stub so downstream CheckV/QUAST receive a file (they gracefully handle it)
        printf ">ASSEMBLY_FAILED mode=${mode} exit=\$SPADES_EXIT\\nN\\n" > ${sample_id}__${mode}_assembly.fasta
        echo "MODE_STATUS: failed (exit \$SPADES_EXIT, no usable output)" >> ${sample_id}__${mode}_spades.log
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: "3.15.5"
        mode: "${mode}"
    END_VERSIONS
    """
}
