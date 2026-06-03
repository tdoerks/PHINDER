/*
========================================================================================
    PHINDER MAIN WORKFLOW
========================================================================================
*/

include { DOWNLOAD_SRA } from '../modules/sra_download'
include { FASTQC } from '../modules/fastqc'
include { FASTP } from '../modules/fastp'
include { UNICYCLER } from '../modules/unicycler'
include { QUAST } from '../modules/quast'
include { CHECKV } from '../modules/checkv'
include { PHAROKKA } from '../modules/pharokka'
include { VIBRANT } from '../modules/vibrant'
include { DIAMOND_PROPHAGE } from '../modules/diamond_prophage'
include { PHANOTATE } from '../modules/phanotate'
include { MULTIQC } from '../modules/multiqc'
include { PHINDER_SUMMARY } from '../modules/phinder_summary'

workflow PHINDER_PIPELINE {

    // Parse input based on mode
    if (params.input_mode == 'sra') {
        // SRA accession list mode
        ch_srr_list = Channel
            .fromPath(params.input, checkIfExists: true)
            .splitText()
            .map { it.trim() }

        DOWNLOAD_SRA(ch_srr_list)
        ch_input = DOWNLOAD_SRA.out.reads
        ch_versions = DOWNLOAD_SRA.out.versions.first()
    } else {
        // Regular samplesheet mode
        ch_input = parse_samplesheet(params.input, params.input_mode)
        ch_versions = Channel.empty()
    }

    // Initialize channels
    ch_multiqc_files = Channel.empty()

    // STEP 1: Quality Control (if starting from reads or SRA)
    if ((params.input_mode == 'reads' || params.input_mode == 'sra') && !params.skip_fastqc) {
        FASTQC(ch_input)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip)
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    // STEP 2: Read Trimming (if starting from reads or SRA)
    if ((params.input_mode == 'reads' || params.input_mode == 'sra') && !params.skip_fastp) {
        // Transform input for fastp (expects sample_id, read1, read2)
        ch_fastp_input = ch_input.map { sample_id, reads ->
            [sample_id, reads[0], reads[1]]
        }
        FASTP(ch_fastp_input)
        ch_trimmed = FASTP.out.reads
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json)
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
    } else {
        ch_trimmed = ch_input
    }

    // STEP 3: Assembly (if starting from reads or SRA)
    if ((params.input_mode == 'reads' || params.input_mode == 'sra') && !params.skip_assembly) {
        if (params.assembler == 'unicycler') {
            // Transform input for unicycler (expects sample_id, read1, read2)
            ch_unicycler_input = ch_trimmed.map { sample_id, reads ->
                [sample_id, reads[0], reads[1]]
            }
            UNICYCLER(ch_unicycler_input)
            ch_assemblies = UNICYCLER.out.assembly
            ch_versions = ch_versions.mix(UNICYCLER.out.versions.first())
        }
        // Add SPAdes option here if needed in future
    } else if (params.input_mode == 'assembly') {
        ch_assemblies = ch_input
    }

    // STEP 4: Assembly Quality Assessment
    if (!params.skip_assembly) {
        QUAST(ch_assemblies)
        ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.tsv)
        ch_versions = ch_versions.mix(QUAST.out.versions.first())
    }

    // STEP 5: CheckV Quality Assessment
    if (!params.skip_checkv) {
        CHECKV(ch_assemblies)
        ch_versions = ch_versions.mix(CHECKV.out.versions.first())
    }

    // STEP 6: Pharokka Annotation
    if (!params.skip_pharokka) {
        PHAROKKA(ch_assemblies)
        ch_versions = ch_versions.mix(PHAROKKA.out.versions.first())
    }

    // STEP 7: VIBRANT Lifestyle Prediction
    if (!params.skip_vibrant) {
        VIBRANT(ch_assemblies)
        ch_versions = ch_versions.mix(VIBRANT.out.versions.first())
    }

    // STEP 8: DIAMOND Prophage Comparison
    if (!params.skip_diamond) {
        DIAMOND_PROPHAGE(ch_assemblies)
        ch_versions = ch_versions.mix(DIAMOND_PROPHAGE.out.versions.first())
    }

    // STEP 9: PHANOTATE Gene Prediction
    if (!params.skip_phanotate) {
        PHANOTATE(ch_assemblies)
        ch_versions = ch_versions.mix(PHANOTATE.out.versions.first())
    }

    // STEP 10: MultiQC Report
    MULTIQC(ch_multiqc_files.collect().ifEmpty([]))
    ch_versions = ch_versions.mix(MULTIQC.out.versions)

    // STEP 11: PHINDER Summary Report
    // Wait for all analyses to complete, then generate summary
    ch_all_complete = Channel.empty()
    if (!params.skip_checkv) {
        ch_all_complete = ch_all_complete.mix(CHECKV.out.quality)
    }
    if (!params.skip_pharokka) {
        ch_all_complete = ch_all_complete.mix(PHAROKKA.out.functions)
    }
    if (!params.skip_vibrant) {
        ch_all_complete = ch_all_complete.mix(VIBRANT.out.quality)
    }
    if (!params.skip_assembly) {
        ch_all_complete = ch_all_complete.mix(QUAST.out.tsv)
    }

    // Generate summary when all samples are done
    PHINDER_SUMMARY(ch_all_complete.collect())
    ch_versions = ch_versions.mix(PHINDER_SUMMARY.out.versions)

    emit:
    versions = ch_versions
    multiqc_report = MULTIQC.out.report
    phinder_summary = PHINDER_SUMMARY.out.html
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

def parse_samplesheet(input_file, input_mode) {
    if (input_mode == 'reads') {
        return Channel
            .fromPath(input_file, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def sample_id = row.sample
                def read1 = file(row.read1, checkIfExists: true)
                def read2 = file(row.read2, checkIfExists: true)
                [sample_id, [read1, read2]]
            }
    } else if (input_mode == 'assembly') {
        return Channel
            .fromPath(input_file, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def sample_id = row.sample
                def assembly = file(row.assembly, checkIfExists: true)
                [sample_id, assembly]
            }
    } else {
        error "Invalid input_mode: ${input_mode}. Must be 'reads', 'assembly', or 'sra'"
    }
}
