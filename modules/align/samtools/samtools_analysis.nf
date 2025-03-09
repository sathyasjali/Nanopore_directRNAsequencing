nextflow.enable.dsl=2

include { SAMTOOLS } from '../samtools/main.nf'

workflow SAMTOOLS_ANALYSIS {
    take:
        sam_output
        reference_ch

    main:
        processed_bam = sam_output
            .combine(reference_ch)
            | SAMTOOLS

    sorted_bam_ch = processed_bam.sorted_bam
    bam_index_ch = processed_bam.bam_index

    emit:
        bam_output = processed_bam.bam_output
        sorted_bam = sorted_bam_ch
        bam_index = bam_index_ch
        reference_index = processed_bam.reference_index
}
