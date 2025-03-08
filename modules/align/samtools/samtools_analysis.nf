nextflow.enable.dsl=2

include { SAMTOOLS } from '../samtools/main.nf'

workflow SAMTOOLS_ANALYSIS {
    take:
        sam_output
        reference_ch

    main:
        processed_bam = sam_output
            .combine(reference_ch) // Ensure proper pairing with reference
            | SAMTOOLS  // Pass to SAMTOOLS process

    emit:
        bam_output = processed_bam.bam_output
        sorted_bam = processed_bam.sorted_bam
        bam_index = processed_bam.bam_index
        reference_index = processed_bam.reference_index
}