nextflow.enable.dsl=2

include { BCFTOOLS } from '../bcftools/main.nf'

workflow BCFTOOLS_ANALYSIS {
    take:
        sorted_bam
        reference_ch

    main:
        variant_calling_results = sorted_bam
            .map { tuple(it[0], it[1]) }  // Ensure the tuples
            .combine(reference_ch)
            | BCFTOOLS

    emit:
        bcf_output = variant_calling_results.bcf_output
        vcf_output = variant_calling_results.vcf_output
        consensus_sequence = variant_calling_results.consensus_sequence
}