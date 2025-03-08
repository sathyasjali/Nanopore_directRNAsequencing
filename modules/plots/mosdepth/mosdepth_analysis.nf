nextflow.enable.dsl=2

include { MOSDEPTH } from './main.nf'

workflow MOSDEPTH_ANALYSIS {
    take:
        sorted_bam_ch
        bam_index_ch

    main:
        paired_inputs = sorted_bam_ch
            .combine(bam_index_ch) // Ensure correct pairing of BAM and index

        coverage_results = paired_inputs | MOSDEPTH

    emit:
        coverage_summary = coverage_results.map { it[1] } // Summary output
        coverage_data = coverage_results.map { it[2] } // Coverage distribution data
}
