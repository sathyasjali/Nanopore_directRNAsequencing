nextflow.enable.dsl=2

include { ERROR_PLOT } from './main.nf'

workflow ERROR_PLOT_ANALYSIS {
    take:
        error_vcf_ch
        coverage_ch

    main:
        error_plot_results = error_vcf_ch.combine(coverage_ch)
            .map { sample_id, vcf_file, coverage_file -> 
                tuple(sample_id, vcf_file, coverage_file) 
            } | ERROR_PLOT

    emit:
        error_plot_results
}