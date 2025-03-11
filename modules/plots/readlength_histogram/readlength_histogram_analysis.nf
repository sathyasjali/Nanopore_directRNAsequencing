nextflow.enable.dsl=2

// Include the module
include { READLENGTH_HISTOGRAM } from './main.nf'

workflow READLENGTH_HISTOGRAM_ANALYSIS {
    
    take:
    filtered_fastq_ch

    main:
    histogram_results = filtered_fastq_ch | READLENGTH_HISTOGRAM

    emit:
    histogram_results
}
