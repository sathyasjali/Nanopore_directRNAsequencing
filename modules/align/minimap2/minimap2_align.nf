nextflow.enable.dsl=2

include { MINIMAP2 } from '../minimap2/main.nf'

workflow MINIMAP2_ANALYSIS {
    take:
    nanofilt_out
    reference_ch

    main:
    sam_output = nanofilt_out
                .map { sample -> tuple(sample[0], sample[1]) }  // Ensure tuples
                .combine(reference_ch)  // Use combine instead of cross
                | MINIMAP2    
    emit:
        sam_output
}