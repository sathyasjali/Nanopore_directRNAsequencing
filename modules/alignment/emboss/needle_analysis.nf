nextflow.enable.dsl=2
nextflow.enable.dsl=2

include { EMBOSS } from '../emboss/main.nf'

workflow NEEDLE_ANALYSIS {
    take:
        reference_ch
        consensus_sequence_ch

    main:
        needle_results = consensus_sequence_ch
            .map { tuple(it[0], it[1]) }  // Ensure tuples
            .combine(reference_ch)
            | EMBOSS

    emit:
        needle_output = needle_results
        // if you want to extract like tuple use below
        // needle_outpu = needle_results.map { result -> result[0], result[1] }
}
