nextflow.enable.dsl=2

include { MINIMAP2 } from '../minimap2/main.nf'

workflow MINIMAP2_ANALYSIS {
    take:
    filtered_results_ch
    reference_ch

main:
    sam_output = filtered_results_ch.combine(reference_ch).map { tuple(it[0], it[1], it[2]) } | MINIMAP2
    
emit:
    sam_output
}