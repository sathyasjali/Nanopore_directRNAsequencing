nextflow.enable.dsl=2

include { MOSDEPTH } from './main.nf'
workflow MOSDEPTH_ANALYSIS {
    take:
        sorted_bam_ch
        bam_index_ch

    main:
        mosdepth_results = sorted_bam_ch
            .combine(bam_index_ch, by: 0)  // Ensures pairing by sample_id
            .map { sample_id, sorted_bam, bam_index -> 
                tuple(sample_id, sorted_bam, bam_index)
            }
            | MOSDEPTH  // Correct tuple structure

    emit:
        global_dist = mosdepth_results.global_dist
        region_dist = mosdepth_results.region_dist
        summary = mosdepth_results.summary
        per_base_coverage = mosdepth_results.per_base_coverage
        per_region_coverage = mosdepth_results.per_region_coverage
        threshold_coverage = mosdepth_results.threshold_coverage
}
