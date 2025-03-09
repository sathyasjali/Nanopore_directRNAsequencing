nextflow.enable.dsl=2

include { PLOTTING } from '../plots/main.nf'
workflow PLOTTING_ANALYSIS {
    take:
        sample_id          // ✅ Add missing sample_id
        global_dist
        region_dist
        summary
        per_base_coverage
        per_region_coverage
        threshold_coverage

    main:
        plot_results = PLOTTING(
            sample_id,      // ✅ Ensure sample_id is passed
            global_dist,
            region_dist,
            summary,
            per_base_coverage,
            per_region_coverage,
            threshold_coverage
        )

    emit:
        plot_output = plot_results.plot_output
}