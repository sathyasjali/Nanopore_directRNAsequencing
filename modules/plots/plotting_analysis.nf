nextflow.enable.dsl=2

include { PLOTTING } from '../plots/main.nf'
workflow PLOTTING_ANALYSIS {
    take:
        global_dist
    main:
        plot_results = PLOTTING(
            global_dist
        )

    emit:
        plot_output = plot_results.plot_output
}