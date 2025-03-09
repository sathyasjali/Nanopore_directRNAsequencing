
// PROCESS: Generate Genome Coverage Plots using Python

process PLOTTING {
    tag "$sample_id"
    container 'python:3.9-slim'

    input:
        val sample_id
        path(global_dist)
        path(region_dist)
        path(summary)
        path(per_base_coverage)
        path(per_region_coverage)
        path(threshold_coverage)

    output:
        path("${sample_id}_plot.png"), emit: plot_output

    script:
    """
    chmod +x plot_script.py
    python plot_script.py \
        --global_dist ${global_dist} \
        --region_dist ${region_dist} \
        --summary ${summary} \
        --per_base ${per_base_coverage} \
        --per_region ${per_region_coverage} \
        --threshold ${threshold_coverage} \
        --output ${sample_id}_plot.png
    """
}
