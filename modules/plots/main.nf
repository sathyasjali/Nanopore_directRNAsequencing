
// PROCESS: Generate Genome Coverage Plots using Python

process PLOTTING {
    container 'python:3.9-slim'
    tag "$sample_id"
    publishDir "${params.out_dir}/plotting", mode: 'copy'


    input:
        path(global_dist)

    output:
        path("coverage_plot.html"), emit: plot_output

    script:
    """
    # Copy the script into the working directory
    cp ${params.parseScript} ./plot-dist.py

    python3 plot-dist.py \
        ${global_dist} \
        --output coverage_plot.html
    """
}