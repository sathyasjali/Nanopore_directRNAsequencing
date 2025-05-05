
// PROCESS: Generate Genome Coverage Plots using Python
process PLOTTING {
    tag "Genome Coverage Plot"

    container "python:3.9"

    input:
    path coverage_data

    output:
    path "coverage_plot.png", emit: coverage_plot

    script:
    """
    python3 plot_genome_coverage.py ${coverage_data} coverage_plot.png
    """
}