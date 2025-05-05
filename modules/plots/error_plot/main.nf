process ERROR_PLOT {
    tag "$sample_id"
    publishDir "${params.out_dir}/error_plots", mode: 'copy'

    input:
        tuple val(sample_id), path(vcf_file), path(coverage_file)

    output:
        path("${sample_id}_error_plot.png")

    script:
    def output_file = "${sample_id}_error_plot.png"

    """
    # Ensure required Python packages are installed
    # Run the Python script to generate the error plot
    python3 ${params.error_script} \\
        --vcf ${vcf_file} \\
        --coverage ${coverage_file} \\
        --output "\$PWD/${output_file}" \\
        --sample_id ${sample_id}

    # Verify the output file was created
    [ -f "${output_file}" ] || { echo "Error: Output file ${output_file} was not created!" >&2; exit 1; }
    """
}
