nextflow.enable.dsl=2

process MULTIQC2 {
    tag "MultiQC Report2"

    publishDir "${params.out_dir}/multiqc2", mode: 'copy'
    
    input:
    path fastqc_filtered_reports_zip

    output:
    path "multiqc_filtered_report.html", emit: multiqc_filtered_html
    path "multiqc_filtered_data", emit: multiqc_filtered_data

    script:
    """
    multiqc ${fastqc_filtered_reports_zip} -o ./
    """
}