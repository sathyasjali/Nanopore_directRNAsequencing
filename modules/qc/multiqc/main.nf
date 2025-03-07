nextflow.enable.dsl=2

process MULTIQC {
    tag "MultiQC Report"

    publishDir "${params.out_dir}/multiqc", mode: 'copy'
    
    input:
    path fastqc_zip

    output:
    path "multiqc_report.html", emit: multiqc_html
    path "multiqc_data", emit: multiqc_data

    script:
    """
    multiqc ${fastqc_zip} -o ./
    """
}