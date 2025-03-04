nextflow.enable.dsl=2

process MULTIQC {
    tag "multiqc"
    publishDir "${params.out_dir}/multiqc", mode: 'copy'

    input:
        path fastqc_results
        path nanofilt_results

    output:
        path 'multiqc_report.html'

    script:
    """
    multiqc ${fastqc_results} ${nanofilt_results} -o ./
    """
}