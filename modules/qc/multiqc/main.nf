nextflow.enable.dsl=2

process MULTIQC {
    tag "MultiQC Report"
    publishDir "${params.out_dir}/multiqc", mode: 'copy'
    
    input:
    path fastqc_reports  // Accepts multiple files as a list
    
    output:
    path "multiqc_report.html", emit: multiqc_html
    
    script:
    """
    fastqc_reports.view { println "Files going into MultiQC: $it" }

    // multiqc ${fastqc_reports} -o ./
    """
}