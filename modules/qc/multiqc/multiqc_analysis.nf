nextflow.enable.dsl=2

include { MULTIQC } from '../multiqc/main.nf'

workflow MULTIQC_ANALYSIS {
take:
fastqc_reports_zip

main:
fastqc_reports_zip.collect { it[1] } | MULTIQC 
    
emit:
multiqc_report_html = MULTIQC.out.multiqc_html
multiqc_data_folder = MULTIQC.out.multiqc_data
}