nextflow.enable.dsl=2

include { MULTIQC2 } from '../multiqc/main2.nf'

workflow MULTIQC_ANALYSIS2 {
take:
fastqc_filtered_reports_zip

main:
fastqc_filtered_reports_zip.collect { it[1] } | MULTIQC2 
    
emit:
multiqc_filtered_report_html = MULTIQC2.out.multiqc_filtered_html
multiqc_filtered_data_folder = MULTIQC2.out.multiqc_filtered_data
}