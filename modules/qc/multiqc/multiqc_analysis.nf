nextflow.enable.dsl=2

include { MULTIQC } from '../multiqc/main.nf'

workflow MULTIQC_ANALYSIS {
    take:
    fastqc_html
    fastqc_zip

    main:
    MULTIQC(fastqc_html, fastqc_zip)

    emit:
    multiqc_html = MULTIQC.out.multiqc_html
}