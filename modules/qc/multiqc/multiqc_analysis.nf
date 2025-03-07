nextflow.enable.dsl=2

include { MULTIQC } from '../multiqc/main.nf'

workflow MULTIQC_ANALYSIS {
    take:
    fastqc_reports  // Now takes a merged channel of all FASTQC outputs

    main:
    MULTIQC(fastqc_reports)

    emit:
    multiqc_html = MULTIQC.out.multiqc_html
}