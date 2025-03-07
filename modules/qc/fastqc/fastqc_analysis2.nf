nextflow.enable.dsl=2
include { FASTQC2 } from '../fastqc/main2.nf'
workflow FASTQC_ANALYSIS2 {
    take:
    nanofilt_out

    main:
    nanofilt_out.map { tuple(it[0], it[1]) } | FASTQC2

    emit:
    fastqc_filtered_reports_zip = FASTQC2.out.fastqc_filtered_zip
    fastqc_filtered_reports_html = FASTQC2.out.fastqc_filtered_html
}