nextflow.enable.dsl=2
include { FASTQC } from '../fastqc/main.nf'
workflow FASTQC_ANALYSIS {
    take:
    fastq_ch

    main:
    fastq_ch.map { tuple(it.getName().replaceAll(/\.fastq$/, ''), it) } | FASTQC
    
    emit:
    fastqc_reports_zip = FASTQC.out.fastqc_zip
    fastqc_reports_html = FASTQC.out.fastqc_html
} 