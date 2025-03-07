nextflow.enable.dsl=2
include { NANOFILT } from '../nanofilt/main.nf'
workflow NANOFILT_ANALYSIS {
    take:
    fastq_ch

    main:
    filtered_results = fastq_ch.map { tuple(it.baseName, it) } | NANOFILT


    emit:
    filtered_results

}