#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process NANOFILT {
    tag "$sample_id"
    publishDir "${params.out_dir}/nanofilt", mode: 'copy'

    input:
        tuple val(sample_id), path(fastq_file)

    output:
        tuple val(sample_id), path("${sample_id}_filtered.fastq")

    script:
    """
    cat ${fastq_file} | NanoFilt -q 5 -l 50 > ${sample_id}_filtered.fastq
    """
}