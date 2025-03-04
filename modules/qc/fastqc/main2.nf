#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FASTQC2 {
    tag "$sample_id"
    publishDir "${params.out_dir}/fastqc2", mode: 'copy'

    input:
        tuple val(sample_id), path("${sample_id}_filtered.fastq")

    output:
        tuple val(sample_id), path("${sample_id}_filtered_fastqc.zip"), path("${sample_id}_filtered_fastqc.html")

    script:
    """
    fastqc ${sample_id}_filtered.fastq --outdir .
    """
}