#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FASTQC {
    tag "$sample_id"
    publishDir "${params.out_dir}/fastqc", mode: 'copy'

    input:
        tuple val(sample_id), path(fastq_file)

    output:
        tuple val(sample_id), path("${sample_id}_fastqc.zip"), path("${sample_id}_fastqc.html")

    script:
    """
    fastqc ${fastq_file} --outdir .
    """
}