#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// MINIMAP2 process
process MINIMAP2 {
    tag "$sample_id"
    publishDir "${params.out_dir}/minimap2", mode: 'copy'

    input:
        tuple val(sample_id), path(filtered_fastq), path(reference_genome)

    output:
        tuple val(sample_id), path("${sample_id}.sam")

    script:
    """
    minimap2 -ax map-ont --secondary=no -k5 ${reference_genome} ${filtered_fastq} > ${sample_id}.sam
    """
}