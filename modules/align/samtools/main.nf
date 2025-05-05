#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// SAMTOOLS_BCFTOOLS process
process SAMTOOLS {
    tag "$sample_id"
    publishDir "${params.out_dir}/samtools", mode: 'copy'

    input:
        tuple val(sample_id), path(sam_file), path(reference)
    
    output:
        tuple val(sample_id), path("${sample_id}.bam"), emit: bam_output
        tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: sorted_bam
        tuple val(sample_id), path("${sample_id}_sorted.bam.bai"), emit: bam_index
        tuple val(sample_id), path("${sample_id}_flagstat.txt"), emit: flagstat_output
        tuple val(sample_id), path("${reference}.fai"), emit: reference_index

    script:
    """
    samtools view -bS ${sam_file} > ${sample_id}.bam
    samtools sort -o ${sample_id}_sorted.bam ${sample_id}.bam
    samtools index ${sample_id}_sorted.bam
    samtools flagstat ${sample_id}_sorted.bam > ${sample_id}_flagstat.txt
    samtools faidx ${reference}
    """
}