#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SAMTOOLS_BCFTOOLS {
    tag "$sample_id"
    publishDir "${params.out_dir}/samtools_bcftools", mode: 'copy'

    input:
        tuple val(sample_id), path(sam_file), 
        path(reference_genome)

    output:
        tuple val(sample_id), path("${sample_id}.bam"), path("${reference_genome}.fa.fai")

    script:
    """"
    # Check if input paths are not null
    if [ ! -f ${sam_file} ]; then
        echo "Error: SAM file path is null or file does not exist"
        exit 1
    fi

    if [ ! -f ${reference_genome} ]; then
        echo "Error: Reference genome path is null or file does not exist"
        exit 1
    fi

    # Convert SAM to BAM
    samtools view -bS ${sam_file} > ${sample_id}.bam

    # Index the reference genome
    samtools faidx ${reference_genome}
    """
}

