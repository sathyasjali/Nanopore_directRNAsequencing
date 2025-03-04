#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// SAMTOOLS_BCFTOOLS process
process SAMTOOLS_BCFTOOLS {
    tag "$sample_id"
    publishDir "${params.out_dir}/samtools_bcftools", mode: 'copy'

    input:
        tuple val(sample_id), path(sam_file), path(reference_genome)

    output:
        tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

    script:
    """
    # Convert SAM to BAM
    samtools view -Sb '${sam_file}' > '${sample_id}.bam'
    
    # Sort BAM file
    samtools sort '${sample_id}.bam' -o '${sample_id}.sorted.bam'
    
    # Index BAM file
    samtools index '${sample_id}.sorted.bam'
    """
}

// SAMTOOLS_INDEX process
process SAMTOOLS_INDEX {
    tag "$sample_id"
    publishDir "${params.out_dir}/samtools_bcftools", mode: 'copy'

    input:
        tuple val(sample_id), path(sam_file), path(reference_genome)

    output:
        path("${sample_id}.txt")
        path("${sample_id}.bcf")
        path("${sample_id}.vcf")
        path("${sample_id}.consensus.fa")

    script:
    """
    # Debugging: Print paths
    echo "Processing sample: ${sample_id}"
    echo "SAM file: ${sam_file}"
    echo "Reference genome: ${reference_genome}"
    
    # Generate flagstat report
    samtools flagstat '${sam_file}' > '${sam_file.baseName}.txt'
    
    # Ensure the reference genome is indexed
    samtools faidx '${reference_genome}'

    # Generate BCF using samtools mpileup and bcftools call
    samtools mpileup -f '${reference_genome}' '${sam_file}' | \
    bcftools call -mv -Ob -o '${sam_file.baseName}.bcf'

    # Convert BCF to VCF
    bcftools view '${sam_file.baseName}.bcf' -Ov -o '${sam_file.baseName}.vcf'

    # Generate consensus sequence
    bcftools consensus -f '${reference_genome}' '${sam_file.baseName}.vcf' > '${sam_file.baseName}.consensus.fa'
    """
}