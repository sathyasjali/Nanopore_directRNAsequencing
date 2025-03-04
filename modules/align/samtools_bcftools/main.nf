#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// SAMTOOLS_BCFTOOLS process
process SAMTOOLS_BCFTOOLS {
    tag "$sample_id"
    publishDir "${params.out_dir}/samtools_bcftools", mode: 'copy'

    input:
        tuple val(sample_id), path(sam_file)  // removed reference_genome

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
// process SAMTOOLS_INDEX {
//     tag "$sample_id"
//     publishDir "${params.out_dir}/samtools_index", mode: 'copy'

//     input:
//         // Updated to take the sorted BAM file from SAMTOOLS_BCFTOOLS output
//         tuple val(sample_id), path(sorted_bam), path(reference_genome)

//     output:
//         path("${sample_id}.txt")
//         path("${sample_id}.bcf")
//         path("${sample_id}.vcf")
//         path("${sample_id}.consensus.fa")

//     script:
//     """
//     # Generate flagstat report using the sorted BAM
//     samtools flagstat '${sorted_bam}' > '${sorted_bam.baseName}.txt'
    
//     # Ensure the reference genome is indexed
//     samtools faidx '${reference_genome}'

//     # Generate BCF using samtools mpileup and bcftools call
//     samtools mpileup -f '${reference_genome}' '${sorted_bam}' | \
//     bcftools call -mv -Ob -o '${sorted_bam.baseName}.bcf'

//     # Convert BCF to VCF
//     bcftools view '${sorted_bam.baseName}.bcf' -Ov -o '${sorted_bam.baseName}.vcf'

//     # Generate consensus sequence
//     bcftools consensus -f '${reference_genome}' '${sorted_bam.baseName}.vcf' > '${sorted_bam.baseName}.consensus.fa'
//     """
// }

process SAMTOOLS_INDEX {
    tag "$sample_id"
    publishDir "${params.out_dir}/samtools_index", mode: 'copy'

    input:
        tuple val(sample_id), path(sorted_bam), path(reference_genome)

    output:
        path("${sample_id}.flagstat.txt")
        path("${sample_id}.bcf")
        path("${sample_id}.vcf")
        path("${sample_id}.consensus.fa")

    script:
    """
    echo "Processing sample: ${sample_id}"

    # Generate flagstat report
    samtools flagstat '${sorted_bam}' > '${sample_id}.flagstat.txt'
    
    # Check if reference genome index exists, otherwise create it
    if [ ! -f '${reference_genome}.fai' ]; then
        samtools faidx '${reference_genome}'
    fi

    # Generate BCF file
    samtools mpileup -f '${reference_genome}' '${sorted_bam}' | \
    bcftools call -mv -Ob -o '${sample_id}.bcf'

    # Convert BCF to VCF
    bcftools view '${sample_id}.bcf' -Ov -o '${sample_id}.vcf'

    # Generate consensus sequence
    bcftools consensus -f '${reference_genome}' '${sample_id}.vcf' > '${sample_id}.consensus.fa'
    """
}
