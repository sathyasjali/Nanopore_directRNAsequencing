nextflow.enable.dsl=2

process SAMTOOLS_BCFTOOLS {
    tag "$sample_id"
    publishDir "${params.out_dir}/samtools_bcftools", mode: 'copy'

    input:
        tuple val(sample_id), path(sam_file), path(reference_genome)

    output:
        tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), path("${sample_id}.txt"), path("${reference_genome}")

    script:
    """
    echo "Processing BAM file: ${sam_file}"
    echo "Reference genome: ${reference_genome}"

    # Convert SAM to BAM
    samtools view -Sb '${sam_file}' > '${sample_id}.bam'
    
    # Sort BAM file
    samtools sort '${sample_id}.bam' -o '${sample_id}.sorted.bam'
    
    # Index Sorted BAM file
    samtools index '${sample_id}.sorted.bam'

    #Generate flagstat report
    samtools flagstat "${sam_file}" > "${sample_id}.txt"
    """  
}