process SAMTOOLS_INDEX {
    tag "$sample_id"
    publishDir "${params.out_dir}/samtools_bcftools", mode: 'copy'

    input:
        tuple val(sample_id), path(sorted_bam), path(reference_genome)

    output:
        path("${sample_id}.mpileup")
        

    script:
    """
    echo "DEBUG: Processing sample: ${sample_id}"
    echo "DEBUG: Sorted BAM file: ${sorted_bam}"
    echo "DEBUG: Reference genome: ${reference_genome}"

    # Ensure required files exist
    if [[ ! -f "${sorted_bam}" ]]; then
        echo "ERROR: Sorted BAM file does not exist!" >&2
        exit 1
    fi

    if [[ ! -f "${reference_genome}" ]]; then
        echo "ERROR: Reference genome does not exist!" >&2
        exit 1
    fi

    # Index the BAM file if index is missing
    if [[ ! -f "${sorted_bam}.bai" ]]; then
        echo "Indexing BAM file..."
        samtools index "${sorted_bam}"
    fi

    # Ensure reference genome is indexed
    if [[ ! -f "${reference_genome}.fai" ]]; then
        echo "Indexing reference genome..."
        samtools faidx "${reference_genome}"
    fi

    # Generate BCF using samtools mpileup and bcftools call
    samtools mpileup -f "${reference_genome}" "${sorted_bam}" > "${sample_id}.mpileup"
    
    """
}