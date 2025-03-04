process REFERENCE_INDEX {
    tag "index_ref"
    publishDir "${params.out_dir}/reference", mode: 'copy'

    input:
        path(reference_genome)

    output:
        path("${reference_genome}.fai")  // Indexed reference genome

    script:
    """
    echo "Indexing reference genome: ${reference_genome}"
    samtools faidx '${reference_genome}'
    """
}