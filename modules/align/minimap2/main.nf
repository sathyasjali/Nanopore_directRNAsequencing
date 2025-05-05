nextflow.enable.dsl=2

process MINIMAP2 {
    tag "$sample_id"
    publishDir "${params.out_dir}/minimap2", mode: 'copy'

    input:
        tuple val(sample_id), path(filtered_fastq), path(reference)

    output:
        tuple val(sample_id), path("${sample_id}.sam"), emit: sam_output

    script:
    """
    minimap2 -ax map-ont --secondary=no -k5 ${reference} ${filtered_fastq} > ${sample_id}.sam
    """
}