nextflow.enable.dsl=2
process EMBOSS {
    tag "$sample_id"
    publishDir "${params.out_dir}/needle", mode: 'copy'

    input:
        tuple val(sample_id), path(reference), path(consensus_sequence)

    output:
        tuple val(sample_id), path("${sample_id}_needle.txt"), emit: needle_output

    script:
    """
    needle -asequence ${reference} -bsequence ${consensus_sequence} \
           -gapopen 10.0 -gapextend 0.5 \
           -outfile ${sample_id}_needle.txt
    """
}
