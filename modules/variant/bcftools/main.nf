process BCFTOOLS {
    tag "$sample_id"
    publishDir "${params.out_dir}/bcftools", mode: 'copy'

    input:
        tuple val(sample_id), path(sorted_bam), path(reference)

    output:
        tuple val(sample_id), path("${sample_id}.bcf"), emit: bcf_output
        tuple val(sample_id), path("${sample_id}.vcf.gz"), emit: vcf_output
        tuple val(sample_id), path("${sample_id}_consensus.fa"), emit: consensus_sequence

    script:
    """
    # Generate BCF file
    bcftools mpileup -Ou -f ${reference} ${sorted_bam} | bcftools call -mv -Ob -o ${sample_id}.bcf

    # Convert BCF to VCF and compress it
    bcftools view ${sample_id}.bcf > ${sample_id}.vcf
    bgzip -c ${sample_id}.vcf > ${sample_id}.vcf.gz

    # Index the VCF file
    bcftools index ${sample_id}.vcf.gz

    # Generate consensus sequence using the indexed VCF file
    bcftools consensus -f ${reference} ${sample_id}.vcf.gz > ${sample_id}_consensus.fa
    """
}
