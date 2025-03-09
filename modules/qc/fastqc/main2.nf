process FASTQC2 {
    tag "$sample_id"
    publishDir "${params.out_dir}/fastqc2", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_file)

    output:
    tuple val(sample_id), path("${sample_id}_fastqc.zip"), emit: fastqc_filtered_zip
    tuple val(sample_id), path("${sample_id}_fastqc.html"), emit: fastqc_filtered_html

    script:
    """
    fastqc $fastq_file --outdir .
    BASENAME=\$(basename ${fastq_file} .fastq)
    mv \${BASENAME}_fastqc.zip ${sample_id}_fastqc.zip
    mv \${BASENAME}_fastqc.html ${sample_id}_fastqc.html
    """
}