process FASTQC {
    tag "$sample_id"
    publishDir "${params.out_dir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_file)

    output:
    tuple val(sample_id), path("${sample_id}_fastqc.zip"), emit: fastqc_zip
    tuple val(sample_id), path("${sample_id}_fastqc.html"), emit: fastqc_html

    script:
    """
    fastqc $fastq_file --outdir .
    mv ${fastq_file.baseName}_fastqc.zip ${sample_id}_fastqc.zip
    mv ${fastq_file.baseName}_fastqc.html ${sample_id}_fastqc.html
    """
}