process FASTQC2 {
    tag "$sample_id"
    publishDir "${params.out_dir}/fastqc2", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_file)

    output:
    tuple val(sample_id), path("${sample_id}_fastqc_filtered.zip"), emit: fastqc_filtered_zip
    tuple val(sample_id), path("${sample_id}_fastqc_filtered.html"), emit: fastqc_filtered_html

    script:
    """
    fastqc $fastq_file --outdir .
    mv ${fastq_file.baseName}_fastqc.zip ${sample_id}_fastqc_filtered.zip
    mv ${fastq_file.baseName}_fastqc.html ${sample_id}_fastqc_filtered.html
    """
}