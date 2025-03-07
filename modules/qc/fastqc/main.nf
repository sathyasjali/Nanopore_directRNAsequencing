nextflow.enable.dsl=2

process FASTQC {
    tag "${sample_id}"
    publishDir "${params.out_dir}/fastqc", mode: 'copy'
  

    input:
    tuple val(sample_id), path(fastq_file)

    output:
    tuple val(sample_id), path("*.zip"), emit: fastqc_zip
    tuple val(sample_id), path("*.html"), emit: fastqc_html

    script:
    """
    fastqc $fastq_file --outdir .
    BASENAME=\$(basename ${fastq_file} .fastq)
    mv \${BASENAME}_fastqc.zip ${sample_id}_fastqc.zip
    mv \${BASENAME}_fastqc.html ${sample_id}_fastqc.html
    """
}