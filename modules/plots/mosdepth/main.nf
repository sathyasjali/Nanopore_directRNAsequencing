nextflow.enable.dsl=2

process MOSDEPTH {
    tag "$sample_id"
    publishDir "${params.out_dir}/mosdepth", mode: 'copy'

    cpus 4
    memory '8GB'

    input:
        tuple val(sample_id), path(sorted_bam), path(bam_index)

    output:
        tuple val(sample_id), path("${sample_id}.mosdepth.global.dist.txt"), emit: global_dist
        tuple val(sample_id), path("${sample_id}.mosdepth.region.dist.txt"), emit: region_dist
        tuple val(sample_id), path("${sample_id}.mosdepth.summary.txt"), emit: summary
        tuple val(sample_id), path("${sample_id}.mosdepth.thresholds.bed.gz"), emit: threshold_coverage
        tuple val(sample_id), path("${sample_id}.mosdepth.per-base.bed.gz"), emit: per_base_coverage
        tuple val(sample_id), path("${sample_id}.mosdepth.per-region.bed.gz"), emit: per_region_coverage

    script:
    """
    #!/bin/bash -ue
    echo "Running MOSDEPTH for sample: ${sample_id}"
    echo "Using sorted BAM file: ${sorted_bam}"
    echo "Using BAM index: ${bam_index}"

    mosdepth -t ${task.cpus} --by 1000 "${sample_id}_mosdepth" "${sorted_bam}"

    mv "${sample_id}_mosdepth.summary.txt" "${sample_id}_mosdepth_summary.txt"
    mv "${sample_id}_mosdepth.global.dist.txt" "${sample_id}_mosdepth_global.dist.txt"
    """
}