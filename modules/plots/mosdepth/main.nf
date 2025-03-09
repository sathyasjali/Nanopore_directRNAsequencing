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
    echo "Checking input files:"
    ls -lh "${sorted_bam}" || { echo "ERROR: Sorted BAM file missing: ${sorted_bam}"; exit 1; }
    ls -lh "${bam_index}" || { echo "ERROR: BAM index missing: ${bam_index}"; exit 1; }

    echo "Using sorted BAM file: ${sorted_bam}"
    echo "Using BAM index: ${bam_index}"

    # Run mosdepth with thresholds enabled
    mosdepth -t ${task.cpus} --by 1000 --thresholds 1,5,10,20,30,50 "${sample_id}_mosdepth" "${sorted_bam}" || { echo "ERROR: mosdepth failed"; exit 1; }

    # Rename outputs to match expected Nextflow output
    mv "${sample_id}_mosdepth.mosdepth.summary.txt" "${sample_id}.mosdepth.summary.txt"
    mv "${sample_id}_mosdepth.mosdepth.global.dist.txt" "${sample_id}.mosdepth.global.dist.txt"
    mv "${sample_id}_mosdepth.mosdepth.region.dist.txt" "${sample_id}.mosdepth.region.dist.txt"
    mv "${sample_id}_mosdepth.per-base.bed.gz" "${sample_id}.mosdepth.per-base.bed.gz"
    mv "${sample_id}_mosdepth.regions.bed.gz" "${sample_id}.mosdepth.per-region.bed.gz"
    mv "${sample_id}_mosdepth.thresholds.bed.gz" "${sample_id}.mosdepth.thresholds.bed.gz" || echo "WARNING: No thresholds file generated."

    """
}