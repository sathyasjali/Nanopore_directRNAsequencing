nextflow.enable.dsl=2

process MOSDEPTH {
    tag "$sample_id"
    publishDir "${params.out_dir}/mosdepth", mode: 'copy'

    input:
        tuple val(sample_id), path(sorted_bam), path(bam_index)

    output:
        path "${sample_id}_mosdepth_summary.txt", emit: summary
        path "${sample_id}_mosdepth.global.dist.txt", emit: dist_data

    script:
    """
    mosdepth --by 1000 ${sample_id}_mosdepth ${sorted_bam}
    mv ${sample_id}_mosdepth.summary.txt ${sample_id}_mosdepth_summary.txt
    mv ${sample_id}_mosdepth.global.dist.txt ${sample_id}_mosdepth.global.dist.txt
    """
}
