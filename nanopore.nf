#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Include modules
include { FASTQC } from './modules/qc/fastqc/main.nf'
include { NANOFILT } from './modules/qc/nanofilt/main.nf'
include { FASTQC2 } from './modules/qc/fastqc/main2.nf'
include { MULTIQC } from './modules/qc/multiqc/main.nf'
include { MULTIQC2 } from './modules/qc/multiqc/main2.nf'
include { MINIMAP2 } from './modules/align/minimap2/main.nf'
include { SAMTOOLS } from './modules/align/samtools/main.nf'
include { BCFTOOLS } from './modules/variant/bcftools/main.nf'
include { EMBOSS} from './modules/alignment/emboss/main.nf'
include { MOSDEPTH } from './modules/plots/mosdepth/main.nf'
include { PLOTTING } from './modules/plots/main.nf'


// Include subworkflows
include { FASTQC_ANALYSIS } from './modules/qc/fastqc/fastqc_analysis.nf'
include { NANOFILT_ANALYSIS } from './modules/qc/nanofilt/nanofilt_analysis.nf'
include { FASTQC_ANALYSIS2 } from './modules/qc/fastqc/fastqc_analysis2.nf'
include { MULTIQC_ANALYSIS } from './modules/qc/multiqc/multiqc_analysis.nf'
include { MULTIQC_ANALYSIS2 } from './modules/qc/multiqc/multiqc_analysis2.nf'
include { MINIMAP2_ANALYSIS } from './modules/align/minimap2/minimap2_analysis.nf'
include { SAMTOOLS_ANALYSIS } from './modules/align/samtools/samtools_analysis.nf'
include { BCFTOOLS_ANALYSIS } from './modules/variant/bcftools/bcftools_analysis.nf'
include { NEEDLE_ANALYSIS } from './modules/alignment/emboss/needle_analysis.nf'
include { MOSDEPTH_ANALYSIS } from './modules/plots/mosdepth/mosdepth_analysis.nf'
include { PLOTTING_ANALYSIS } from './modules/plots/plotting_analysis.nf'

workflow {
    // Creating channels for inputs
    fastq_ch = Channel.fromPath(params.fastq_files, checkIfExists: true).ifEmpty {
        error "No FASTQ files found in ${params.fastq_files}"
    }

    reference_ch = Channel
        .fromPath("${params.reference_dir}/*.fa", checkIfExists: true)
        .ifEmpty { error "No reference genome files found in ${params.reference_dir}" }

    // Run FASTQC analysis
    fastqc_results = FASTQC_ANALYSIS(fastq_ch)

    // Run MultiQC on FASTQC reports
    multiqc_results = MULTIQC_ANALYSIS(fastqc_results.fastqc_reports_zip)

    // Run Nanofilt analysis
    NANOFILT_ANALYSIS(fastq_ch)

    // Retrieve outputs from NANOFILT_ANALYSIS
    nanofilt_out = NANOFILT_ANALYSIS.out.filtered_results

    // Run FASTQC again on filtered results
    fastqc_filtered_results = FASTQC_ANALYSIS2(nanofilt_out)

    // MULTIQC_ANALYSIS2(fastqc_filtered_results.fastqc_filtered_reports_zip)

    // Run Minimap2 alignment on filtered FASTQ files
    minimap2_results = MINIMAP2_ANALYSIS(nanofilt_out, reference_ch)

    // Samtools Processing
    samtools_results = SAMTOOLS_ANALYSIS(minimap2_results.sam_output, reference_ch)

    // BCFtools Variant Calling
    bcftools_results = BCFTOOLS_ANALYSIS(samtools_results.sorted_bam, reference_ch)

    // Create BAM index channel
    bam_index_ch = samtools_results.bam_index

    // Run EMBOSS Needle alignment
    needle_results = NEEDLE_ANALYSIS(reference_ch, bcftools_results.consensus_sequence)

    
    // ✅ Pass correct input channels
    mosdepth_results = MOSDEPTH_ANALYSIS(samtools_results.sorted_bam, samtools_results.bam_index)
    // ✅ Now pass results to PLOTTING_ANALYSIS
    plot_results = PLOTTING_ANALYSIS(
        mosdepth_results.global_dist.map { it[0] },  // Extract sample_id (First input)
        mosdepth_results.global_dist.map { it[1] },  // Extract file path
        mosdepth_results.region_dist.map { it[1] },
        mosdepth_results.summary.map { it[1] },
        mosdepth_results.per_base_coverage.map { it[1] },
        mosdepth_results.per_region_coverage.map { it[1] },
        mosdepth_results.threshold_coverage.map { it[1] },
        )

    // Use plot_results to save or view the plots
    plot_results.view { it -> "Plotting results: ${it}" }
}