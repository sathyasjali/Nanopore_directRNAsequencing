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
include { READLENGTH_HISTOGRAM } from './modules/plots/readlength_histogram/main.nf'
include { ERROR_PLOT } from './modules/plots/error_plot/main.nf'



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
include { READLENGTH_HISTOGRAM_ANALYSIS } from './modules/plots/readlength_histogram/readlength_histogram_analysis.nf'
include { ERROR_PLOT_ANALYSIS } from './modules/plots/error_plot/error_plot_analysis.nf'
include { MAPPING_STATS_ANALYSIS } from './modules/stats/mapping_stats_analysis.nf'
include { MERGE_CSV_FILES } from './modules/stats/merge_csv_files.nf'

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

     // Extract VCF output from BCFTOOLS
    error_vcf_ch = bcftools_results.vcf_output

    // Create BAM index channel
    bam_index_ch = samtools_results.bam_index

    // Run EMBOSS Needle alignment
    needle_results = NEEDLE_ANALYSIS(reference_ch, bcftools_results.consensus_sequence)

    // Pass correct input channels
    mosdepth_results = MOSDEPTH_ANALYSIS(samtools_results.sorted_bam, samtools_results.bam_index)

    // Now pass results to PLOTTING_ANALYSIS
    plot_results = PLOTTING_ANALYSIS(
        mosdepth_results.global_dist.map { it[1] }.collect() // Extract file path
        )

    // Inside workflow block
    histogram_plots = READLENGTH_HISTOGRAM_ANALYSIS(nanofilt_out)

    histogram_plots.view { it -> "Read Length Histogram Generated: ${it}" } .println()

    // Extract Coverage File
    coverage_ch = mosdepth_results.global_dist.map { it[1] }

    // Ensure error_vcf_ch is properly formatted with sample_id and vcf_file
    error_vcf_ch = error_vcf_ch.map { sample_id, vcf_file -> tuple(sample_id, vcf_file) }

    // Run ERROR_PLOT_ANALYSIS with **two separate input channels**
    error_plot_results = ERROR_PLOT_ANALYSIS(error_vcf_ch, coverage_ch)

    // Display error plot outputs
    error_plot_results.view { it -> "Error Plot Generated: ${it}" }

    // Extract necessary channels
    sam_ch = minimap2_results.map { id, sam_file -> tuple(id, sam_file) }
    needle_ch = needle_results.map { id, needle_file -> tuple(id, needle_file) }

    // Ensure correct pairing by sample ID using join()
    combined_stats_ch = sam_ch.join(needle_ch)

    // Debug: Verify the output before running the process
    combined_stats_ch.view()
    sam_ch.view()
    needle_ch.view()
    // Run mapping statistics analysis
    mapping_stats_results = MAPPING_STATS_ANALYSIS(combined_stats_ch)

    // Output mapping statistics results
    mapping_stats_results.view { it -> "Mapping Stats Generated: ${it}" }


    // Collect all generated CSVs
    merged_csv = MERGE_CSV_FILES(mapping_stats_results)

    // Save final output
    merged_csv.view()
}