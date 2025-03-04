#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//Include modules
include { FASTQC } from './modules/qc/fastqc/main.nf'
include { FASTQC2 } from './modules/qc/fastqc/main2.nf'
include { NANOFILT } from './modules/qc/nanofilt/main.nf'
include { MULTIQC } from './modules/qc/multiqc/main.nf'
include { MINIMAP2 } from './modules/align/minimap2/main.nf'
include { SAMTOOLS_BCFTOOLS } from './modules/align/samtools_bcftools/main.nf'
include { SAMTOOLS_INDEX } from './modules/align/samtools_bcftools/main.nf'


workflow {
    // creating channels for inputs input
    fastq_ch = Channel.fromPath("${params.fastq_files}")
    reference_ch = Channel.fromPath("${params.reference_dir}/*.fa").ifEmpty { error "No reference genome files found in ${params.reference_dir}" }

    // Step1: Run FastQC on raw reads.
    fastqc_results = fastq_ch.map { file -> tuple(file.baseName, file) } | FASTQC

    // Step2: Run NanoFilt.
    nanofilt_results = fastq_ch.map { file -> tuple(file.baseName, file) } | NANOFILT

    // Step3: Run FastQC on NanoFilt filtered results
    fastqc_nanofilt_results = nanofilt_results.map { sample -> tuple(sample[0], sample[1]) } | FASTQC2
    
    // Step4: Run MultiQC on FastQC results and NanoFilt results
    multiqc_report = MULTIQC(fastqc_results.flatMap { it[1] }, fastqc_nanofilt_results.flatMap { it[1] })

    // Step5: Run Minimap2 with correct channel combination
    minimap2_results = nanofilt_results
            .map { sample -> tuple(sample[0], sample[1]) }  // Ensure tuples
            .combine(reference_ch)  // Use combine instead of cross
            | MINIMAP2
    
    // Step6: Convert SAM to BAM, sort, and index using Samtools and BCFtools
    samtools_bcftools_results = minimap2_results
        .map { tuple(it[0], it[1]) }  // Removed extra element
        | SAMTOOLS_BCFTOOLS

    // // Step6: Run Samtools Index
    // samtools_index_results = samtools_bcftools_results
    //     .map { tuple(it[0], it[1], reference_ch.first()) }  // Ensure correct tuple structure
    //     | SAMTOOLS_INDEX

    // Step7: Run Samtools Index with correct reference handling
    samtools_index_results = samtools_bcftools_results
        .combine(reference_ch)  // Ensure reference genome is included properly
        | SAMTOOLS_INDEX
}