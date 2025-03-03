#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//Include modules
include { FASTQC } from './modules/qc/fastqc/main.nf'
include { NANOFILT } from './modules/qc/nanofilt/main.nf'
include { MULTIQC } from './modules/qc/multiqc/main.nf'
include { MINIMAP2 } from './modules/align/minimap2/main.nf'
include { SAMTOOLS_BCFTOOLS } from './modules/align/samtools_bcftools/main.nf'
include { SAMTOOLS_INDEX } from './modules/align/samtools_bcftools/main.nf'


workflow {
    // creating channels for inputs input
    fastq_ch = Channel.fromPath("${params.fastq_files}")
    reference_ch = Channel.fromPath("${params.reference_dir}/*.fa").ifEmpty { error "No reference genome files found in ${params.reference_dir}" }

    // Step1: Run FastQC
    fastqc_results = fastq_ch.map { file -> tuple(file.baseName, file) } | FASTQC

    // Step2: Run NanoFilt
    nanofilt_results = fastq_ch.map { file -> tuple(file.baseName, file) } | NANOFILT

    // Step3: Run MultiQC
    multiqc_report = MULTIQC(fastqc_results.map { it[1] }, nanofilt_results.map { it[1] })

    // Step 4: Run Minimap2 with correct channel combination
    //minimap2_results = nanofilt_results
    //.map { sample -> tuple(sample[0], sample[1]) }  // Ensure tuples
    // .combine(reference_ch)  // Use combine instead of cross
    //| MINIMAP2

     // Step 4: Run Minimap2 with correct channel combination
    minimap2_results = nanofilt_results
        .map { sample -> tuple(sample[0], sample[1]) }  // Ensure tuples
        .combine(reference_ch.map { ref -> tuple(ref) })  // Ensure reference is a tuple
        | MINIMAP2
    
    // Step5: Run Samtools and Bcftools
    samtools_bcftools_results = minimap2_results
        .map { tuple(it[0], it[1], it[2]) }  // Ensure correct tuple structure
        | SAMTOOLS_BCFTOOLS

    // Step6: Run Samtools Index
    samtools_index_results = samtools_bcftools_results
        .map { tuple(it[0], it[1], reference_ch.first()) }  // Ensure correct tuple structure
        | SAMTOOLS_INDEX
}