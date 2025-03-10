#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//Include modules
include { FASTQC } from './modules/qc/fastqc/main.nf'
include { NANOFILT } from './modules/qc/nanofilt/main.nf'
include { MULTIQC } from './modules/qc/multiqc/main.nf'
include { MINIMAP2 } from './modules/align/minimap2/main.nf'
include { SAMTOOLS_BCFTOOLS } from './modules/align/samtools_bcftools/main.nf'
include { SAMTOOLS_INDEX } from './modules/align/samtools_bcftools/main2.nf'
include { REFERENCE_INDEX } from './modules/align/reference_index/main.nf'  // NEW: Reference genome indexing


workflow {
    // creating channels for inputs input
    fastq_ch = Channel.fromPath("${params.fastq_files}")
    reference_ch = Channel.fromPath("${params.reference_dir}/*.fa").ifEmpty { error "No reference genome files found in ${params.reference_dir}" }

    // Step 1: Index Reference Genome
    indexed_reference_ch = reference_ch | REFERENCE_INDEX
    
    // Step 2: Run FastQC
    fastqc_results = fastq_ch.map { file -> tuple(file.baseName, file) } | FASTQC

    // Step 3: Run NanoFilt
    nanofilt_results = fastq_ch.map { file -> tuple(file.baseName, file) } | NANOFILT

    // Step 4: Run MultiQC
    multiqc_report = MULTIQC(fastqc_results.map { it[1] }, nanofilt_results.map { it[1] })
    
    // Step 5: Run Minimap2 with correct channel combination
    minimap2_results = nanofilt_results
        .map { sample -> tuple(sample[0], sample[1]) }  // Ensure tuples
        .combine(reference_ch.map { ref -> tuple(ref) })  // Ensure reference is a tuple
        | MINIMAP2

    // Step 6: Run Samtools/Bcftools
    samtools_bcftools_results = minimap2_results
        .combine(indexed_reference_ch)  // Pair BAM with FASTA
        | SAMTOOLS_BCFTOOLS

    // Step 7: Run Samtools Index (correct tuple structure)
    samtools_index_results = samtools_bcftools_results
        .map { tuple(it[0], it[2]) }  // Pass only sample_id and BAM
        .combine(reference_ch.map { tuple(it) }) // Pass reference separately
        | SAMTOOLS_INDEX 
}