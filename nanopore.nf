#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Include modules
include { FASTQC } from './modules/qc/fastqc/main.nf'
include { NANOFILT } from './modules/qc/nanofilt/main.nf'
include { FASTQC2 } from './modules/qc/fastqc/main2.nf'
include { MULTIQC } from './modules/qc/multiqc/main.nf'

// Include subworkflows
include { FASTQC_ANALYSIS } from './modules/qc/fastqc/fastqc_analysis.nf'
include { NANOFILT_ANALYSIS } from './modules/qc/nanofilt/nanofilt_analysis.nf'
include { FASTQC_ANALYSIS2 } from './modules/qc/fastqc/fastqc_analysis2.nf'
include { MULTIQC_ANALYSIS } from './modules/qc/multiqc/multiqc_analysis.nf'


workflow {
    // Creating channels for inputs
    fastq_ch = Channel.fromPath(params.fastq_files, checkIfExists: true).ifEmpty {
        error "No FASTQ files found in ${params.fastq_files}"
    }

    //reference_ch = Channel.fromPath(params.reference_dir + "/*.fa", checkIfExists: true).ifEmpty { 
        //error "No reference genome files found in ${params.reference_dir}"

    // Run FASTQC analysis
    fastqc_results = FASTQC_ANALYSIS(fastq_ch)

    // Run MultiQC on FASTQC reports
    multiqc_results = MULTIQC_ANALYSIS(fastqc_results.fastqc_reports_zip)

    // Run Nanofilt analysis
    NANOFILT_ANALYSIS(fastq_ch)

     // Retrieve outputs from NANOFILT_ANALYSIS
    nanofilt_out = NANOFILT_ANALYSIS.out.filtered_results

    // Run FASTQC again on filtered results
    FASTQC_ANALYSIS2(nanofilt_out)
    }