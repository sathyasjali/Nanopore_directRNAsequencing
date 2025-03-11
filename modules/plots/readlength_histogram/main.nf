nextflow.enable.dsl=2

process READLENGTH_HISTOGRAM {
    tag "$sample_id"
    publishDir "${params.out_dir}/readlength", mode: 'copy'

    input:
        tuple val(sample_id), path(filtered_fastq)

    output:
        path "read_lengths.txt"
        //path "read_length_histogram.png"

    script:
    """
    # Ensure the file exists before processing
    if [ ! -f "${filtered_fastq}" ]; then
        echo "Error: File ${filtered_fastq} does not exist!" >&2
        exit 1
    fi

    echo "Processing FASTQ file: ${filtered_fastq}"

    awk '(NR%4==2){print length(\$1)}' "${filtered_fastq}" > read_lengths.txt
    # Run the Python script for plotting
    #python3 ${params.script_dir}/readlength_histogram.py --input read_lengths.txt --output read_length_histogram.png
    """
}
