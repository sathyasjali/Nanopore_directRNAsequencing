process READLENGTH_HISTOGRAM {
    tag "$sample_id"
    publishDir "${params.out_dir}/readlength", mode: 'copy'

    input:
        tuple val(sample_id), path(filtered_fastq)

    output:
        tuple val(sample_id), path("${sample_id}_read_lengths.txt"), path("${sample_id}_read_length_histogram.png")

    script:
    """
    # Ensure the file exists before processing
    if [ ! -f "${filtered_fastq}" ]; then
        echo "Error: File ${filtered_fastq} does not exist!" >&2
        exit 1
    fi

    echo "Processing FASTQ file: ${filtered_fastq}"

    # Extract read lengths
    awk '(NR%4==2){print length($1)}' "${filtered_fastq}" > ${sample_id}_read_lengths.txt

    # Ensure required Python packages are installed
    python3 -m pip install --no-cache-dir --upgrade pip
    python3 -m pip install --no-cache-dir matplotlib pandas

    # Run the Python script for plotting
    python3 ${params.script_dir} --input ${sample_id}_read_lengths.txt 
                                 --output ${sample_id}_read_length_histogram.png
                                 --sample_id ${sample_id}
    """
}
