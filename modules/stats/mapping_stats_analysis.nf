nextflow.enable.dsl=2

process MAPPING_STATS_ANALYSIS {
    tag "$sample_id"
    publishDir "${params.out_dir}/mapping_stats", mode: 'copy'

    input:
        tuple val(sample_id), path(sam_file), path(needle_file)

    output:
        path("${sample_id}_mapping_stats.csv")

    script:
    def output_file = "${sample_id}_mapping_stats.csv"

    """
    python3 -m pip install --no-cache-dir argparse pysam pandas
    python3 ${params.mapping_script} --sam ${sam_file} --needle ${needle_file} --output ${output_file}
    
    # Check if the output file was created successfully
    if [ ! -f "${output_file}" ]; then
        echo "❌ Error: Output file ${output_file} was not created!" >&2
        exit 1
    fi

    echo "✅ Mapping statistics successfully generated: ${output_file}"
    """
    }
