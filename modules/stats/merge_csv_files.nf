nextflow.enable.dsl=2

process MERGE_CSV_FILES {
    tag "merge_csv"
    publishDir "${params.out_dir}/mapping_stats", mode: 'copy'

    input:
        path(csv_files)

    output:
        path("combined_mapping_stats.csv")

    script:
    """
    python3 -m pip install --no-cache-dir pandas
    python3 -c '
import pandas as pd
import glob

# Get all CSV files
csv_files = glob.glob("*.csv")
if not csv_files:
    print("❌ No CSV files found to combine.")
    exit(1)

# Merge all CSVs
df_list = [pd.read_csv(f) for f in csv_files]
df = pd.concat(df_list, ignore_index=True)

df.to_csv("combined_mapping_stats.csv", index=False)
print("✅ Combined CSV file generated: combined_mapping_stats.csv")
    '
    """
}
