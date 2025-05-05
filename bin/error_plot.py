import argparse
import matplotlib.pyplot as plt
import pandas as pd

def plot_error_rate(vcf_file, coverage_file, output_file, sample_id):
    # Load error rate (VCF parsing - assuming simple error metric)
    vcf_data = pd.read_csv(vcf_file, comment='#', sep='\t', header=None)
    error_rate = vcf_data.shape[0]  # Number of variants as a proxy for error

    # Load coverage data
    coverage_data = pd.read_csv(coverage_file, header=None, sep='\t', names=["region", "coverage"])

    # Plot
    plt.figure(figsize=(8, 5))
    plt.scatter(coverage_data["coverage"], [error_rate] * len(coverage_data), alpha=0.6, color="red")
    plt.xlabel("Coverage")
    plt.ylabel("Error Rate (Variant Count)")
    plt.title(f"Error Rate vs Coverage for {sample_id}")

    # Save plot
    plt.savefig(output_file, dpi=300)
    print(f"✅ Error Plot saved: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate Error Rate vs Coverage Plot")
    parser.add_argument("--vcf", required=True, help="Input VCF file")
    parser.add_argument("--coverage", required=True, help="Coverage file")
    parser.add_argument("--output", required=True, help="Output plot file")
    parser.add_argument("--sample_id", required=True, help="Sample ID for labeling")

    args = parser.parse_args()
    plot_error_rate(args.vcf, args.coverage, args.output, args.sample_id)
