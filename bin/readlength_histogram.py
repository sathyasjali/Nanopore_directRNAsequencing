import argparse
import matplotlib.pyplot as plt
import pandas as pd

def plot_histogram(input_file, output_file):
    """Generates and saves histogram plot for read lengths."""
    read_lengths = pd.read_csv(input_file, header=None, names=["length"])

    plt.figure(figsize=(8, 5))
    plt.hist(read_lengths["length"], bins=50, color="blue", alpha=0.7, edgecolor="black")
    plt.xlabel("Read Length")
    plt.ylabel("Frequency")
    plt.title(f"Read Length Distribution: {input_file}")
    plt.grid(axis="y", linestyle="--", alpha=0.7)

    plt.savefig(output_file, dpi=300)
    print(f"✅ Histogram saved: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate Read Length Histograms for Multiple FASTQ Files")
    parser.add_argument("--input", nargs="+", required=True, help="List of input files containing read lengths")
    parser.add_argument("--output", nargs="+", required=True, help="List of output image files for histograms")

    args = parser.parse_args()

    # Ensure input and output lists match in length
    if len(args.input) != len(args.output):
        raise ValueError("Number of input files must match number of output files")

    for input_file, output_file in zip(args.input, args.output):
        plot_histogram(input_file, output_file)
