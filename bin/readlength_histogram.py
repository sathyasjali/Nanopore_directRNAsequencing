import argparse
import matplotlib.pyplot as plt
import pandas as pd

def plot_histogram(input_file, output_file):
    # Load read lengths
    read_lengths = pd.read_csv(input_file, header=None, names=["length"])

    # Plot histogram
    plt.figure(figsize=(8, 5))
    plt.hist(read_lengths["length"], bins=50, color="blue", alpha=0.7, edgecolor="black")
    plt.xlabel("Read Length")
    plt.ylabel("Frequency")
    plt.title("Read Length Distribution")
    plt.grid(axis="y", linestyle="--", alpha=0.7)

    # Save plot
    plt.savefig(output_file, dpi=300)
    print(f"✅ Histogram saved: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate Read Length Histogram")
    parser.add_argument("--input", required=True, help="Input file containing read lengths")
    parser.add_argument("--output", required=True, help="Output image file for histogram")
    
    args = parser.parse_args()
    plot_histogram(args.input, args.output)
