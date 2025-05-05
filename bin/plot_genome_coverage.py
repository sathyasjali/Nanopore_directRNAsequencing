import sys
import matplotlib.pyplot as plt
import pandas as pd

def plot_genome_coverage(coverage_file, output_file):
    """Reads mosdepth output and generates a genome coverage plot."""
    try:
        df = pd.read_csv(coverage_file, sep="\t", names=["Depth", "Fraction"])
    except Exception as e:
        print(f"Error reading coverage file: {e}")
        sys.exit(1)

    plt.figure(figsize=(8, 5))
    plt.plot(df["Depth"], df["Fraction"], marker="o", linestyle="-", color="red")
    plt.xlabel("Coverage Depth")
    plt.ylabel("Fraction of Genome")
    plt.title("Genome Coverage Distribution")
    plt.grid()
    plt.savefig(output_file, dpi=300)
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 plot_genome_coverage.py <coverage_file> <output_file>")
        sys.exit(1)

    plot_genome_coverage(sys.argv[1], sys.argv[2])