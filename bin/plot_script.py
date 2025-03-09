#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
import pandas as pd

def plot_genome_coverage(input_files, output_file):
    """Reads multiple mosdepth outputs and generates genome coverage plots."""
    
    plt.figure(figsize=(10, 6))
    
    for file in input_files:
        try:
            df = pd.read_csv(file, sep="\t", names=["Depth", "Fraction"])
        except Exception as e:
            print(f"Error reading coverage file {file}: {e}")
            sys.exit(1)
        
        plt.plot(df["Depth"], df["Fraction"], marker="o", linestyle="-", label=file)

    plt.xlabel("Coverage Depth")
    plt.ylabel("Fraction of Genome")
    plt.title("Genome Coverage Distribution")
    plt.legend()
    plt.grid()
    plt.savefig(output_file, dpi=300)
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 plot_script.py <output_file> <input_files...>")
        sys.exit(1)

    output_file = sys.argv[1]
    input_files = sys.argv[2:]

    plot_genome_coverage(input_files, output_file)
