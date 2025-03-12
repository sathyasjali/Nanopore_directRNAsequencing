import pysam
import argparse
import pandas as pd
import re

def parse_sam(sam_file):
    """Parses the SAM file to extract alignment statistics."""
    try:
        samfile = pysam.AlignmentFile(sam_file, "r")
    except Exception as e:
        print(f"❌ Error opening SAM file: {sam_file} - {str(e)}")
        return {
            "Total Reads": 0,
            "Mapped Reads": 0,
            "Unmapped Reads": 0,
            "Mapping Percentage": 0.0
        }

    total_reads = 0
    mapped_reads = 0
    unmapped_reads = 0

    for read in samfile:
        total_reads += 1
        if read.is_unmapped:
            unmapped_reads += 1
        else:
            mapped_reads += 1

    mapping_percentage = (mapped_reads / total_reads) * 100 if total_reads > 0 else 0.0

    return {
        "Total Reads": total_reads,
        "Mapped Reads": mapped_reads,
        "Unmapped Reads": unmapped_reads,
        "Mapping Percentage": mapping_percentage
    }

def parse_needle(needle_file):
    """Parses the Needle output file to extract sequence identity."""
    identity_values = []
    
    with open(needle_file, 'r') as file:
        for line in file:
            match = re.search(r"(\d+\.\d+)%", line)  # Extract numbers like 97.8%
            if match:
                try:
                    identity = float(match.group(1))  # Extract just the number
                    identity_values.append(identity)
                except ValueError:
                    print(f"⚠️ Warning: Could not convert {match.group(1)} to float")

    if identity_values:
        return {"Average Identity": sum(identity_values) / len(identity_values)}
    else:
        print(f"⚠️ Warning: No valid identity percentage found in {needle_file}")
        return {"Average Identity": 0.0}  # Prevent errors in main()

def main():
    parser = argparse.ArgumentParser(description="Generate Mapping Statistics from SAM and Needle Files")
    parser.add_argument("--sam", required=True, help="Input SAM file")
    parser.add_argument("--needle", required=True, help="Input Needle alignment file")
    parser.add_argument("--output", required=True, help="Output CSV file")
    
    args = parser.parse_args()
    
    sam_stats = parse_sam(args.sam)
    needle_stats = parse_needle(args.needle)

    # Merge results into a single dataframe
    stats = {**sam_stats, **needle_stats}
    df = pd.DataFrame([stats])

    # Save to CSV
    df.to_csv(args.output, index=False)
    print(f"✅ Mapping Statistics saved to: {args.output}")

if __name__ == "__main__":
    main()
