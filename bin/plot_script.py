import argparse
import os
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate genome coverage plots.")
    parser.add_argument("--global_dist", required=True, help="Path to global distribution file")
    parser.add_argument("--region_dist", required=True, help="Path to region distribution file")
    parser.add_argument("--summary", required=True, help="Path to summary file")
    parser.add_argument("--per_base", required=True, help="Path to per-base coverage file")
    parser.add_argument("--per_region", required=True, help="Path to per-region coverage file")
    parser.add_argument("--threshold", required=True, help="Path to threshold coverage file")
    parser.add_argument("--output", required=True, help="Output image file")

    args = parser.parse_args()

    # Debug: Print received arguments
    print("Received arguments:")
    for arg, value in vars(args).items():
        print(f"{arg}: {value}")

    # Check if files exist
    for file_path in [args.global_dist, args.region_dist, args.summary, args.per_base, args.per_region, args.threshold]:
        if not os.path.exists(file_path):
            print(f"ERROR: Missing input file -> {file_path}", file=sys.stderr)
            sys.exit(1)

    return args

if __name__ == "__main__":
    parse_arguments()
