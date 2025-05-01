#!/usr/bin/env python3
import argparse
from utility import get_ensemble_defect
from collections import defaultdict


# example usage: python3 eval_ed.py '((((((((((((((.(((((((.....(((((((.............))))..))).......)))))))..))))..((((((((...))))))))............)))))))))).' ./test_data/5S_Bacteria_B00612/bpp.txt
def main():
    parser = argparse.ArgumentParser(
        description="Evaluate ensemble defect to quantify the deviation of an RNA ensemble from a target structure."
    )
    parser.add_argument("struc", help="Target structure (dot-bracket or CT file).")
    parser.add_argument("bpp", help="Predicted BPP file.")
    parser.add_argument(
        "-z",
        "--zero-index",
        action="store_true",
        help="Use zero-indexed BPP files.",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.01,
        help="Probability threshold for BPPs.",
    )
    args = parser.parse_args()
    defect = get_ensemble_defect(
        args.struc, args.bpp, args.zero_index, args.threshold
    )
    print(f"Ensemble Defect: {defect:.4f}")


if __name__ == "__main__":
    main()
