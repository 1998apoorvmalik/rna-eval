#!/usr/bin/env python3
import argparse
from utility import evaluate, parse_ct_file


def main():
    parser = argparse.ArgumentParser(description="Evaluate the similarity between two RNA secondary structures using precision, sensitivity, F1 score, and structural distance.")
    parser.add_argument("struc1", type=str, help="First structure (dot-bracket or CT file).")
    parser.add_argument("struc2", type=str, help="Second structure (dot-bracket or CT file).")
    parser.add_argument("-s", "--slip", action="store_true", help="Allow one-nucleotide slip while evaluating.")
    args = parser.parse_args()

    is_ct1 = args.struc1.endswith(".ct")
    is_ct2 = args.struc2.endswith(".ct")

    if is_ct1:
        _, struc1_paired, struc1_unpaired = parse_ct_file(args.struc1)
    else:
        struc1_paired = struc1_unpaired = None

    if is_ct2:
        _, struc2_paired, struc2_unpaired = parse_ct_file(args.struc2)
    else:
        struc2_paired = struc2_unpaired = None

    if is_ct1 and is_ct2:
        precision, sensitivity, f1, sd = evaluate(
            struc1_paired_pos_tuple=struc1_paired,
            struc1_unpaired_pos=struc1_unpaired,
            struc2_paired_pos_tuple=struc2_paired,
            struc2_unpaired_pos=struc2_unpaired,
            allow_slip=args.slip
        )
    elif is_ct1:
        precision, sensitivity, f1, sd = evaluate(
            struc1_paired_pos_tuple=struc1_paired,
            struc1_unpaired_pos=struc1_unpaired,
            struc2=args.struc2,
            allow_slip=args.slip
        )
    elif is_ct2:
        precision, sensitivity, f1, sd = evaluate(
            struc1=args.struc1,
            struc2_paired_pos_tuple=struc2_paired,
            struc2_unpaired_pos=struc2_unpaired,
            allow_slip=args.slip
        )
    else:
        precision, sensitivity, f1, sd = evaluate(
            struc1=args.struc1,
            struc2=args.struc2,
            allow_slip=args.slip
        )

    print(f"Precision: {precision:.4f}")
    print(f"Sensitivity: {sensitivity:.4f}")
    print(f"F1 Score: {f1:.4f}")
    print(f"Structural Distance: {sd:.4f}")

if __name__ == "__main__":
    main()
