#!/usr/bin/env python3
import argparse
import os
from utility import calculate_pair_sequence_identity

from utility import get_seq_from_file

def read_sequence(input_str):
    """
    Determines if the input is a filepath or a raw sequence and returns the sequence string.
    Supports plain sequence strings and single-sequence FASTA files.
    """
    if os.path.isfile(input_str):
        with open(input_str, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
            seq = []
            for line in lines:
                if line.startswith('>'):
                    continue
                seq.append(line)
            if not seq:
                raise ValueError(f"No sequence data found in file: {input_str}")
            return ''.join(seq)
    assert set(input_str.upper()).issubset(set("AUCGT-")), f"Invalid sequence input: {input_str}"
    return input_str

def main():
    parser = argparse.ArgumentParser(
        description="Evaluate pairwise sequence identity among multiple sequences."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="Aligned sequences or file paths (dot-bracket files or direct sequences)."
    )
    parser.add_argument("-a", "--avg", action="store_true", help="Only print average sequence identity.")
    parser.add_argument("-m", "--mat", action="store_true", help="Only print the pairwise identity matrix.")
    args = parser.parse_args()

    sequences = []
    if len(args.inputs) == 1 and os.path.isfile(args.inputs[0]):
        with open(args.inputs[0], 'r') as f:
            current_seq = []
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_seq:
                        sequences.append(''.join(current_seq))
                        current_seq = []
                else:
                    current_seq.append(line)
            if current_seq:
                sequences.append(''.join(current_seq))
    else:
        sequences = [read_sequence(inp) for inp in args.inputs]
    if len(sequences) < 2:
        print("[Error] Need at least two sequences for comparison.")
        return
    n = len(sequences)

    total_identity = 0
    count = 0
    if not args.avg:
        print("Pairwise Sequence Identity Matrix (%):")
    for i in range(n):
        for j in range(i + 1, n):
            try:
                identity = calculate_pair_sequence_identity(sequences[i], sequences[j])
                if not args.avg:
                    print(f"(Seq{i+1}, Seq{j+1}): {identity:.2f}%")
                total_identity += identity
                count += 1
            except ValueError as e:
                print(f"[Error] Seq{i+1} vs Seq{j+1}: {e}")
    if count > 0 and not args.mat:
        average_identity = total_identity / count
        print(f"Average Sequence Identity: {average_identity:.2f}%")

if __name__ == "__main__":
    main()
    