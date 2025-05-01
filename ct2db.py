#!/usr/bin/env python3

import argparse
from utility import parse_ct_file
from db_lib import get_db_structure

def main():
    parser = argparse.ArgumentParser(description="Parse a CT file and output the dot-bracket structure.")
    parser.add_argument("ct_file", type=str, help="Path to the input CT file.")
    parser.add_argument("-s", "--seq", action="store_true", help="If set, also output the sequence.")
    args = parser.parse_args()

    try:
        seq, paired_pos, _ = parse_ct_file(args.ct_file)
        seq_length = len(seq)
        db = get_db_structure(seq_length, paired_pos, one_based_index=False)
        print(db)
        if args.seq:
            print(seq)
    except Exception as e:
        print(f"Error parsing CT file: {e}")
        return

if __name__ == "__main__":
    main()
    