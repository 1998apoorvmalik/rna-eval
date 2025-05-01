#!/usr/bin/env python3

import argparse
from utility import parse_ct_file

def main():
    parser = argparse.ArgumentParser(description="Parse a CT file and output the sequence")
    parser.add_argument("ct_file", type=str, help="Path to the input CT file.")
    args = parser.parse_args()

    try:
        seq, _, _ = parse_ct_file(args.ct_file)
        print(seq)
    except Exception as e:
        print(f"Error parsing CT file: {e}")
        return

if __name__ == "__main__":
    main()
    