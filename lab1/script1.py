#!/usr/bin/env python3

###### ChatGPT Code | Vibe Coding Tutorial ######

import sys

def fasta_length(fasta_file):
    seq_len = 0
    with open(fasta_file, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue  # skip header
            seq_len += len(line)
    return seq_len

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} seq.fa")
        sys.exit(1)

    fasta_file = sys.argv[1]
    length = fasta_length(fasta_file)
    print(f"hi")
    print(length)
    