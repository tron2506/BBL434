#!/usr/bin/env python3

import sys
from collections import defaultdict

K = 8
T = 3
L = 1000


def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip().upper())
    return "".join(seq)


def find_clumps(genome, k, t, L):
    clumps = defaultdict(list)
    freq = defaultdict(int)

    for i in range(L - k + 1):
        kmer = genome[i:i + k]
        if "N" not in kmer:
            freq[kmer] += 1

    for kmer, count in freq.items():
        if count >= t:
            clumps[kmer].append(0)

    for start in range(1, len(genome) - L + 1):
        out_kmer = genome[start - 1:start - 1 + k]
        in_kmer = genome[start + L - k:start + L]

        if "N" not in out_kmer:
            freq[out_kmer] -= 1
            if freq[out_kmer] == 0:
                del freq[out_kmer]

        if "N" not in in_kmer:
            freq[in_kmer] += 1

        for kmer, count in freq.items():
            if count >= t:
                clumps[kmer].append(start)

    return clumps


def main(fasta_file):
    genome = read_fasta(fasta_file)
    clumps = find_clumps(genome, K, T, L)

    print(f"Found {len(clumps)} k-mers forming ({L}, {K}, {T})-clumps\n")

    for kmer in sorted(clumps):
        positions = clumps[kmer]
        print(f"{kmer} -> {len(positions)} windows, first at position {positions[0]}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} genome.fasta")
        sys.exit(1)

    with open("output2.txt", "w") as f:
        sys.stdout = f
        main(sys.argv[1])

