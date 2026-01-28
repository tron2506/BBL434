#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

WINDOW_SIZE = 5000
STEP = 500


def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip().upper())
    return "".join(seq)


def gc_skew(window):
    g = window.count("G")
    c = window.count("C")
    if g + c == 0:
        return 0.0
    return (g - c) / (g + c)


def main(fasta_file):
    genome = read_fasta(fasta_file)

    positions = []
    cumulative = []

    cum_sum = 0.0

    for start in range(0, len(genome) - WINDOW_SIZE + 1, STEP):
        window = genome[start:start + WINDOW_SIZE]
        skew = gc_skew(window)

        cum_sum += skew
        positions.append(start + WINDOW_SIZE // 2)
        cumulative.append(cum_sum)

    print(f"Genome length: {len(genome)}")
    print(f"Total windows: {len(cumulative)}")

    plt.figure(figsize=(12, 5))
    plt.plot(positions, cumulative, linewidth=1)
    plt.xlabel("Genome position")
    plt.ylabel("Cumulative GC skew")
    plt.title("Cumulative GC Skew Plot")
    plt.grid(alpha=0.3)
    plt.tight_layout()

    plt.savefig("script3image.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} genome.fasta")
        sys.exit(1)

    with open("output3.txt", "w") as f:
        sys.stdout = f
        main(sys.argv[1])
