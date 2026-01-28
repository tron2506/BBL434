#!/usr/bin/env python3

import sys
import math
from collections import Counter
import matplotlib.pyplot as plt

K = 8
WINDOW_SIZE = 5000
STEP = 500
EPS = 1e-9


def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip().upper())
    return "".join(seq)


def count_kmers(seq, k):
    counts = Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if "N" not in kmer:
            counts[kmer] += 1
    return counts


def normalize(counts):
    total = sum(counts.values())
    return {k: v / total for k, v in counts.items()}


def gc_skew(window):
    g = window.count("G")
    c = window.count("C")
    if g + c == 0:
        return 0.0
    return (g - c) / (g + c)


def main(fasta_file):
    genome = read_fasta(fasta_file)
    n = len(genome)

    global_counts = count_kmers(genome, K)
    global_freq = normalize(global_counts)

    positions = []
    max_enrichment = []
    cumulative_gc = []

    gc_cum = 0.0

    for start in range(0, n - WINDOW_SIZE + 1, STEP):
        window = genome[start:start + WINDOW_SIZE]

        window_counts = count_kmers(window, K)
        window_freq = normalize(window_counts)

        enrich = []
        for kmer, wf in window_freq.items():
            gf = global_freq.get(kmer, 0.0)
            enrich.append(math.log2((wf + EPS) / (gf + EPS)))

        positions.append(start + WINDOW_SIZE // 2)
        max_enrichment.append(max(enrich) if enrich else 0.0)

        skew = gc_skew(window)
        gc_cum += skew
        cumulative_gc.append(gc_cum)

    # ORI prediction
    ori_index = cumulative_gc.index(min(cumulative_gc))
    ori_position = positions[ori_index]

    print(f"Genome length: {n}")
    print(f"Total windows: {len(positions)}")
    print(f"Predicted ORI position: {ori_position}")

    # Plot 1: k-mer enrichment
    plt.figure(figsize=(12, 4))
    plt.plot(positions, max_enrichment)
    plt.xlabel("Genome position")
    plt.ylabel("Max log2 k-mer enrichment")
    plt.title("k-mer Enrichment Signal")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("script4image1.png", dpi=300)
    plt.close()

    # Plot 2: cumulative GC skew
    plt.figure(figsize=(12, 4))
    plt.plot(positions, cumulative_gc)
    plt.axvline(ori_position, linestyle="--")
    plt.xlabel("Genome position")
    plt.ylabel("Cumulative GC skew")
    plt.title("Cumulative GC Skew (ORI at minimum)")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("script4image2.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} genome.fasta")
        sys.exit(1)

    with open("output4.txt", "w") as f:
        sys.stdout = f
        main(sys.argv[1])
