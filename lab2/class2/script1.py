import matplotlib.pyplot as plt
from collections import Counter
import random

# --- Configuration ---
K = 8
WINDOW_SIZE = 5000
STEP = 500
TOP_N_KMERS = 3  # Number of top enriched k-mers to plot/report

def generate_sample_genome(length=50000):
    """Generates a synthetic genome with an artificial signal for demonstration."""
    random.seed(42)
    dna = ''.join(random.choice('ATGC') for _ in range(length))
    # Inject an artificial "DnaA box" signal at position 25,000
    signal = 'ATGATCAT'
    dna_list = list(dna)
    for i in range(25000, 27000, 50): 
        dna_list[i:i+K] = list(signal)
    return "".join(dna_list)

def get_kmer_counts(sequence, k):
    """Returns a Counter of all k-mers in a specific sequence string."""
    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    return Counter(kmers)

def main():
    # 1. Load Genome (using synthetic data for this example)
    genome = generate_sample_genome()
    print(f"Processing genome of length: {len(genome)} bp...")

    # 2. Sliding Window Analysis
    window_stats = [] # Stores (position, Counter_object)
    
    for i in range(0, len(genome) - WINDOW_SIZE + 1, STEP):
        window_seq = genome[i : i + WINDOW_SIZE]
        counts = get_kmer_counts(window_seq, K)
        window_stats.append((i, counts))

    # 3. Identify Globally Enriched K-mers
    # We find which k-mers have the highest total frequency across all windows
    global_counter = Counter()
    for _, counts in window_stats:
        global_counter.update(counts)
    
    top_kmers = [kmer for kmer, count in global_counter.most_common(TOP_N_KMERS)]

    # 4. Save Text Output
    with open("output1.txt", "w") as f:
        f.write("--- ORI Signal Checker Results ---\n")
        f.write(f"Parameters: k={K}, window={WINDOW_SIZE}, step={STEP}\n\n")
        f.write("Top Enriched K-mers:\n")
        for kmer in top_kmers:
            f.write(f"Sequence: {kmer} | Total Occurrences: {global_counter[kmer]}\n")

    # 5. Plotting and Saving Images
    for idx, target_kmer in enumerate(top_kmers, 1):
        positions = [item[0] for item in window_stats]
        frequencies = [item[1].get(target_kmer, 0) for item in window_stats]

        plt.figure(figsize=(10, 5))
        plt.plot(positions, frequencies, color='teal', linewidth=2, marker='o', markersize=4)
        
        plt.title(f"Enrichment Pattern for K-mer: {target_kmer}", fontsize=14)
        plt.xlabel("Genome Position (bp)", fontsize=12)
        plt.ylabel(f"Occurrences (Window={WINDOW_SIZE})", fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.6)
        
        # Highlight the peak
        max_freq = max(frequencies)
        max_pos = positions[frequencies.index(max_freq)]
        plt.annotate(f'Peak: {max_freq}', xy=(max_pos, max_freq), xytext=(max_pos+2000, max_freq),
                     arrowprops=dict(facecolor='black', shrink=0.05))

        image_name = f"output1image{idx}.png"
        plt.tight_layout()
        plt.savefig(image_name)
        plt.close()
        print(f"Generated {image_name}")

    print("Analysis complete. Check output1.txt and image files.")

if __name__ == "__main__":
    main()
