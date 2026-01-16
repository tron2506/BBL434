#run the code in terminal - python3 script2.py fasta_file.fa

import sys

def read_fasta(fasta_file):
    seq = []
    with open(fasta_file, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line)
    return "".join(seq)

def kmer_count(sequence, k=3):
    kmer_dict = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmer_dict[kmer] = kmer_dict.get(kmer, 0) + 1
    return kmer_dict

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} seq.fa")
        sys.exit(1)

    fasta_file = sys.argv[1]
    sequence = read_fasta(fasta_file)

    if len(sequence) < 3:
        print({})
    else:
        kmer_dict = kmer_count(sequence, k=3)
        print(kmer_dict)
