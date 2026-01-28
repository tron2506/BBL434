import sys

def read_fasta(fasta_file):
    count = 0
    with open(fasta_file, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                count = count + 1
    return count


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} seq.fa")
        sys.exit(1)

    fasta_file = sys.argv[1]
    count = read_fasta(fasta_file)
    print(count)