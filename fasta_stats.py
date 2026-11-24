"""Read a FASTA file and report simple statistics for each sequence.

Usage:
    python fasta_stats.py sample.fasta
"""
import sys


def read_fasta(path):
    """Return a dictionary: sequence name -> sequence (upper case letters)."""
    sequences = {}
    name = None
    with open(path) as file:
        for line in file:
            line = line.strip()
            if not line:
                continue  # skip empty lines
            if line.startswith(">"):
                name = line[1:].split()[0]  # first word after ">"
                sequences[name] = ""
            else:
                if name is None:
                    raise ValueError("File does not start with a '>' header line")
                sequences[name] += line.upper()
    return sequences


def gc_content(sequence):
    """Return the percentage of G and C bases."""
    if len(sequence) == 0:
        return 0.0
    gc = sequence.count("G") + sequence.count("C")
    return gc / len(sequence) * 100


def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    pairs = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(pairs.get(base, "N") for base in reversed(sequence))


def main():
    if len(sys.argv) != 2:
        print("Usage: python fasta_stats.py <file.fasta>")
        sys.exit(1)

    sequences = read_fasta(sys.argv[1])
    print(f"Sequences found: {len(sequences)}\n")

    for name, sequence in sequences.items():
        print(f"Name: {name}")
        print(f"  Length: {len(sequence)} bases")
        print(f"  GC content: {gc_content(sequence):.1f}%")
        counts = {base: sequence.count(base) for base in "ACGT"}
        print(f"  Base counts: {counts}")
        print(f"  Reverse complement (first 30 bases): {reverse_complement(sequence)[:30]}")
        print()


if __name__ == "__main__":
    main()
