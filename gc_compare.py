"""Compare the length and GC content of the two genomes.

Needs the two FASTA files in the data folder:
    data/sars_cov_2.fasta
    data/ratg13.fasta

Usage:
    python gc_compare.py
"""
from pathlib import Path

from fasta_stats import read_fasta, gc_content

FILES = {
    "SARS-CoV-2 (NC_045512.2)": Path("data/sars_cov_2.fasta"),
    "RaTG13 (MN996532.2)": Path("data/ratg13.fasta"),
}


def main():
    for label, path in FILES.items():
        if not path.exists():
            print(f"Missing file: {path}")
            return
        sequences = read_fasta(path)
        sequence = next(iter(sequences.values()))
        print(f"{label}")
        print(f"  Length: {len(sequence)} bases")
        print(f"  GC content: {gc_content(sequence):.2f}%\n")


if __name__ == "__main__":
    main()
