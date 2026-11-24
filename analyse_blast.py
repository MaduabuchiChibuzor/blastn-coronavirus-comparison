"""Summarise a pairwise BLASTn result and show where two genomes differ.

Reads the XML file saved from NCBI BLAST ("Align two or more sequences"),
prints a summary, and saves a table and a plot in the results folder.

Usage:
    python analyse_blast.py                    # uses the first .xml file in data/
    python analyse_blast.py data/some_file.xml
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # lets the script save plots without opening a window
import matplotlib.pyplot as plt
from Bio.Blast import NCBIXML

DATA_FOLDER = Path("data")
RESULTS_FOLDER = Path("results")
WINDOW = 2000 # bases of the query genome per window

# Gene positions on the query genome (SARS-CoV-2, NC_045512.2).
# Check these against the annotation on the NCBI record before you rely on them.
GENES = {
    "ORF1ab": (266, 21555),
    "S": (21563, 25384),
    "ORF3a": (25393, 26220),
    "E": (26245, 26472),
    "M": (26523, 27191),
    "ORF6": (27202, 27387),
    "ORF7a": (27394, 27759),
    "ORF7b": (27756, 27887),
    "ORF8": (27894, 28259),
    "N": (28274, 29533),
    "ORF10": (29558, 29674),
}


def read_blast(xml_path):
    """Read the first alignment from a BLAST XML file."""
    with open(xml_path) as file:
        record = NCBIXML.read(file)
    if not record.alignments:
        raise ValueError("No alignments found in this BLAST file")
    alignment = record.alignments[0]
    hsp = alignment.hsps[0]
    return record, alignment, hsp


def summarise(record, alignment, hsp):
    """Work out the main numbers from one alignment."""
    mismatches = hsp.align_length - hsp.identities - hsp.gaps
    return {
        "query": record.query.split()[0],
        "subject": alignment.hit_def.split()[0],
        "query_length": record.query_length,
        "subject_length": alignment.length,
        "alignment_length": hsp.align_length,
        "identities": hsp.identities,
        "mismatches": mismatches,
        "gap_positions": hsp.gaps,
        "percent_identity": round(hsp.identities / hsp.align_length * 100, 2),
        "query_coverage": round((hsp.query_end - hsp.query_start + 1) / record.query_length * 100, 2),
        "subject_coverage": round((hsp.sbjct_end - hsp.sbjct_start + 1) / alignment.length * 100, 2),
        "bit_score": hsp.bits,
        "e_value": hsp.expect,
    }


def column_matches(hsp):
    """For every alignment column, give the query position and whether the bases match."""
    columns = []
    query_position = hsp.query_start - 1
    for query_base, subject_base in zip(hsp.query, hsp.sbjct):
        if query_base != "-":
            query_position += 1
        columns.append((query_position, query_base == subject_base))
    return columns


def identity_by_window(columns, window=WINDOW):
    """Percent identity for each block of the query genome."""
    totals = {}
    for position, is_match in columns:
        block = (position - 1) // window
        matches, count = totals.get(block, (0, 0))
        totals[block] = (matches + is_match, count + 1)
    rows = []
    for block in sorted(totals):
        matches, count = totals[block]
        rows.append({
            "start": block * window + 1,
            "end": (block + 1) * window,
            "percent_identity": round(matches / count * 100, 2),
        })
    return rows


def identity_by_gene(columns, genes=GENES):
    """Percent identity inside each gene, counting query positions only."""
    rows = []
    for name, (start, end) in genes.items():
        inside = [m for position, m in columns if start <= position <= end]
        if inside:
            rows.append({
                "gene": name,
                "start": start,
                "end": end,
                "percent_identity": round(sum(inside) / len(inside) * 100, 2),
            })
    return rows


def save_csv(rows, path):
    with open(path, "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def plot_identity(window_rows, path, genes=GENES):
    middles = [(r["start"] + r["end"]) / 2 for r in window_rows]
    identities = [r["percent_identity"] for r in window_rows]

    fig, (ax, track) = plt.subplots(
        2, 1, figsize=(11, 5.5), sharex=True,
        gridspec_kw={"height_ratios": [4, 1], "hspace": 0.05},
    )
    ax.plot(middles, identities, color="black", linewidth=1.2)
    ax.set_ylabel(f"Identity per {WINDOW} base window (%)")
    ax.set_title("Identity between SARS-CoV-2 and RaTG13 along the genome")
    ax.set_ylim(min(identities) - 2, 100.5)

    # Gene track: alternate two greys. Only the longer genes get a label.
    shades = ["#bdbdbd", "#7a7a7a"]
    for i, (name, (start, end)) in enumerate(genes.items()):
        track.broken_barh([(start, end - start)], (0.2, 0.5), facecolors=shades[i % 2])
        if end - start >= 600:  # short genes are left unlabelled, see identity_by_gene.csv
            track.text((start + end) / 2, 0.85, name, ha="center", va="bottom", fontsize=9)
    track.set_ylim(0, 1.4)
    track.set_yticks([])
    track.set_xlabel("Position on SARS-CoV-2 genome (NC_045512.2)")
    for side in ("left", "right", "top"):
        track.spines[side].set_visible(False)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    if len(sys.argv) > 1:
        xml_path = Path(sys.argv[1])
    else:
        xml_files = sorted(DATA_FOLDER.glob("*.xml"))
        if not xml_files:
            print("No .xml file found in the data folder.")
            sys.exit(1)
        xml_path = xml_files[0]

    RESULTS_FOLDER.mkdir(exist_ok=True)
    record, alignment, hsp = read_blast(xml_path)
    summary = summarise(record, alignment, hsp)

    print(f"File: {xml_path.name}\n")
    for key, value in summary.items():
        print(f"{key}: {value}")

    columns = column_matches(hsp)
    windows = identity_by_window(columns)
    genes = identity_by_gene(columns)

    save_csv([summary], RESULTS_FOLDER / "summary.csv")
    save_csv(windows, RESULTS_FOLDER / "identity_by_window.csv")
    save_csv(genes, RESULTS_FOLDER / "identity_by_gene.csv")
    plot_identity(windows, RESULTS_FOLDER / "identity_plot.png")

    print("\nIdentity by gene:")
    for row in genes:
        print(f"  {row['gene']:7s} {row['percent_identity']:.2f}%")
    lowest = min(windows, key=lambda r: r["percent_identity"])
    print(f"\nLowest identity window: {lowest['start']} to {lowest['end']} "
          f"({lowest['percent_identity']}%)")
    print(f"\nSaved tables and plot in {RESULTS_FOLDER}/")


if __name__ == "__main__":
    main()
