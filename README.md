# BLASTn comparison of SARS-CoV-2 and RaTG13

A small project that compares two coronavirus genomes using NCBI BLASTn, then uses Python to summarise the result and show which parts of the genome differ most.

## Why I made it

I first did this comparison as coursework on my MSc Bioinformatics course. I re-ran it and saved the outputs so the work can be repeated, and I wrote a Python script with help from an AI assistant, to read the BLAST results and analyse them. I ran the tests and changed the settings to see how it works.

## Data

|Role|Accession|Description|Length|
|-|-|-|-|
|Query|NC\_045512.2|SARS-CoV-2 isolate Wuhan-Hu-1, complete genome|29,903 bases|
|Subject|MN996532.2|Bat coronavirus RaTG13, complete genome|29,855 bases|

Both sequences come from NCBI.

## Method

1. Ran nucleotide BLAST on the NCBI website using "Align two or more sequences" on 24 November 2026 (job ID GB09Z26C114).
2. Program: BLASTN 2.17.0+ with the default settings (match 1, mismatch -2, linear gap costs, expect 0.05, low complexity filter).
3. Saved the results as XML, text and a hit table in the `data` folder.
4. `analyse\_blast.py` reads the XML with Biopython, works out the summary, and calculates identity in windows of 500 bases and inside each gene.

## Results

|Measure|Value|
|-|-|
|Aligned length|29,877 columns|
|Identical bases|28,723|
|Percent identity|96.14%|
|Mismatches|1,130|
|Gap positions|24|
|Query coverage|99.91%|
|Subject coverage|100%|
|E value|0.0|

These numbers match the NCBI hit table.

!\[Identity along the genome](results/identity\_plot.png)

Identity is not even across the genome. The spike gene (S) has the lowest identity at 92.89%, and the lowest 500 base window is at positions 22,501 to 23,000 (87.8%). The envelope gene (E) has the highest at 99.56%. The full table is in `results/identity\_by\_gene.csv`. Short genes are not labelled on the plot, so check the table for ORF6, ORF7a, ORF7b, ORF8 and ORF10.

## GC content

SARS-CoV-2 (NC\_045512.2): 29,903 bases, GC content 37.97%.

RaTG13 (MN996532.2): 29,855 bases, GC content 38.03%.

The two genomes have almost the same GC content.

## Limits of this analysis

* This is one pairwise BLAST alignment. It shows where the two sequences differ, but not why, and it does not show evolutionary history.
* The 500 base window size is a choice. A different size gives a slightly different plot.
* Gene positions in `analyse\_blast.py` are taken from the SARS-CoV-2 annotation and should be checked against the NCBI record.

## How to run it

You need Python 3.

```
pip install -r requirements.txt
python analyse\_blast.py
python gc\_compare.py
python -m pytest
```

`gc\_compare.py` needs the two FASTA files saved as `data/sars\_cov\_2.fasta` and `data/ratg13.fasta`.

## What I learned

I learned how to read a BLAST result and what percent identity, mismatches and query coverage mean. I also saw that two genomes can be 96% similar overall but differ more in some regions, like the spike gene. When I changed the window size from 500 to 2000 bases, the dip in the plot became shallower and moved, because each point then averaged more bases. Next I would like to add more genomes to the comparison and write more of the code myself.

## Author

Chibuzor Daniels Maduabuchi, MSc Bioinformatics, Teesside University

