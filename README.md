# Longest gene search tool

## Summary
Tool for finding list of gene-like sequences of max common length for
given input list of sequences (aligned contigs). Input is in FASTA
format, output is in FASTA-2 format.

## Requirements
* Python 3.6+
* BioPython 1.71+

## Parameters
* `-i`, `--input-file`: input FASTA file with aligned contigs.
* `-o`, `--output-file`: path to output FASTA-2 file.
* `-ms`, `--min-similarity`: Minimal allowed similarity of sequences,
in [0, 1] interval. By default is 0.95.

## Example
```bash
python3 ./find_gene_like_in_contigs.py -i examples/gregalis.fa -o output.fa
```
