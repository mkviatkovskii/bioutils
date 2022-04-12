#!/usr/bin/env python
"""
Find list of gene-like sequences of max common length for given
input list of sequences. Input is in fasta format, output is in
fasta-2 format.

Requires Python 3.6+ and BioPython 1.71+
"""
import argparse
import logging
import sys
import traceback
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Union

from Bio import AlignIO  # type: ignore
from Bio.Align import MultipleSeqAlignment  # type: ignore
from Bio.Seq import Seq  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore

logging.basicConfig(
    datefmt="%Y-%m-%d %H:%M:%S",
    format="[%(levelname)s] %(asctime)s %(message)s",
    level=logging.DEBUG,
)


def find_start_codons(string: Union[str, Seq]) -> List[int]:
    """
    Find all positions of start codon triplet in given string-like
    :param string: Input string-like object of triplets
    :return: List of start positions for start codon, empty if none found
    """
    result = []
    for i in range(0, len(string), 3):
        if string[i : i + 3] == "ATG":
            result.append(i)
    return result


def is_stop_codon(string: Union[str, Seq]) -> bool:
    """
    Check if input string-like object starts with stop codon
    Stop codons for DNA are: TAA, TAG, TGA
    :param string: String-like to check (could be str or BioPython Seq)
    :return: True if input string-like starts with stop codon
    """
    if string[0] == "T":
        if string[1] == "A":
            if string[2] in ("A", "G"):
                return True
        elif string[1] == "G" and string[2] == "A":
            return True
    return False


def find_stop_codon(string: Union[str, Seq]) -> Optional[int]:
    """
    Find stop codon position in given string-like object
    :param string: String-like to check (could be str or BioPython Seq)
    :return: If stop codon found, then return position of it's start. Return
             None if stop codon not found.
    """
    for i in range(0, len(string) - 2, 3):
        if is_stop_codon(string[i:]):
            return i
    return None


def expand_sequence(sequence: Seq) -> List[Seq]:
    """
    Clear sequence of gaps and find all possible variants of input sequence.
    There could be 6 options:
    * Forward, 0-, 1- and 2- shifted.
    * Reverse-Complement, 0-, 1-, and 2- shifted.
    :param sequence: input sequence, BioPython Seq.
    :return: List of 6 Seq objects
    """
    output = []
    seq: Seq = sequence.replace("-", "N")
    reverse_seq = seq.reverse_complement()
    for i in range(3):
        output.append(seq[i:].upper())
        output.append(reverse_seq[i:].upper())
    return output


def find_gene_like_in_sequence(input_sequence: Seq) -> Dict[int, Set[Seq]]:
    """
    Find all gene-like sequences in given sequence.
    Gene-like sequence starts with start codon triplet and stops with codon
    triplet. This method takes into account all possible sequences in given
    raw sequence (forward and reverse complement, with 1- and 2- shift).
    :param input_sequence: Input sequence, BioPython Seq object
    :return: Dictionary (sequence length -> set of gene-like sequences)
    """
    all_frames = set()
    result = defaultdict(set)
    for sequence in expand_sequence(input_sequence):
        for start_position in find_start_codons(sequence):
            stop_position = find_stop_codon(sequence[start_position:])
            if not stop_position:
                continue
            gene_sequence = sequence[
                start_position : start_position + stop_position + 3
            ]
            all_frames.add(str(gene_sequence))
            result[len(gene_sequence)].add(gene_sequence)
    return result


def sequences_similarity(seq1: Seq, seq2: Seq) -> float:
    """
    Return ratio of common symbols in two sequences,
    1.0 for identical sequences
    """
    common_symbols = 0
    # pylint: disable=consider-using-enumerate
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            common_symbols += 1
    return common_symbols / len(seq1)


def are_sequences_similar(
    sequences: List[Seq], min_similarity: float = 0.95
) -> bool:
    """Check if list of sequences are all similar enough"""
    seq0 = sequences[0]
    for seq in sequences[1:]:
        similarity = sequences_similarity(seq0, seq)
        if similarity < min_similarity:
            return False
    return True


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input-file", help="Input file in fasta format"
    )
    parser.add_argument(
        "-o", "--output-file", help="Output file in fasta2 format"
    )
    parser.add_argument(
        "-ms",
        "--min-similarity",
        default=0.95,
        type=float,
        help="Minimal similarity of sequences",
    )
    return parser.parse_args(sys.argv[1:])


def main(
    input_file: Path, output_file: Path, min_similarity: float = 0.5
) -> None:
    """
    Parse command-line arguments and run the search procedure.
    First argument: input Fasta file
    Second argument: output Fasta-2 file
    """

    logging.info(f"{input_file}: reading input data...")
    input_fasta_recs = AlignIO.read(input_file, "fasta")

    # Find all gene-like sequences grouped by their length for all records
    gene_like_dicts = [
        find_gene_like_in_sequence(rec.seq) for rec in input_fasta_recs
    ]

    # Find max length of sequence common for all given records
    seq_lengths = sorted(
        set.intersection(
            *[set(seq_lengths) for seq_lengths in gene_like_dicts]
        ),
        reverse=True,
    )

    # Find the longest sequence length for that all
    # sequences are similar enough
    max_seq_length = None
    records = []
    for seq_length in seq_lengths:
        # Prepare output records
        seqs = [
            tuple(gene_like_dict[seq_length])[0]
            for gene_like_dict in gene_like_dicts
        ]
        if are_sequences_similar(seqs, min_similarity):
            max_seq_length = seq_length
            for i, seq in enumerate(seqs):
                records.append(
                    SeqRecord(
                        seq,
                        id=input_fasta_recs[i].id,
                        name=input_fasta_recs[i].name,
                        description="",
                    )
                )
            break

    if not records:
        raise ValueError(
            f"Could not found similar sequences, "
            f"minimal similarity is set to {min_similarity}"
        )

    logging.info(
        f"{input_file}: found gene-like sequence of length {max_seq_length}"
    )

    # Write results
    result = MultipleSeqAlignment(records)
    AlignIO.write(result, output_file, "fasta-2line")

    logging.info(f"{input_file}: output saved to {output_file.resolve()}")


if __name__ == "__main__":
    args = parse_args()
    try:
        main(
            Path(args.input_file).resolve(),
            Path(args.output_file),
            args.min_similarity,
        )
    # pylint: disable=broad-except
    except Exception as e:
        logging.critical(
            f"Fatal error during processing {args.input_file}: {e}"
        )
        logging.debug(traceback.format_exc())
