#!/usr/bin/env python3

from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import re
import argparse
import sys
import os
from Bio.Seq import reverse_complement
import pandas as pd
import numpy as np
from Bio.Seq import Seq

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def parse_fasta_file(file_path):
    """
    Parses a FASTA file and organizes sequences by their key.
    """
    sequences = []
    with open(file_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append(record)
    return sequences


def perform_alignment(
    gene_seq,
    pseudo_seq,
    gene_strandedness,
    pseudo_strandedness,
    match_score,
    mismatch_score,
    gap_open,
    gap_extend,
    use_global=False,
):
    """
    Performs local pairwise alignment between two sequences, considering strandedness.
    """
    if pseudo_strandedness != gene_strandedness:
        pseudo_seq = reverse_complement(pseudo_seq)
    if use_global:
        print("Using global alignment")
        alignments = pairwise2.align.globalxx(gene_seq, pseudo_seq)
    else:
        alignments = pairwise2.align.localms(gene_seq, pseudo_seq, match_score, mismatch_score, gap_open, gap_extend)
    return alignments[0]


def parse_coordinates(coord_str):
    """
    Parses a coordinate string in the format 'chr:start-end' and returns chromosome, start, and end as integers.
    """
    chrom, positions = coord_str.split(":")
    start, end = positions.split("-")[0], positions.split("-")[1]
    return chrom, int(start), int(end)


def extract_differences(
    aligned_gene_seq,
    aligned_pseudo_seq,
    output_file,
    gene_chrom,
    pseudo_chrom,
    gene_start,
    gene_end,
    pseudo_start,
    pseudo_end,
):
    """
    Extracts differences between aligned gene and pseudogene sequences and writes them to a TSV file.
    """
    gene_pos = int(gene_start)
    pseudo_pos = int(pseudo_start)

    differences = []

    for i in range(len(aligned_gene_seq)):
        gene_base = aligned_gene_seq[i].upper()
        pseudo_base = aligned_pseudo_seq[i].upper()

        if gene_base != "-" and pseudo_base != "-":
            if gene_base != pseudo_base:
                # Mismatch
                differences.append([gene_chrom, gene_pos, pseudo_chrom, pseudo_pos, "X", gene_base, pseudo_base])
            else:
                # Match
                differences.append([gene_chrom, gene_pos, pseudo_chrom, pseudo_pos, "M", gene_base, pseudo_base])
            gene_pos += 1
            pseudo_pos += 1
        elif gene_base != "-" and pseudo_base == "-":
            # Deletion in pseudogene
            differences.append([gene_chrom, gene_pos, pseudo_chrom, "-", "D", gene_base, "-"])
            gene_pos += 1
        elif gene_base == "-" and pseudo_base != "-":
            # Insertion in pseudogene
            differences.append([gene_chrom, "-", pseudo_chrom, pseudo_pos, "I", "-", pseudo_base])
            pseudo_pos += 1
        else:
            # Both gaps (unlikely in valid alignments)
            continue

    df = pd.DataFrame(
        differences,
        columns=[
            "Gene_Chromosome",
            "Gene_Position",
            "Pseudogene_Chromosome",
            "Pseudogene_Position",
            "Difference_Type",
            "Gene_Base",
            "Pseudogene_Base",
        ],
    )

    df["Gene_Position"] = df["Gene_Position"].replace("-", np.nan)
    df["Pseudogene_Position"] = df["Pseudogene_Position"].replace("-", np.nan)
    min_gene_position_index = df["Gene_Position"].idxmin()
    max_gene_position_index = df["Gene_Position"].idxmax()
    min_pseudogene_position_index = df["Pseudogene_Position"].idxmin()
    max_pseudogene_position_index = df["Pseudogene_Position"].idxmax()
    length_gene = max_gene_position_index - min_gene_position_index
    length_pseudogene = max_pseudogene_position_index - min_pseudogene_position_index

    df["Gene_Position"] = df["Gene_Position"].fillna(-1).astype(int)
    df["Pseudogene_Position"] = df["Pseudogene_Position"].fillna(-1).astype(int)

    if length_gene < length_pseudogene:
        df = df.loc[min_gene_position_index:max_gene_position_index]
    else:
        df = df.loc[min_pseudogene_position_index:max_pseudogene_position_index]

    # Remove trailing insertion or deletion pairs
    while not df.empty and df.iloc[-1]["Difference_Type"] in ["I", "D"]:
        df = df.iloc[:-1]
    while not df.empty and df.iloc[0]["Difference_Type"] in ["I", "D"]:
        df = df.iloc[1:]
    df.to_csv(output_file, sep="\t", index=False)

    print(f"Differences have been written to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Process sequences.")
    parser.add_argument("--gene", required=True, help="Path to the gene reference FASTA file")
    parser.add_argument("--pseudogene", required=True, help="Path to the pseudogene reference FASTA file")
    parser.add_argument("--output", required=True, help="Path to the output of TSV file")
    parser.add_argument(
        "--strandedness",
        "-s",
        required=False,
        choices=["positive", "negative"],
        help="Strandedness of the sequences (if not specified in sequence headers).",
    )
    parser.add_argument(
        "--use_global",
        "-g",
        action="store_true",
        help="Use global alignment instead of local alignment.",
    )
    args = parser.parse_args()

    gene_reference = args.gene
    pseudogene_reference = args.pseudogene
    output_file = args.output

    # Parse sequences from FASTA files
    gene_sequences = parse_fasta_file(gene_reference)
    pseudogene_sequences = parse_fasta_file(pseudogene_reference)

    # Alignment scoring parameters
    alignment_params = {"match_score": 2, "mismatch_score": -3, "gap_open": -5, "gap_extend": -2}

    # Process each gene sequence individually
    for gene_number, gene_seq_record in enumerate(gene_sequences):
        for pseudogene_number, pseudogene_seq_record in enumerate(pseudogene_sequences):

            gene_sq_coordinates = re.search(r"[0-9XY]+:\d+-\d+", gene_seq_record.description).group()
            pseudogene_sq_coordinates = re.search(r"[0-9XY]+:\d+-\d+", pseudogene_seq_record.description).group()

            gene_chrom, gene_start, gene_end = parse_coordinates(gene_sq_coordinates)
            pseudo_chrom, pseudo_start, pseudo_end = parse_coordinates(pseudogene_sq_coordinates)

            # Extract strandedness from sequence headers or use provided strandedness
            gene_strandedness = args.strandedness
            pseudo_strandedness = args.strandedness

            gene_header_strand = re.search(r"strand:(\+|\-)", gene_seq_record.description)
            if gene_header_strand:
                gene_strandedness = 'positive' if gene_header_strand.group(1) == '+' else 'negative'

            pseudo_header_strand = re.search(r"strand:(\+|\-)", pseudogene_seq_record.description)
            if pseudo_header_strand:
                pseudo_strandedness = 'positive' if pseudo_header_strand.group(1) == '+' else 'negative'

            # Default to '+' if strandedness is not specified
            if not gene_header_strand:
                gene_strandedness = "negative"
            if not pseudo_header_strand:
                pseudo_strandedness = "negative"

            # Perform alignment considering strandedness
            best_alignment = perform_alignment(
                str(gene_seq_record.seq),
                str(pseudogene_seq_record.seq),
                gene_strandedness,
                pseudo_strandedness,
                **alignment_params,
                use_global=args.use_global,
            )
            aligned_exon_seq, aligned_pseudo_seq, score, begin, end = best_alignment

            # Extract and store differences

            extract_differences(
                aligned_gene_seq=aligned_exon_seq,
                aligned_pseudo_seq=aligned_pseudo_seq,
                output_file=output_file,
                gene_chrom=gene_chrom,
                gene_start=gene_start,
                pseudo_chrom=pseudo_chrom,
                pseudo_start=pseudo_start,
                gene_end=gene_end,
                pseudo_end=pseudo_end,
            )


if __name__ == "__main__":
    main()