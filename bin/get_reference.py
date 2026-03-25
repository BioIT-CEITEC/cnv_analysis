#!/usr/bin/env python


import pandas as pd
from typing import Optional
from Bio import SeqIO
import argparse
import os
from bed_reader import open_bed, sample_file
import re


def get_sequence(reference, chromosome, start, end):
    """
    Get the sequence from a reference genome file.

    Args:
        reference (str): Path to the reference genome file.
        chromosome (str): Chromosome name.
        start (int): Start position of the sequence (1-based).
        end (int): End position of the sequence (1-based).

    Returns:
        str: The sequence from the reference genome.
    """
    with open(reference) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if not chromosome.startswith("chr"):
                chromosome = "chr" + chromosome
            if record.id == chromosome:
                start = int(start) - 1
                return str(record.seq[start:end])


def validate_region(region):
    if re.match(r"^\d+$", region):
        return f"chr{region}"
    return region


def save_sequences_as_fasta(sequences_df, output_fasta):
    with open(output_fasta, "w") as fasta_file:
        for _, row in sequences_df.iterrows():
            coordinates = f"{row['chromosome']}:{row['start']}-{row['end']}"
            header_parts = [coordinates]
            if "gene_name" in row and row["gene_name"]:
                header_parts.append(row["gene_name"])
            # Add other columns excluding specified ones
            header_parts.extend(
                str(row[col])
                for col in row.index
                if col not in ["sequence", "chromosome", "start", "end", "gene_name", "strand", "score"]
            )
            if "strand" in row:
                header_parts.append(f"strand:{row['strand']}")
            # Construct header
            header = "|".join(header_parts).strip("|")
            fasta_file.write(f">{header}\n")
            fasta_file.write(f"{row['sequence']}\n")
    print(f"Sequences saved to {output_fasta}")


def read_bed(bed_file, strandedness="-"):
    """
    Read a BED file and return it as a pandas DataFrame.

    Args:
      bed_file (str): Path to the BED file.
    Returns:
      pd.DataFrame: DataFrame containing the BED file data.
    """
    bed = pd.read_csv(bed_file, sep="\s+", header=None, comment="#")
    BED_COLUMNS = [
        "chromosome",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "thick_start",
        "thick_end",
        "item_rgb",
        "block_count",
        "block_sizes",
        "block_starts",
    ]
    n_columns = bed.shape[1]
    bed.columns = BED_COLUMNS[:n_columns]
    if "strand" not in bed.columns:
        bed["strand"] = strandedness
    return bed


def main():
    args = parse_arguments()

    if args.bed:
        if args.strand:
            bed_gene_df = read_bed(args.bed, args.strand)
        bed_gene_df = read_bed(args.bed)
        bed_gene_df["sequence"] = bed_gene_df.apply(
            lambda row: get_sequence(args.reference, row["chromosome"], row["start"], row["end"]),
            axis=1,
        )
        save_sequences_as_fasta(bed_gene_df, args.output)

    elif args.region:
        region = re.split("[:\-]", args.region)
        sequence = get_sequence(
            args.reference,
            region[0],
            int(region[1]),
            int(region[2]),
        )
        strand = args.strand if args.strand else "-"
        save_sequences_as_fasta(
            pd.DataFrame(
                {
                    "chromosome": [region[0]],
                    "start": [int(region[1])],
                    "end": [int(region[2])],
                    "sequence": [sequence],
                    "strand": [strand],
                }
            ),
            args.output,
        )


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Extract sequences from GTF (with name and attributes) or BED files and save as FASTA."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument("--reference", required=True, help="Path to the reference genome FASTA file.")

    group.add_argument("--bed", help="Path to the BED file containing gene regions.")
    parser.add_argument("--output", help="Path to the output FASTA file for gene sequences.")
    parser.add_argument("--strand", help="Known strand of the input sequence")

    group.add_argument("--region", help="One life of BED file with form chr start end")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()