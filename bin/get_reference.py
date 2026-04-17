#!/usr/bin/env python


import pandas as pd
from typing import Optional
from Bio import SeqIO
import argparse
import os
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
    chrom_query = str(chromosome)
    chrom_query_with_prefix = chrom_query if chrom_query.startswith("chr") else "chr" + chrom_query
    chrom_query_without_prefix = chrom_query.lstrip("chr") if chrom_query.startswith("chr") else chrom_query

    with open(reference) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            record_id = record.id
            if record_id in (chrom_query, chrom_query_with_prefix, chrom_query_without_prefix):
                start = int(start) - 1  # Convert to 0-based
                return str(record.seq[start:int(end)])

    raise ValueError(
        f"Chromosome '{chromosome}' not found in reference. "
        f"Tried: '{chrom_query}', '{chrom_query_with_prefix}', '{chrom_query_without_prefix}'"
    )


def validate_region(region):
    """Ensure chromosome has 'chr' prefix."""
    if re.match(r"^\d+$", region):
        return f"chr{region}"
    return region


def save_sequences_as_fasta(sequences_df, output_fasta):
    with open(output_fasta, "w") as fasta_file:
        for _, row in sequences_df.iterrows():
            coordinates = f"{row['chromosome']}:{row['start']}-{row['end']}"
            header_parts = [coordinates]
            if "gene_name" in row and pd.notna(row["gene_name"]):
                header_parts.append(str(row["gene_name"]))
            header_parts.extend(
                str(row[col])
                for col in row.index
                if col not in ["sequence", "chromosome", "start", "end", "gene_name", "strand", "score"]
                and pd.notna(row[col])
            )
            if "strand" in row and pd.notna(row["strand"]):
                header_parts.append(f"strand:{row['strand']}")
            header = "|".join(part for part in header_parts if part).strip("|")
            fasta_file.write(f">{header}\n")
            fasta_file.write(f"{row['sequence']}\n")
    print(f"Sequences saved to {output_fasta}")


def read_bed(bed_file, strandedness="-"):
    """
    Read a BED file and return it as a pandas DataFrame.

    Args:
      bed_file (str): Path to the BED file.
      strandedness (str): Default strand if not present in BED file.
    Returns:
      pd.DataFrame: DataFrame containing the BED file data.
    """
    bed = pd.read_csv(bed_file, sep=r"\s+", header=None, comment="#")
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
        strand = args.strand if args.strand else "-"
        bed_gene_df = read_bed(args.bed, strand)
        bed_gene_df["sequence"] = bed_gene_df.apply(
            lambda row: get_sequence(args.reference, row["chromosome"], row["start"], row["end"]),
            axis=1,
        )
        save_sequences_as_fasta(bed_gene_df, args.output)

    elif args.region:
        region = re.split(r"[:\-]", args.region)
        chrom = validate_region(region[0])
        sequence = get_sequence(
            args.reference,
            chrom,
            int(region[1]),
            int(region[2]),
        )
        strand = args.strand if args.strand else "-"
        save_sequences_as_fasta(
            pd.DataFrame(
                {
                    "chromosome": [chrom],
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
    parser.add_argument("--output", required=True, help="Path to the output FASTA file for gene sequences.")
    parser.add_argument("--strand", help="Known strand of the input sequence")

    group.add_argument("--region", help="Genomic region in the form chr:start-end")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()