#!/usr/bin/env python

import pysam
import argparse
import re


def parse_arguments():
    parser = argparse.ArgumentParser(description="Reindex a remapped BAM file without using the .fai index.")
    parser.add_argument("--bam", type=str, required=True, help="Path to the input BAM file.")
    parser.add_argument("--output", type=str, required=True, help="Path to the output BAM file.")
    parser.add_argument("--header", type=str, required=False, help="Path to the header file.")
    return parser.parse_args()


def reindex_bam(args):
    """
    Reindex a BAM file using a new header and remap the reference sequences.

    Args:
        args (Namespace): A namespace object containing the following attributes:
            - bam (str): Path to the input BAM file.
            - header (str): Path to the BAM file containing the new header.
            - output (str): Path to the output BAM file.

    This function performs the following steps:
        1. Opens the input BAM file and the header BAM file.
        2. Extracts the header from the header BAM file.
        3. Creates an output BAM file using the extracted header.
        4. Iterates through each read in the input BAM file and remaps the reference sequences.
        5. Writes the remapped reads to the output BAM file.

    Note:
        - Unmapped reads are skipped.
        - The reference sequence information is extracted from the read's reference name.
        - The reference start position is adjusted based on the extracted reference sequence information.
        - If the read is paired, the mate reference information is also remapped.
    """
    # Open the input BAM and the header BAM
    tmp_bam = pysam.AlignmentFile(args.bam, "rb")
    header_bam = pysam.AlignmentFile(args.header, "rb")
    new_header = header_bam.header
    header_bam.close()

    # Create output BAM using the same header
    with pysam.AlignmentFile(args.output, "wb", header=new_header) as outf:
        # Copy reads from input to output while printing the reference sequence for each read
        dict_ref_seq = {}
        for read in tmp_bam:
            if read.is_unmapped:
                continue
            ref_info = read.reference_name.split("|")[0]
            name, pos = ref_info.split(":")
            start, end = map(int, pos.split("-"))
            new_read = pysam.AlignedSegment()
            new_read.query_name = read.query_name
            new_read.query_sequence = read.query_sequence
            new_read.flag = read.flag
            ref_id = outf.get_tid(str(name))
            new_read.reference_id = ref_id

            new_read.reference_start = read.reference_start + start - 1
            new_read.mapping_quality = read.mapping_quality
            new_read.cigar = read.cigar
            if read.is_paired:
                mate_ref_info = read.next_reference_name.split("|")[0]  # e.g., "chr1:15001-25000"
                mate_name, mate_pos = mate_ref_info.split(":")  # mate_name = "chr1", mate_pos = "15001-25000"
                mate_start, mate_end = map(int, mate_pos.split("-"))
                new_read.next_reference_id = mate_name
                new_read.next_reference_start = start + mate_start - 1
            new_read.template_length = read.template_length
            new_read.query_qualities = read.query_qualities
            new_read.tags = read.tags

            outf.write(new_read)
        tmp_bam.close()


def main():
    args = parse_arguments()
    reindex_bam(args)


if __name__ == "__main__":
    main()