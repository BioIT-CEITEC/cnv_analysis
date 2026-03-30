#!/usr/bin/env python


import pysam
import argparse
import pandas as pd
from Bio import pairwise2
from Bio.Seq import reverse_complement
from collections import defaultdict
import sys
from tqdm import tqdm


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process BAM file and realign reads to pseudogene if necessary.")
    parser.add_argument("--tsv", required=True, help="Path to the ref_diffs TSV file containing reference differences.")
    parser.add_argument("--bam", required=True, help="Path to the BAM file with read alignments.")
    parser.add_argument("--output", default="detailed_read_table.tsv", help="Path for the output detailed TSV file.")
    parser.add_argument("--min_mapq", type=int, default=20, help="Minimum mapping quality to consider a read.")
    parser.add_argument("--flank_size", type=int, default=3, help="Flanking region size for realignment.")
    return parser.parse_args()


def load_ref_diffs(ref_diffs_file):
    # Load the ref_diffs file into a DataFrame
    df = pd.read_csv(ref_diffs_file, sep="\t")
    return df


def clean_chromosome_name(chromosome):
    # Clean the chromosome name
    if isinstance(chromosome, str) and chromosome.startswith("chr"):
        return chromosome[3:].strip()
    return chromosome


def get_insertion_intervals(ref_diffs_df):
    # Filter insertions into the pseudogene
    insertions = ref_diffs_df[ref_diffs_df["Difference_Type"] == "I"].copy()
    insertions["Index"] = insertions.index

    # Initialize lists to store results
    insertion_intervals = []

    # Loop over each insertion
    for idx, row in insertions.iterrows():
        insertion_pos = row["Pseudogene_Position"]

        # Find closest smaller gene position
        smaller_positions = ref_diffs_df.loc[: idx - 1]
        smaller_positions = smaller_positions[smaller_positions["Gene_Position"] != -1]
        if not smaller_positions.empty:
            closest_smaller = smaller_positions.iloc[-1]["Gene_Position"]
        else:
            closest_smaller = None

        # Find closest greater gene position
        greater_positions = ref_diffs_df.loc[idx + 1 :]
        greater_positions = greater_positions[greater_positions["Gene_Position"] != -1]
        if not greater_positions.empty:
            closest_greater = greater_positions.iloc[0]["Gene_Position"]
        else:
            closest_greater = None

        # Store the interval
        insertion_intervals.append(
            (
                idx,
                insertion_pos,
                closest_smaller,
                closest_greater,
            )
        )

    return insertion_intervals


def realign_read(read_seq_for_alignment, reference):
    """
    Realigns a read sequence to a reference sequence using local alignment.

    This function performs a local alignment between the provided read sequence and
    the reference sequence using specified scoring parameters. It then generates a CIGAR
    string from the best alignment and compares it with a provided CIGAR string.

    Parameters:
        read_seq_for_alignment (str): The sequence of the read to be aligned.
        reference (str): The reference sequence to align against.
        cigar_1 (str): The CIGAR string to compare with the generated alignment.

    Returns:
        tuple:
            bool: True if the CIGAR strings match, False otherwise.
            tuple or None: The best alignment if found, else None.
            str or None: The generated CIGAR string if alignment is found, else None.
    """
    # Perform local alignment
    best_alignment, cigar_2 = None, None
    alignments = pairwise2.align.localms(read_seq_for_alignment, reference, 2, -2, -5, -2)
    if alignments:
        best_alignment = alignments[0]
        cigar_2 = generate_cigar_from_alignment(best_alignment)
    return best_alignment, cigar_2


def process_reads(bam_file, ref_diffs_df, output, min_mapq, flank_size):
    bam = pysam.AlignmentFile(bam_file, "rb")
    pseudogene_insertion_intervals = get_insertion_intervals(ref_diffs_df)
    # Extract insertion regions from ref_diffs
    read_sequences = {}
    read_positions = {}
    alignment_results = {}
    all_reads = []
    start = ref_diffs_df[ref_diffs_df["Gene_Position"] != -1]["Gene_Position"].min()
    end = ref_diffs_df[ref_diffs_df["Gene_Position"] != -1]["Gene_Position"].max()
    chromosome = clean_chromosome_name(ref_diffs_df["Gene_Chromosome"].iloc[0])
    start_idx = ref_diffs_df[ref_diffs_df["Gene_Position"] == start].index[0]
    end_idx = ref_diffs_df[ref_diffs_df["Gene_Position"] == end].index[0]
    ref_diffs_df = ref_diffs_df.iloc[start_idx:end_idx + 1]
    print(start, end)

    # Iterate over each read in the specified region
    for read in bam.fetch(str(chromosome), start=start, end=end):
        if read.is_unmapped or read.mapping_quality < min_mapq or read.is_duplicate or read.is_supplementary:
            continue
        if read.is_reverse:
            read.query_sequence = reverse_complement(read.query_sequence)
            read.query_name = read.query_name + ".1"
        read_start = read.reference_start
        read_end = read.reference_end + 1
        read_name = read.query_name
        read_seq = read.query_sequence
        all_reads.append(read_name)

        # Retrieve the CIGAR string for the read
        cigar_tuples = read.cigartuples  # List of (operation, length) tuples
        cigar_1 = summarize_cigar(cigar_tuples)

        # Check if the read overlaps with an insertion region
        overlaps_insertion = False
        sum_insertions = 0
        # start, end = last non-null gene positions
        for _, _, ins_start, ins_end in pseudogene_insertion_intervals:
            # Adjusted logic to check for overlap with insertion positions
            if read_start <= ins_start and ins_end <= read_end:
                overlaps_insertion = True
                # insert_size is the length of the insertion
                insertion_start_idx = ref_diffs_df.index[ref_diffs_df["Gene_Position"] == ins_start]
                insertion_end_idx = ref_diffs_df.index[ref_diffs_df["Gene_Position"] == ins_end]
                if not insertion_start_idx.empty and not insertion_end_idx.empty:
                    insert_size = insertion_end_idx[0] - insertion_start_idx[0] - 1
                    sum_insertions += insert_size
                break

        has_insertion = any(op == 1 for op, _ in cigar_tuples)

        # Initialize variables for realignment
        realign_to_pseudogene = False
        cigar_2 = None

        if overlaps_insertion or has_insertion:
            read_start_idx = ref_diffs_df[ref_diffs_df["Gene_Position"] == read_start].index
            read_end_idx = ref_diffs_df[ref_diffs_df["Gene_Position"] == read_end].index
            if not read_start_idx.empty and not read_end_idx.empty:
                read_start_idx = read_start_idx[0]
                read_end_idx = read_end_idx[0]
                pseudogene_seq, new_reference_start_idx = get_sequence_reference(
                    ref_diffs_df, read_start_idx, read_end_idx, flank_size, sum_insertions, sequence_type="pseudogene"
                )
                read_seq_for_alignment = read_seq

                # Use the standalone realignment function
                best_alignment, cigar_2 = realign_read(read_seq_for_alignment, pseudogene_seq)
                if best_alignment is None:
                    realign_to_pseudogene = False
                elif compare_cigars(cigar_1, cigar_2) >= -1:
                    realign_to_pseudogene = True
            else:
                # Unable to map read positions to ref_diffs indices
                realign_to_pseudogene = False
        else:
            realign_to_pseudogene = False

        # Handle high CIGAR scores (here, we define a threshold arbitrarily)
        cigar_score_threshold = len(read_seq) // 2
        cigar_score = calculate_cigar_score(cigar_tuples)
        if cigar_score < cigar_score_threshold and not realign_to_pseudogene:
            realign_to_pseudogene = False

        # Store results
        if realign_to_pseudogene:
            # Update read sequence and positions
            read_sequences[read_name] = best_alignment.seqA  # Realigned read sequence
            read_positions[read_name] = generate_reference_positions(
                best_alignment, ref_diffs_df.iloc[new_reference_start_idx]["Gene_Position"]
            )
            alignment_results[read_name] = {
                "Original_CIGAR": cigar_1,
                "Realigned_CIGAR": cigar_2,
                "Alignment_Score": best_alignment.score,
            }
        else:
            # Keep original alignment
            read_sequences[read_name] = read_seq
            # Use the original reference positions
            read_positions[read_name] = read.get_reference_positions()
            # Convert positions to 1-based
            read_positions[read_name] = [pos + 1 if pos is not None else None for pos in read_positions[read_name]]
            alignment_results[read_name] = {
                "Original_CIGAR": cigar_1,
                "Alignment_Score": None,
            }

    bam.close()

    # Output the results to the detailed file
    create_detailed_table(ref_diffs_df, read_sequences, read_positions, alignment_results, all_reads, output)


def summarize_cigar(cigar_tuples):
    # Summarize the CIGAR string into counts of operations
    summary = {"M": 0, "I": 0, "D": 0, "X": 0, "S": 0}
    for op, length in cigar_tuples:
        if op == 0:
            summary["M"] += length
        elif op == 1:
            summary["I"] += length
        elif op == 2:
            summary["D"] += length
        elif op == 8:
            summary["X"] += length
        elif op == 4:
            summary["S"] += length
    return summary


def get_sequence_reference(ref_diffs_df, start_idx, end_idx, flank_size, insert_size, sequence_type="pseudogene"):
    # Construct the sequence based on ref_diffs
    if sequence_type not in ["gene", "pseudogene"]:
        raise ValueError("sequence_type must be 'gene' or 'pseudogene'")

    start_idx = max(0, start_idx - flank_size - insert_size)
    end_idx = min(len(ref_diffs_df) - 1, end_idx + flank_size + insert_size)
    ref_diffs_subset = ref_diffs_df.iloc[start_idx : end_idx + 1]

    if sequence_type == "gene":
        bases = ref_diffs_subset["Gene_Base"].tolist()
    else:
        bases = ref_diffs_subset["Pseudogene_Base"].tolist()
    # Remove gaps or missing bases represented by '-'
    sequence = "".join([base for base in bases if base != "-"])
    return sequence, start_idx


def generate_cigar_from_alignment(alignment):
    # Generate a CIGAR string from the alignment
    aligned_seqA = alignment.seqA  # Read sequence
    aligned_seqB = alignment.seqB  # Reference sequence
    cigar = []
    match = 0
    for a, b in zip(aligned_seqA, aligned_seqB):
        if a == "-" and b != "-":
            if match > 0:
                cigar.append(("M", match))
                match = 0
            cigar.append(("D", 1))
        elif a != "-" and b == "-":
            if match > 0:
                cigar.append(("M", match))
                match = 0
            cigar.append(("I", 1))
        elif a != b:
            if match > 0:
                cigar.append(("M", match))
                match = 0
            cigar.append(("X", 1))
        else:
            match += 1
    if match > 0:
        cigar.append(("M", match))
    # Summarize the CIGAR operations
    cigar_summary = {"M": 0, "I": 0, "D": 0, "X": 0}
    for op, length in cigar:
        cigar_summary[op] += length
    return cigar_summary


def cigar_dict_to_string(cigar_dict):
    if not cigar_dict:
        return None
    cigar_order = ["M", "I", "D", "X", "S"]
    cigar_str = ""
    for op in cigar_order:
        if op in cigar_dict and cigar_dict[op] > 0:
            cigar_str += f"{cigar_dict[op]}{op}"
    return cigar_str


def compare_cigars(cigar_1, cigar_2):
    # Compare two CIGAR summaries to decide which alignment is better
    # Returns 0 if they are the same, -1 if cigar_2 is better, and 1 if cigar_1 is better
    if cigar_1 is None or cigar_2 is None:
        return 1
    mismatches_1 = cigar_1.get("X", 0) + cigar_1.get("I", 0) + cigar_1.get("D", 0)
    mismatches_2 = cigar_2.get("X", 0) + cigar_2.get("I", 0) + cigar_2.get("D", 0)
    if mismatches_1 == mismatches_2:
        return 0
    elif mismatches_2 < mismatches_1:
        return -1
    else:
        return 1


def calculate_cigar_score(cigar_tuples):
    # Calculate a score for the CIGAR string (arbitrary for demonstration)
    score = 0
    for op, length in cigar_tuples:
        if op == 0:
            score += length * 1
        elif op == 1:
            score -= length * 2
        elif op == 2:
            score -= length * 2
        elif op == 8:
            score -= length * 1
    return score


def generate_reference_positions(alignment, start_pos):
    # Generate reference positions based on the alignment and starting position
    positions = []
    pos = start_pos
    aligned_seqA = alignment.seqA  # Read sequence (with gaps)
    aligned_seqB = alignment.seqB  # Reference sequence (with gaps)

    for a, b in zip(aligned_seqA, aligned_seqB):
        if b != "-":
            positions.append(pos)
            pos += 1
        else:
            positions.append(-1)  # Insertion in read; no reference position
    return positions


def create_detailed_table(ref_diffs_df, read_sequences, read_positions, alignment_results, all_reads, output_file):
    # Prepare the detailed table
    detailed_rows = []

    for idx, row in tqdm(ref_diffs_df.iterrows(), total=ref_diffs_df.shape[0], desc="Processing rows"):
        gene_pos = row["Gene_Position"]
        pseudo_pos = row["Pseudogene_Position"]
        diff_type = row["Difference_Type"]
        gene_pos_str = f"{row['Gene_Chromosome']}:{gene_pos}" if gene_pos != -1 else "-"
        pseudo_pos_str = f"{row['Pseudogene_Chromosome']}:{pseudo_pos}" if pseudo_pos != -1 else "-"
        data_row = {
            "Gene_Position": gene_pos_str,
            "Gene_Base": row["Gene_Base"],
            "Pseudogene_Position": pseudo_pos_str,
            "Pseudogene_Base": row["Pseudogene_Base"],
            "Difference_Type": diff_type,
        }
        for read in all_reads:
            base = "-"
            read_seq = read_sequences.get(read, "")
            read_pos_list = read_positions.get(read, [])
            if read_pos_list:
                for read_base, read_pos in zip(read_seq, read_pos_list):
                    if read_pos == gene_pos or read_pos == pseudo_pos:
                        base = read_base
            data_row[read] = base
        detailed_rows.append(data_row)

    # Create DataFrame
    detailed_df = pd.DataFrame(detailed_rows)
    # Arrange columns
    columns = ["Gene_Position", "Gene_Base", "Pseudogene_Position", "Pseudogene_Base", "Difference_Type"] + all_reads
    detailed_df = detailed_df[columns]
    # Save to file
    detailed_df.to_csv(output_file, sep="\t", index=False)
    print(f"Detailed table saved to {output_file}")
    # Output alignment results to a second file with suffix _alignment
    alignment_file = output_file.replace(".tsv", "_aln_results.tsv")
    alignment_rows = []
    for read, result in alignment_results.items():
        row = {
            "Read_Name": read,
            "Original_CIGAR": result.get("Original_CIGAR"),
            "Realigned_CIGAR": result.get("Realigned_CIGAR"),
            "Alignment_Score": result.get("Alignment_Score"),
        }
        alignment_rows.append(row)
    alignment_df = pd.DataFrame(
        alignment_rows,
        columns=["Read_Name", "Original_CIGAR", "Realigned_CIGAR", "Alignment_Score"],
    )
    alignment_df.to_csv(alignment_file, sep="\t", index=False)
    print(f"Alignment results saved to {alignment_file}")


if __name__ == "__main__":
    args = parse_arguments()
    ref_diffs_df = load_ref_diffs(args.tsv)
    process_reads(
        bam_file=args.bam,
        ref_diffs_df=ref_diffs_df,
        output=args.output,
        min_mapq=args.min_mapq,
        flank_size=0,
    )