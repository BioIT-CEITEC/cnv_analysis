#!/usr/bin/env python3

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile


VALID_DNA = set("ACGTN")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Align gene vs pseudogene FASTA using minimap2 and export base-level differences with genomic coordinates."
    )
    parser.add_argument("--gene_fasta", required=True, help="Reference FASTA for gene region")
    parser.add_argument("--pseudogene_fasta", required=True, help="Query FASTA for pseudogene region")
    parser.add_argument("--gene_region", required=True, help="Gene region, e.g. 7:5970925-6009130")
    parser.add_argument("--pseudogene_region", required=True, help="Pseudogene region, e.g. 7:6735305-6751392")
    parser.add_argument("--output", required=True, help="Output TSV")
    parser.add_argument("--preset", default="asm5", help="minimap2 preset (default: asm5)")
    parser.add_argument("--threads", type=int, default=1, help="Threads for minimap2")
    parser.add_argument("--keep_sam", action="store_true", help="Keep intermediate SAM file")
    return parser.parse_args()


def parse_region(region_str):
    m = re.match(r"^(?:GRCh38_)?(?:chr)?([^:]+):(\d+)-(\d+)$", region_str)
    if not m:
        raise ValueError(f"Invalid region format: {region_str}")
    chrom = m.group(1)
    start = int(m.group(2))
    end = int(m.group(3))
    if start > end:
        raise ValueError(f"Region start > end: {region_str}")
    return chrom, start, end


def read_single_fasta(path):
    header = None
    seq_parts = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    raise ValueError(f"{path} contains more than one FASTA record")
                header = line[1:].split()[0]
            else:
                seq_parts.append(line.upper())

    if header is None:
        raise ValueError(f"{path} does not look like a FASTA file")

    seq = "".join(seq_parts)
    invalid = sorted(set(seq) - VALID_DNA)
    if invalid:
        raise ValueError(
            f"{path} contains invalid sequence characters: {invalid}. Expected only A/C/G/T/N."
        )

    return header, seq


def revcomp(seq):
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]


def cigar_ref_bases(cigar):
    total = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length)
        if op in ("M", "D", "N", "=", "X"):
            total += length
    return total


def cigar_query_bases(cigar):
    total = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length)
        if op in ("M", "I", "S", "=", "X"):
            total += length
    return total


def choose_best_alignment(sam_path):
    best = None

    with open(sam_path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue

            flag = int(fields[1])

            # skip unmapped, secondary, supplementary
            if flag & 0x4 or flag & 0x100 or flag & 0x800:
                continue

            cigar = fields[5]
            mapq = int(fields[4])
            ref_consumed = cigar_ref_bases(cigar)

            score = (mapq, ref_consumed)
            if best is None or score > best["score"]:
                best = {
                    "line": fields,
                    "score": score
                }

    if best is None:
        raise RuntimeError("No primary mapped alignment found in SAM.")

    return best["line"]


def run_minimap2(gene_fasta, pseudo_fasta, preset, threads, sam_path):
    cmd = [
        "minimap2",
        "-a",
        "--eqx",
        "--secondary=no",
        "-N", "1",
        "-x", preset,
        "-t", str(threads),
        gene_fasta,
        pseudo_fasta
    ]

    with open(sam_path, "w") as out:
        result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        raise RuntimeError("minimap2 failed")

    return cmd


def main():
    args = parse_args()

    gene_id, gene_seq = read_single_fasta(args.gene_fasta)
    pseudo_id, pseudo_seq = read_single_fasta(args.pseudogene_fasta)

    gene_chr, gene_start, gene_end = parse_region(args.gene_region)
    pseudo_chr, pseudo_start, pseudo_end = parse_region(args.pseudogene_region)

    if len(gene_seq) != gene_end - gene_start + 1:
        raise ValueError(
            f"Gene FASTA length ({len(gene_seq)}) does not match region length ({gene_end - gene_start + 1})"
        )

    if len(pseudo_seq) != pseudo_end - pseudo_start + 1:
        raise ValueError(
            f"Pseudogene FASTA length ({len(pseudo_seq)}) does not match region length ({pseudo_end - pseudo_start + 1})"
        )

    if args.keep_sam:
        sam_path = os.path.abspath(args.output + ".sam")
        run_minimap2(args.gene_fasta, args.pseudogene_fasta, args.preset, args.threads, sam_path)
    else:
        with tempfile.NamedTemporaryFile(suffix=".sam", delete=False) as tmp:
            sam_path = tmp.name
        run_minimap2(args.gene_fasta, args.pseudogene_fasta, args.preset, args.threads, sam_path)

    fields = choose_best_alignment(sam_path)

    qname = fields[0]
    flag = int(fields[1])
    rname = fields[2]
    pos_1based = int(fields[3])
    mapq = int(fields[4])
    cigar = fields[5]
    seq_from_sam = fields[9].upper()

    is_reverse = bool(flag & 0x10)

    # Use original pseudogene FASTA sequence, but orient it according to alignment strand
    query_seq_oriented = revcomp(pseudo_seq) if is_reverse else pseudo_seq

    # In principle SEQ in SAM should match the query, but to avoid ambiguity across tools,
    # we rely on the FASTA we provided.
    ref_idx = pos_1based - 1  # 0-based index within gene_seq
    query_idx = 0

    rows = []

    cigar_ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar)

    for length_str, op in cigar_ops:
        length = int(length_str)

        if op == "H":
            continue

        if op == "S":
            query_idx += length
            continue

        if op in ("=", "X", "M"):
            for _ in range(length):
                gbase = gene_seq[ref_idx]
                pbase = query_seq_oriented[query_idx]

                if op == "=":
                    dtype = "M"
                elif op == "X":
                    dtype = "X"
                else:
                    dtype = "M" if gbase == pbase else "X"

                gene_pos = f"{gene_chr}:{gene_start + ref_idx}"
                if is_reverse:
                    pseudo_pos = f"{pseudo_chr}:{pseudo_end - query_idx}"
                else:
                    pseudo_pos = f"{pseudo_chr}:{pseudo_start + query_idx}"

                rows.append({
                    "Gene_Position": gene_pos,
                    "Gene_Base": gbase,
                    "Pseudogene_Position": pseudo_pos,
                    "Pseudogene_Base": pbase,
                    "Difference_Type": dtype
                })

                ref_idx += 1
                query_idx += 1

        elif op == "I":
            for _ in range(length):
                pbase = query_seq_oriented[query_idx]
                if is_reverse:
                    pseudo_pos = f"{pseudo_chr}:{pseudo_end - query_idx}"
                else:
                    pseudo_pos = f"{pseudo_chr}:{pseudo_start + query_idx}"

                rows.append({
                    "Gene_Position": ".",
                    "Gene_Base": "-",
                    "Pseudogene_Position": pseudo_pos,
                    "Pseudogene_Base": pbase,
                    "Difference_Type": "I"
                })
                query_idx += 1

        elif op in ("D", "N"):
            # deletion in pseudogene relative to gene
            for _ in range(length):
                gbase = gene_seq[ref_idx]
                gene_pos = f"{gene_chr}:{gene_start + ref_idx}"

                rows.append({
                    "Gene_Position": gene_pos,
                    "Gene_Base": gbase,
                    "Pseudogene_Position": ".",
                    "Pseudogene_Base": "-",
                    "Difference_Type": "D"
                })
                ref_idx += 1

        else:
            raise RuntimeError(f"Unsupported CIGAR op: {op}")

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=[
                "Gene_Position",
                "Gene_Base",
                "Pseudogene_Position",
                "Pseudogene_Base",
                "Difference_Type"
            ],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"[INFO] gene fasta: {args.gene_fasta}")
    print(f"[INFO] pseudogene fasta: {args.pseudogene_fasta}")
    print(f"[INFO] alignment strand: {'reverse' if is_reverse else 'forward'}")
    print(f"[INFO] MAPQ: {mapq}")
    print(f"[INFO] CIGAR: {cigar}")
    print(f"[INFO] rows written: {len(rows)}")
    print(f"[INFO] output: {args.output}")

    if not args.keep_sam and os.path.exists(sam_path):
        os.unlink(sam_path)


if __name__ == "__main__":
    main()
