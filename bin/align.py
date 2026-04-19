#!/usr/bin/env python3

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile

import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="For each paired gene/pseudogene region in a BED file, extract FASTAs from the "
                    "reference, align with minimap2, and write a per-pair base-level diff TSV."
    )
    parser.add_argument("--reference", required=True, help="Full-genome reference FASTA (must have .fai index)")
    parser.add_argument("--bed", required=True,
                        help="BED file with paired regions; columns: chrom start end name pair_id [strand ...]")
    parser.add_argument("--output_dir", required=True, help="Directory for per-pair TSV (and optional SAM) output")
    parser.add_argument("--preset", default="asm5", help="minimap2 preset (default: asm5)")
    parser.add_argument("--threads", type=int, default=1, help="Threads for minimap2")
    parser.add_argument("--keep_sam", action="store_true", help="Keep intermediate SAM files")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# BED parsing
# ---------------------------------------------------------------------------

def _infer_region_type(name, extra_fields):
    """Infer region type from the name or optional BED annotation columns.
    
    Recognizes pseudogene patterns: explicit labels, CL/B/_dup/_v suffixes.
    Recognizes gene: explicit 'gene' label.
    Falls back to None if ambiguous.
    """
    candidates = [name] + list(extra_fields)
    pseudo_keys = ("pseudogene", "pseudo", "_ps", "ps_")
    pseudo_suffixes = ("cl", "_cl", "_b", "_dup", "_v", "_v2", "_v3")  # Copy-like paralogs

    for value in candidates:
        lv = value.lower()
        if any(k in lv for k in pseudo_keys):
            return "pseudogene"
        if any(lv.endswith(s) for s in pseudo_suffixes):
            return "pseudogene"
        if "gene" in lv and "pseudogene" not in lv:
            return "gene"
    return None


def parse_bed(bed_file):
    """Return {pair_id: {"gene": {...}, "pseudogene": {...}}} from a BED file.

    Required columns: chrom(0) start(1) end(2) name(3) pair_id(4)
    Optional columns: strand(5), region-type annotations (6+)

    If labels are missing but exactly two regions are present for a pair, the
    first region is treated as gene and the second as pseudogene.
    """
    parsed_pairs = {}

    with open(bed_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 5:
                raise ValueError(f"Invalid BED line (expected >=5 columns): {line}")

            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            name = parts[3]
            pair_id = parts[4]
            strand = parts[5] if len(parts) > 5 else "-"
            extra_fields = parts[6:] if len(parts) > 6 else []
            region_type = _infer_region_type(name, extra_fields)

            region = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "name": name,
                "strand": strand,
            }

            if pair_id not in parsed_pairs:
                parsed_pairs[pair_id] = {"gene": None, "pseudogene": None, "unknown": []}

            if region_type == "gene":
                parsed_pairs[pair_id]["gene"] = region
            elif region_type == "pseudogene":
                parsed_pairs[pair_id]["pseudogene"] = region
            else:
                parsed_pairs[pair_id]["unknown"].append(region)

    resolved_pairs = {}
    for pair_id, entry in parsed_pairs.items():
        gene = entry["gene"]
        pseudogene = entry["pseudogene"]
        unknown = entry["unknown"]

        # Fill empty slots from unknowns
        if gene is None and unknown:
            gene = unknown.pop(0)
        if pseudogene is None and unknown:
            pseudogene = unknown.pop(0)

        # Warn if we had to use positional fallback
        if (gene is not None or pseudogene is not None) and (entry["gene"] is None or entry["pseudogene"] is None):
            print(
                f"[WARNING] {pair_id}: region types not explicit in BED; "
                "assigned by position (gene first, pseudogene second)",
                file=sys.stderr,
            )

        resolved_pairs[pair_id] = {"gene": gene, "pseudogene": pseudogene}

    return resolved_pairs


# ---------------------------------------------------------------------------
# Reference FASTA extraction
# ---------------------------------------------------------------------------

def revcomp(seq):
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]


def extract_fasta(reference, region, out_path):
    """Fetch region from reference and write a single-record FASTA.

    region dict: chrom, start (0-based), end (exclusive), name, strand

    FIX: If the region is on the minus strand, the raw forward-strand sequence
    is reverse-complemented so that minimap2 receives both sequences in their
    biological 5'→3' orientation. This prevents the query (pseudogene) from
    appearing reverse-complemented in the alignment, which would cause
    pseudogene coordinates to be reported in decreasing order.
    """
    with pysam.FastaFile(reference) as fa:
        seq = fa.fetch(region["chrom"], region["start"], region["end"]).upper()

    if region["strand"] == "-":
        seq = revcomp(seq)
        print(f"[INFO]   reverse-complemented {region['name']} (strand -)", file=sys.stderr)

    with open(out_path, "w") as f:
        f.write(f">{region['name']} {region['chrom']}:{region['start']}-{region['end']} strand:{region['strand']}\n")
        f.write(seq + "\n")
    return seq


# ---------------------------------------------------------------------------
# minimap2 helpers
# ---------------------------------------------------------------------------

def cigar_ref_bases(cigar):
    total = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op in ("M", "D", "N", "=", "X"):
            total += int(length)
    return total


def choose_best_alignment(sam_path):
    best = None
    total_records = 0
    mapped_records = 0
    primary_mapped_records = 0
    with open(sam_path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            total_records += 1
            flag = int(fields[1])
            if flag & 0x4:
                continue
            mapped_records += 1
            if flag & 0x100 or flag & 0x800:
                continue
            primary_mapped_records += 1
            cigar = fields[5]
            mapq = int(fields[4])
            score = (mapq, cigar_ref_bases(cigar))
            if best is None or score > best["score"]:
                best = {"line": fields, "score": score}
    if best is None:
        raise RuntimeError(
            "No primary mapped alignment found in SAM "
            f"(records={total_records}, mapped={mapped_records}, primary_mapped={primary_mapped_records})."
        )
    return best["line"]


def run_minimap2(gene_fasta, pseudo_fasta, preset, threads, sam_path):
    cmd = [
        "minimap2", "-a", "--eqx", "--secondary=no", "-N", "1",
        "-x", preset, "-t", str(threads),
        gene_fasta, pseudo_fasta,
    ]
    with open(sam_path, "w") as out:
        result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        raise RuntimeError("minimap2 failed")


# ---------------------------------------------------------------------------
# CIGAR walking → diff rows
# ---------------------------------------------------------------------------

def build_diff_rows(gene_seq, pseudo_seq, gene_info, pseudo_info, sam_path):
    """Walk the CIGAR and emit one row per aligned base.

    FIX: Because extract_fasta() now reverse-complements minus-strand sequences
    before alignment, minimap2 should always produce a forward alignment
    (is_reverse=False). The coordinate reconstruction is updated accordingly:

    - Gene (minus strand): the sequence passed to minimap2 was the RC of the
      genomic forward strand, so base index 0 corresponds to genomic position
      gene_end and index i corresponds to gene_end - i (1-based).
    - Pseudogene (plus strand): unchanged — index i → pseudo_start + i + 1.

    If for any reason minimap2 still returns a reverse alignment (is_reverse=True),
    a warning is emitted and the original fallback coordinate logic is used so
    the output is never silently wrong.
    """
    fields = choose_best_alignment(sam_path)

    flag = int(fields[1])
    pos_1based = int(fields[3])
    mapq = int(fields[4])
    cigar = fields[5]

    is_reverse = bool(flag & 0x10)

    if is_reverse:
        print(
            "[WARNING] Alignment is reverse — gene sequence may not have been "
            "reverse-complemented correctly. Coordinates may be decreasing.",
            file=sys.stderr,
        )

    # After the fix, gene_seq is already in 5'→3' orientation (RC if strand=="-")
    # and pseudo_seq is likewise in its correct orientation.
    query_seq_oriented = revcomp(pseudo_seq) if is_reverse else pseudo_seq

    gene_chr = gene_info["chrom"]
    gene_start = gene_info["start"]
    gene_end = gene_info["end"]
    gene_strand = gene_info["strand"]

    pseudo_chr = pseudo_info["chrom"]
    pseudo_start = pseudo_info["start"]
    pseudo_end = pseudo_info["end"]
    pseudo_strand = pseudo_info["strand"]

    ref_idx = pos_1based - 1
    query_idx = 0
    rows = []

    for length_str, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
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
                dtype = "M" if (op == "=" or gbase == pbase) else "X"

                # Gene coordinate:
                # - minus strand: index 0 of the RC sequence = genomic gene_end (1-based)
                # - plus  strand: index 0 = genomic gene_start + 1 (1-based)
                if gene_strand == "-":
                    gene_pos = f"{gene_chr}:{gene_end - ref_idx}"
                else:
                    gene_pos = f"{gene_chr}:{gene_start + ref_idx + 1}"

                # Pseudogene coordinate:
                # - plus  strand: straightforward
                # - minus strand (unusual): mirror logic
                if pseudo_strand == "-":
                    pseudo_pos = f"{pseudo_chr}:{pseudo_end - query_idx}"
                else:
                    pseudo_pos = f"{pseudo_chr}:{pseudo_start + query_idx + 1}"

                rows.append({
                    "Gene_Position": gene_pos, "Gene_Base": gbase,
                    "Pseudogene_Position": pseudo_pos, "Pseudogene_Base": pbase,
                    "Difference_Type": dtype,
                })
                ref_idx += 1
                query_idx += 1

        elif op == "I":
            for _ in range(length):
                pbase = query_seq_oriented[query_idx]
                if pseudo_strand == "-":
                    pseudo_pos = f"{pseudo_chr}:{pseudo_end - query_idx}"
                else:
                    pseudo_pos = f"{pseudo_chr}:{pseudo_start + query_idx + 1}"
                rows.append({
                    "Gene_Position": ".", "Gene_Base": "-",
                    "Pseudogene_Position": pseudo_pos, "Pseudogene_Base": pbase,
                    "Difference_Type": "I",
                })
                query_idx += 1

        elif op in ("D", "N"):
            for _ in range(length):
                gbase = gene_seq[ref_idx]
                if gene_strand == "-":
                    gene_pos = f"{gene_chr}:{gene_end - ref_idx}"
                else:
                    gene_pos = f"{gene_chr}:{gene_start + ref_idx + 1}"
                rows.append({
                    "Gene_Position": gene_pos, "Gene_Base": gbase,
                    "Pseudogene_Position": ".", "Pseudogene_Base": "-",
                    "Difference_Type": "D",
                })
                ref_idx += 1

        else:
            raise RuntimeError(f"Unsupported CIGAR op: {op}")

    print(f"[INFO]   strand: {'reverse' if is_reverse else 'forward'} | MAPQ: {mapq} | CIGAR: {cigar} | rows: {len(rows)}")
    return rows


def write_tsv(rows, output_tsv):
    with open(output_tsv, "w", newline="") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=["Gene_Position", "Gene_Base", "Pseudogene_Position", "Pseudogene_Base", "Difference_Type"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    pairs = parse_bed(args.bed)

    for pair_id, regions in pairs.items():
        gene_info = regions.get("gene")
        pseudo_info = regions.get("pseudogene")
        if gene_info is None or pseudo_info is None:
            print(f"[WARNING] {pair_id}: missing gene or pseudogene region (gene={gene_info is not None}, pseudo={pseudo_info is not None}), skipping", file=sys.stderr)
            continue

        print(f"\n[INFO] Pair: {pair_id}")
        print(f"[INFO]   gene:       {gene_info['chrom']}:{gene_info['start']}-{gene_info['end']} ({gene_info['strand']})")
        print(f"[INFO]   pseudogene: {pseudo_info['chrom']}:{pseudo_info['start']}-{pseudo_info['end']} ({pseudo_info['strand']})")

        output_tsv = os.path.join(args.output_dir, f"{pair_id}.tsv")

        with tempfile.TemporaryDirectory() as tmpdir:
            gene_fasta = os.path.join(tmpdir, "gene.fasta")
            pseudo_fasta = os.path.join(tmpdir, "pseudogene.fasta")

            gene_seq = extract_fasta(args.reference, gene_info, gene_fasta)
            pseudo_seq = extract_fasta(args.reference, pseudo_info, pseudo_fasta)

            sam_path = (
                os.path.join(args.output_dir, f"{pair_id}.sam")
                if args.keep_sam
                else os.path.join(tmpdir, f"{pair_id}.sam")
            )

            run_minimap2(gene_fasta, pseudo_fasta, args.preset, args.threads, sam_path)
            try:
                rows = build_diff_rows(gene_seq, pseudo_seq, gene_info, pseudo_info, sam_path)
            except RuntimeError as e:
                if "No primary mapped alignment found in SAM" in str(e):
                    print(f"[WARNING] {pair_id}: {e} Skipping pair.", file=sys.stderr)
                    continue
                raise
            write_tsv(rows, output_tsv)
            print(f"[INFO]   output: {output_tsv}")

    print(f"\n[INFO] Done. TSVs in: {args.output_dir}/")


if __name__ == "__main__":
    main()