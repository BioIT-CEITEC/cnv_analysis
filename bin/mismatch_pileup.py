#!/usr/bin/env python3
"""
For each mismatch position (Difference_Type=X) in the gene/pseudogene diff TSV,
count per-nucleotide read coverage from the realigned gene.bam and pseudogene.bam.

Compatible with 3_realign_specific.py after the liftover fix:
  - BAM RNAME = real chromosome name (e.g. chr7), not the gene name
  - BAM POS   = true genomic coordinates (1-based in SAM, 0-based in pysam)
  - diff TSV positions are 1-based genomic coords

Output TSV columns per mismatch position:
  sample, pair_id,
  gene_chrom, gene_pos, gene_base,
  pseudo_chrom, pseudo_pos, pseudo_base,
  g_bam_{A,C,G,T,N}, g_bam_depth, g_bam_gene_frac, g_bam_pseudo_frac,
  p_bam_{A,C,G,T,N}, p_bam_depth, p_bam_gene_frac, p_bam_pseudo_frac
"""

import argparse
import csv
import os
from collections import defaultdict

import pysam


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--diff_dir",         required=True,
                   help="Directory with {pair_id}.tsv files from ALIGN_REGIONS")
    p.add_argument("--realigned_dir",    required=True,
                   help="Directory with {pair_id}_gene.bam and {pair_id}_pseudogene.bam")
    p.add_argument("--bed",              required=True,
                   help="BED file (chrom start end name pair_id [strand ...])")
    p.add_argument("--sample",           required=True,  help="Sample name")
    p.add_argument("--output",           required=True,  help="Output TSV path")
    p.add_argument("--min_base_quality", type=int, default=20,
                   help="Minimum base quality for pileup (default: 20)")
    p.add_argument("--min_map_quality",  type=int, default=10,
                   help="Minimum mapping quality for pileup (default: 10)")
    return p.parse_args()


def _infer_region_type(name, extra_fields=()):
    """Mirror align.py logic: detect gene vs pseudogene from name/annotations."""
    candidates = [name] + list(extra_fields)
    pseudo_keys     = ("pseudogene", "pseudo", "_ps", "ps_")
    pseudo_suffixes = ("cl", "_cl", "_b", "_dup", "_v", "_v2", "_v3")
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
    """Return {pair_id: {"gene": region_dict, "pseudogene": region_dict}}.

    Uses the same region-type detection and positional fallback as align.py.
    """
    raw = defaultdict(lambda: {"gene": None, "pseudogene": None, "unknown": []})
    with open(bed_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            chrom, start, end, name, pair_id = (
                parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4]
            )
            extra = parts[6:] if len(parts) > 6 else []
            rtype = _infer_region_type(name, extra)
            region = {"chrom": chrom, "start": start, "end": end, "name": name}
            if rtype == "gene" and raw[pair_id]["gene"] is None:
                raw[pair_id]["gene"] = region
            elif rtype == "pseudogene" and raw[pair_id]["pseudogene"] is None:
                raw[pair_id]["pseudogene"] = region
            else:
                raw[pair_id]["unknown"].append(region)

    pairs = {}
    for pair_id, entry in raw.items():
        gene, pseudo, unknown = entry["gene"], entry["pseudogene"], list(entry["unknown"])
        if gene is None and unknown:
            gene = unknown.pop(0)
        if pseudo is None and unknown:
            pseudo = unknown.pop(0)
        pairs[pair_id] = {"gene": gene, "pseudogene": pseudo}
    return pairs


def parse_position(pos_str):
    """'chr7:1001' -> ('chr7', 1001); '.' -> None."""
    if pos_str == ".":
        return None
    chrom, pos = pos_str.rsplit(":", 1)
    return chrom, int(pos)


def pileup_at(bam_path, chrom, pos_1based, min_bq, min_mq):
    """
    Return nucleotide counts and depth at a genomic position.

    Parameters
    ----------
    bam_path  : path to the liftover-corrected BAM (true genomic coordinates)
    chrom     : real chromosome name, e.g. 'chr7'  (RNAME in the BAM header)
    pos_1based: 1-based genomic position (as reported in the diff TSV)
    min_bq    : minimum base quality
    min_mq    : minimum mapping quality
    """
    counts = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
    if not os.path.exists(bam_path):
        return counts | {"depth": 0}

    # pysam pileup uses 0-based half-open coordinates
    pos_0based = pos_1based - 1

    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for col in bam.pileup(
                chrom,
                pos_0based,
                pos_0based + 1,
                min_base_quality=min_bq,
                truncate=True,
                stepper="nofilter",
            ):
                if col.reference_pos != pos_0based:
                    continue
                for pread in col.pileups:
                    if pread.is_del or pread.is_refskip:
                        continue
                    if pread.alignment.mapping_quality < min_mq:
                        continue
                    base = pread.alignment.query_sequence[pread.query_position].upper()
                    counts[base] = counts.get(base, 0) + 1
                break
    except (ValueError, KeyError):
        pass

    return counts | {"depth": sum(counts.values())}


def frac(count, depth):
    return round(count / depth, 4) if depth > 0 else 0.0


def main():
    args = parse_args()
    pairs = parse_bed(args.bed)
    out_rows = []

    for pair_id, regions in sorted(pairs.items()):
        gene_info   = regions["gene"]
        pseudo_info = regions["pseudogene"]
        if gene_info is None or pseudo_info is None:
            print(f"[WARNING] {pair_id}: incomplete pair definition, skipping")
            continue

        diff_tsv   = os.path.join(args.diff_dir,      f"{pair_id}.tsv")
        gene_bam   = os.path.join(args.realigned_dir, f"{pair_id}_gene.bam")
        pseudo_bam = os.path.join(args.realigned_dir, f"{pair_id}_pseudogene.bam")

        if not os.path.exists(diff_tsv):
            print(f"[WARNING] {pair_id}: diff TSV not found ({diff_tsv}), skipping")
            continue

        # After the liftover fix, RNAME in the BAM is the real chromosome,
        # not the gene name. Positions are true genomic coords (1-based in TSV).
        gene_chrom_rname   = gene_info["chrom"]
        pseudo_chrom_rname = pseudo_info["chrom"]

        with open(diff_tsv) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if row["Difference_Type"] != "X":
                    continue

                gene_coord   = parse_position(row["Gene_Position"])
                pseudo_coord = parse_position(row["Pseudogene_Position"])
                if gene_coord is None or pseudo_coord is None:
                    continue

                gene_chrom,   gene_pos   = gene_coord
                pseudo_chrom, pseudo_pos = pseudo_coord
                gene_base   = row["Gene_Base"].upper()
                pseudo_base = row["Pseudogene_Base"].upper()

                # Sanity check: chrom in TSV should match BED
                if gene_chrom != gene_chrom_rname:
                    print(
                        f"[WARNING] {pair_id}: gene chrom mismatch — "
                        f"TSV={gene_chrom}, BED={gene_chrom_rname}"
                    )
                if pseudo_chrom != pseudo_chrom_rname:
                    print(
                        f"[WARNING] {pair_id}: pseudo chrom mismatch — "
                        f"TSV={pseudo_chrom}, BED={pseudo_chrom_rname}"
                    )

                # Query BAM directly with genomic coordinates — no offset arithmetic needed
                g = pileup_at(gene_bam,   gene_chrom_rname,   gene_pos,
                              args.min_base_quality, args.min_map_quality)
                p = pileup_at(pseudo_bam, pseudo_chrom_rname, pseudo_pos,
                              args.min_base_quality, args.min_map_quality)

                out_rows.append({
                    "sample":            args.sample,
                    "pair_id":           pair_id,
                    "gene_chrom":        gene_chrom,
                    "gene_pos":          gene_pos,
                    "gene_base":         gene_base,
                    "pseudo_chrom":      pseudo_chrom,
                    "pseudo_pos":        pseudo_pos,
                    "pseudo_base":       pseudo_base,
                    "g_bam_A":           g["A"],
                    "g_bam_C":           g["C"],
                    "g_bam_G":           g["G"],
                    "g_bam_T":           g["T"],
                    "g_bam_N":           g.get("N", 0),
                    "g_bam_depth":       g["depth"],
                    "g_bam_gene_frac":   frac(g.get(gene_base,   0), g["depth"]),
                    "g_bam_pseudo_frac": frac(g.get(pseudo_base, 0), g["depth"]),
                    "p_bam_A":           p["A"],
                    "p_bam_C":           p["C"],
                    "p_bam_G":           p["G"],
                    "p_bam_T":           p["T"],
                    "p_bam_N":           p.get("N", 0),
                    "p_bam_depth":       p["depth"],
                    "p_bam_gene_frac":   frac(p.get(gene_base,   0), p["depth"]),
                    "p_bam_pseudo_frac": frac(p.get(pseudo_base, 0), p["depth"]),
                })

    fieldnames = [
        "sample", "pair_id",
        "gene_chrom", "gene_pos", "gene_base",
        "pseudo_chrom", "pseudo_pos", "pseudo_base",
        "g_bam_A", "g_bam_C", "g_bam_G", "g_bam_T", "g_bam_N",
        "g_bam_depth", "g_bam_gene_frac", "g_bam_pseudo_frac",
        "p_bam_A", "p_bam_C", "p_bam_G", "p_bam_T", "p_bam_N",
        "p_bam_depth", "p_bam_gene_frac", "p_bam_pseudo_frac",
    ]

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"[INFO] {args.sample}: wrote {len(out_rows)} mismatch positions -> {args.output}")


if __name__ == "__main__":
    main()