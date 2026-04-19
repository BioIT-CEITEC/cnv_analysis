#!/usr/bin/env python3
#
# Usage: python 4_classify_reads.py --bam_original sample.bam
#                                  --realigned_dir realigned/
#                                  --bed regions.bed --out classifications.tsv

import pysam
import argparse
import pandas as pd
import os
from collections import defaultdict


def parse_bed(bed_file):
    pairs = defaultdict(list)
    with open(bed_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            chrom, start, end, name, pair_id = (
                parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4]
            )
            region_type = "pseudogene" if any(
                k in name.lower() for k in ["pseudogene", "pseudo", "_ps", "ps_"]
            ) else "gene"
            pairs[pair_id].append({
                "chrom": chrom, "start": start, "end": end,
                "name": name, "type": region_type, "pair_id": pair_id
            })
    return pairs


def get_original_mapping(bam_original, pairs):
    """
    Record original mapping metadata for each read.

    Only record reads whose OWN alignment overlaps the BED region,
    not reads where only the mate overlaps.
    """
    seen = {}
    with pysam.AlignmentFile(bam_original, "rb") as bam:
        for pair_id, regions in pairs.items():
            for region in regions:
                chrom = region["chrom"]
                start = region["start"]
                end   = region["end"]
                for read in bam.fetch(chrom, start, end):
                    if read.is_unmapped or read.query_name is None:
                        continue
                    # Filter: the read itself must overlap the region
                    if read.reference_name != chrom:
                        continue
                    if not (read.reference_start < end and read.reference_end > start):
                        continue
                    rname = read.query_name
                    if rname not in seen:
                        seen[rname] = {
                            "original_chrom":         read.reference_name,
                            "original_pos":           read.reference_start,
                            "original_mapq":          read.mapping_quality,
                            "original_mapped_region": region["name"],
                            "pair_id":                pair_id,
                        }
    return seen


def get_as_tag(read):
    try:
        return read.get_tag("AS")
    except KeyError:
        return -1


def collect_best_as_from_bam(bam_path):
    """
    Read a BAM aligned against ONE reference only (gene or pseudogene)
    and return the best AS per read_name.

    Since the reference is unique, AS is always available for each read.
    We use the maximum AS across all alignments of the read (primary + secondary,
    if present) to be conservative.

    Returns: {read_name: {"as": int, "mapq": int}}
    """
    best = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.query_name is None:
                continue
            rname  = read.query_name
            score  = get_as_tag(read)
            mapq   = read.mapping_quality
            is_primary = not read.is_secondary and not read.is_supplementary

            if rname not in best:
                best[rname] = {"as": score, "mapq": 0}
            else:
                if score > best[rname]["as"]:
                    best[rname]["as"] = score

            if is_primary:
                best[rname]["mapq"] = max(best[rname]["mapq"], mapq)

    return best


def classify_pair(pair_id, realigned_dir, original_mapping, mapq_threshold=20, min_as_delta=10):
    """
    Classify each fragment by comparing:
    - AS from the BAM aligned against ONLY the gene       (pair1_gene.bam)
    - AS from the BAM aligned against ONLY the pseudogene (pair1_pseudogene.bam)

    Since each read was aligned separately against each reference,
    AS is always available for both and should never be -1.

    Rules:
    - as_delta = as_gene - as_pseudo
    - as_delta >  min_as_delta -> confident gene
    - as_delta < -min_as_delta -> confident pseudogene
    - |as_delta| <= min_as_delta -> ambiguous (too similar)
    - mapq_gene=0 and mapq_pseudo=0 -> ambiguous (aligner cannot decide)

    misassigned = the read was in gene coordinates in the original BAM,
                  but as_pseudo > as_gene by a significant margin,
                  and the original alignment was reliable (original_mapq > 0)
    """
    bam_gene   = os.path.join(realigned_dir, f"{pair_id}_gene.bam")
    bam_pseudo = os.path.join(realigned_dir, f"{pair_id}_pseudogene.bam")

    if not os.path.exists(bam_gene) or not os.path.exists(bam_pseudo):
        print(f"[4_classify] WARNING: Missing BAMs for {pair_id}, skipping")
        return []

    scores_gene   = collect_best_as_from_bam(bam_gene)
    scores_pseudo = collect_best_as_from_bam(bam_pseudo)

    all_reads = set(scores_gene.keys()) | set(scores_pseudo.keys())

    records = []
    for rname in all_reads:
        orig = original_mapping.get(rname, {})

        as_gene   = scores_gene.get(rname,   {}).get("as",   -1)
        as_pseudo = scores_pseudo.get(rname, {}).get("as",   -1)
        mapq_gene = scores_gene.get(rname,   {}).get("mapq",  0)
        mapq_pseudo = scores_pseudo.get(rname, {}).get("mapq", 0)
        as_delta  = as_gene - as_pseudo

        original_mapq          = orig.get("original_mapq", -1)
        original_mapped_region = orig.get("original_mapped_region", "unknown")
        originally_in_gene     = "pseudogene" not in original_mapped_region.lower() and \
                                 not any(k in original_mapped_region.lower()
                                         for k in ["pseudo", "_ps", "ps_"])
        originally_in_pseudo   = not originally_in_gene
        originally_ambiguous   = (original_mapq == 0)

        both_mapq_zero = (mapq_gene == 0 and mapq_pseudo == 0)

        if both_mapq_zero and abs(as_delta) <= min_as_delta:
            best_ref    = "ambiguous"
            misassigned = False
            ambiguous   = True

        elif as_gene == -1 and as_pseudo == -1:
            best_ref    = "unmapped"
            misassigned = False
            ambiguous   = True

        elif as_gene == -1:
            best_ref    = "pseudogene"
            ambiguous   = mapq_pseudo < mapq_threshold
            misassigned = originally_in_gene and not ambiguous and not originally_ambiguous

        elif as_pseudo == -1:
            best_ref    = "gene"
            ambiguous   = mapq_gene < mapq_threshold
            misassigned = originally_in_pseudo and not ambiguous and not originally_ambiguous

        elif as_delta > min_as_delta:
            best_ref    = "gene"
            ambiguous   = originally_ambiguous
            misassigned = originally_in_pseudo and not originally_ambiguous

        elif as_delta < -min_as_delta:
            best_ref    = "pseudogene"
            ambiguous   = originally_ambiguous
            misassigned = originally_in_gene and not originally_ambiguous

        else:
            best_ref    = "ambiguous"
            misassigned = False
            ambiguous   = True

        records.append({
            "read_name":             rname,
            "pair_id":               pair_id,
            "original_chrom":        orig.get("original_chrom", "unknown"),
            "original_pos":          orig.get("original_pos", -1),
            "original_mapq":         original_mapq,
            "original_mapped_region": original_mapped_region,
            "as_gene":               as_gene,
            "mapq_gene":             mapq_gene,
            "as_pseudogene":         as_pseudo,
            "mapq_pseudogene":       mapq_pseudo,
            "as_delta":              as_delta,
            "best_ref":              best_ref,
            "misassigned":           misassigned,
            "ambiguous":             ambiguous,
        })

    return records


def classify_all_pairs(bam_original, realigned_dir, bed_path, out_tsv,
                       mapq_threshold=20, min_as_delta=10):
    pairs = parse_bed(bed_path)

    print("[4_classify] Recording original mapping...")
    original_mapping = get_original_mapping(bam_original, pairs)
    print(f"[4_classify] Recorded reads: {len(original_mapping)}")

    all_records = []
    for pair_id in pairs:
        print(f"\n[4_classify] Classifying pair {pair_id}...")
        records = classify_pair(
            pair_id, realigned_dir, original_mapping,
            mapq_threshold, min_as_delta
        )
        all_records.extend(records)

        total       = len(records)
        if total == 0:
            continue
        misassigned = sum(r["misassigned"] for r in records)
        ambiguous   = sum(r["ambiguous"]   for r in records)
        pct_mis = 100 * misassigned / total
        pct_amb = 100 * ambiguous   / total
        print(f"[4_classify]   total={total} | "
              f"misassigned={misassigned} ({pct_mis:.1f}%) | "
              f"ambiguous={ambiguous} ({pct_amb:.1f}%)")

    if not all_records:
        print("\n[4_classify] ERROR: No reads were classified.")
        bams = [f for f in os.listdir(realigned_dir)
                if f.endswith(".bam")] if os.path.isdir(realigned_dir) else []
        print(f"  BAMs in --realigned_dir: {bams}")
        return

    df = pd.DataFrame(all_records)
    df.to_csv(out_tsv, sep="\t", index=False)

    total       = len(df)
    misassigned = df["misassigned"].sum()
    ambiguous   = df["ambiguous"].sum()
    print(f"\n[4_classify] ══════════════════════════════════════════")
    print(f"[4_classify] TOTAL classified reads            : {total}")
    print(f"[4_classify] Misassigned (gene<->pseudogene)   : {misassigned} ({100*misassigned/total:.1f}%)")
    print(f"[4_classify] Ambiguous                         : {ambiguous} ({100*ambiguous/total:.1f}%)")
    print(f"[4_classify] Done -> {out_tsv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Classify reads by comparing AS gene vs pseudogene (two separate BAMs)"
    )
    parser.add_argument("--bam_original",  required=True)
    parser.add_argument("--realigned_dir", required=True)
    parser.add_argument("--bed",           required=True)
    parser.add_argument("--out",           required=True)
    parser.add_argument("--mapq",          type=int, default=20,
                        help="Minimum MAPQ threshold for confident classification (default: 20)")
    parser.add_argument("--min_as_delta",  type=int, default=10,
                        help="Minimum AS difference to consider one alignment better (default: 10)")
    args = parser.parse_args()

    classify_all_pairs(
        args.bam_original, args.realigned_dir, args.bed, args.out,
        mapq_threshold=args.mapq, min_as_delta=args.min_as_delta
    )
