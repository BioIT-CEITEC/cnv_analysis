#!/usr/bin/env python3

import argparse
import csv
import os
from collections import defaultdict

import pysam

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def complement(base):
    return base.translate(_COMPLEMENT)


def normalize_base(base, gene_base_plus, pseudo_base_plus):
    """
    Map complement(gene_base_plus) → gene_base_plus and
    complement(pseudo_base_plus) → pseudo_base_plus when unambiguous.
    Handles BWA FR-pairing where R2 reads store complement(true_base).
    """
    if base == ".":
        return base
    c_gene   = complement(gene_base_plus)
    c_pseudo = complement(pseudo_base_plus)
    if base == c_gene and c_gene != pseudo_base_plus:
        return gene_base_plus
    if base == c_pseudo and c_pseudo != gene_base_plus:
        return pseudo_base_plus
    return base


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--diff_dir",         required=True,
                   help="Directory with {pair_id}.tsv files from ALIGN_REGIONS")
    p.add_argument("--realigned_dir",    required=True,
                   help="Directory with {pair_id}_gene.bam and {pair_id}_pseudogene.bam")
    p.add_argument("--bam_original",     required=True,
                   help="Original/main BAM file (true genomic coordinates)")
    p.add_argument("--bed",              required=True,
                   help="BED file (chrom start end name pair_id [strand ...])")
    p.add_argument("--sample",           required=True, help="Sample name")
    p.add_argument("--output",           required=True, help="Output TSV path")
    p.add_argument("--min_base_quality", type=int, default=20,
                   help="Minimum base quality (default: 20)")
    p.add_argument("--min_map_quality",  type=int, default=10,
                   help="Minimum mapping quality (default: 10)")
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
    """Return {pair_id: {"gene": region_dict, "pseudogene": region_dict}}."""
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
            strand = parts[5] if len(parts) > 5 else "+"
            extra  = parts[6:] if len(parts) > 6 else []
            rtype  = _infer_region_type(name, extra)
            region = {"chrom": chrom, "start": start, "end": end, "name": name, "strand": strand}
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


def bases_at(bam_path, chrom, pos_1based, min_bq, min_mq):
    """
    Return {read_name: base} for all reads spanning a genomic position.

    Uses the same coordinate conventions as pileup_at() in mismatch_pileup.py.
    Skips deletions, ref-skips, and reads below quality thresholds.
    When a read appears multiple times (supplementary alignments), the first
    observed base is kept.
    """
    result = {}
    if not os.path.exists(bam_path):
        return result

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
                    rname = pread.alignment.query_name
                    base  = pread.alignment.query_sequence[pread.query_position].upper()
                    if rname not in result:
                        result[rname] = base
                break
    except (ValueError, KeyError):
        pass

    return result


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

        gene_chrom_rname   = gene_info["chrom"]
        pseudo_chrom_rname = pseudo_info["chrom"]
        gene_neg_strand   = gene_info.get("strand", "+") == "-"
        pseudo_neg_strand = pseudo_info.get("strand", "+") == "-"

        n_positions = 0
        n_reads     = 0

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
                gene_base_tsv   = row["Gene_Base"].upper()
                pseudo_base_tsv = row["Pseudogene_Base"].upper()

                gene_base_plus   = complement(gene_base_tsv)   if gene_neg_strand   else gene_base_tsv
                pseudo_base_plus = complement(pseudo_base_tsv) if pseudo_neg_strand else pseudo_base_tsv

                if gene_base_plus == pseudo_base_plus:
                    continue

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

                orig_bases   = bases_at(args.bam_original,  gene_chrom_rname,   gene_pos,
                                        args.min_base_quality, args.min_map_quality)
                gene_bases   = bases_at(gene_bam,            gene_chrom_rname,   gene_pos,
                                        args.min_base_quality, args.min_map_quality)
                pseudo_bases = bases_at(pseudo_bam,          pseudo_chrom_rname, pseudo_pos,
                                        args.min_base_quality, args.min_map_quality)

                all_reads = set(orig_bases) | set(gene_bases) | set(pseudo_bases)
                n_positions += 1
                n_reads     += len(all_reads)

                for rname in sorted(all_reads):
                    b_orig   = normalize_base(orig_bases.get(rname,   "."), gene_base_plus, pseudo_base_plus)
                    b_gene   = normalize_base(gene_bases.get(rname,   "."), gene_base_plus, pseudo_base_plus)
                    b_pseudo = normalize_base(pseudo_bases.get(rname, "."), gene_base_plus, pseudo_base_plus)
                    if b_gene != b_pseudo:
                        misassigned = "ambiguous"
                    elif b_orig == ".":
                        if b_gene == gene_base_plus:
                            misassigned = "gene"
                        elif b_gene == pseudo_base_plus:
                            misassigned = "pseudogene"
                        else:
                            misassigned = "ambiguous"
                    else:
                        if b_orig == b_gene:
                            if b_orig == gene_base_plus:
                                misassigned = "gene"
                            elif b_orig == pseudo_base_plus:
                                misassigned = "pseudogene"
                            else:
                                misassigned = "ambiguous"
                        else:
                            misassigned = "misassigned"
                    out_rows.append({
                        "sample":          args.sample,
                        "pair_id":         pair_id,
                        "read_name":       rname,
                        "gene_chrom":      gene_chrom,
                        "gene_pos":        gene_pos,
                        "gene_base":       gene_base_plus,
                        "pseudo_chrom":    pseudo_chrom,
                        "pseudo_pos":      pseudo_pos,
                        "pseudo_base":     pseudo_base_plus,
                        "base_original":   b_orig,
                        "base_gene_bam":   b_gene,
                        "base_pseudo_bam": b_pseudo,
                        "supports_gene":   b_gene   == gene_base_plus,
                        "supports_pseudo": b_pseudo == pseudo_base_plus,
                        "misassigned":     misassigned,
                    })

        print(
            f"[INFO] {pair_id}: {n_positions} mismatch positions, "
            f"{n_reads} total read×position observations"
        )

    fieldnames = [
        "sample", "pair_id", "read_name",
        "gene_chrom", "gene_pos", "gene_base",
        "pseudo_chrom", "pseudo_pos", "pseudo_base",
        "base_original", "base_gene_bam", "base_pseudo_bam",
        "supports_gene", "supports_pseudo", "misassigned",
    ]

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"[INFO] {args.sample}: wrote {len(out_rows)} read×position rows -> {args.output}")


if __name__ == "__main__":
    main()
