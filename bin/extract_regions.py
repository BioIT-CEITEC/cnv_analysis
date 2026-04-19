
#
# Usage: python 2_extract_regions.py --bam sample.bam --bed regions.bed --outdir extracted/

import pysam
import argparse
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
            if len(parts) < 5:
                raise ValueError(
                    f"The BED must have 5 columns (chrom, start, end, name, pair_id). Line: {line}"
                )
            chrom, start, end, name, pair_id = (
                parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4]
            )
            region_type = "pseudogene" if any(k in name.lower() for k in ["pseudogene", "pseudo", "_ps", "ps_"]) else "gene"
            pairs[pair_id].append({
                "chrom": chrom, "start": start, "end": end,
                "name": name, "type": region_type, "pair_id": pair_id
            })
    return pairs


def extract_by_pairs(bam_path, bed_path, outdir):
    """
    For each homologous pair, extract all reads that map to
    ANY of its regions (gene or pseudogene) into a single BAM.

    The origin is not labeled because:
    - Pseudogene reads can be in gene coordinates
    - Gene reads can be in pseudogene coordinates
    The true origin is determined by re-alignment in step 3/4.
    """
    os.makedirs(outdir, exist_ok=True)
    pairs = parse_bed(bed_path)

    with pysam.AlignmentFile(bam_path, "rb") as bam_in:
        for pair_id, regions in pairs.items():
            out_bam = os.path.join(outdir, f"{pair_id}.bam")
            tmp_bam = out_bam + ".tmp.bam"

            gene_regions   = [r for r in regions if r["type"] == "gene"]
            pseudo_regions = [r for r in regions if r["type"] == "pseudogene"]
            print(f"[2_extract] Par {pair_id}: "
                f"{len(gene_regions)} gene region(s) + "
                f"{len(pseudo_regions)} pseudogene region(s)")

            read_ids_written = set()
            total = 0

            with pysam.AlignmentFile(tmp_bam, "wb", header=bam_in.header) as bam_out:
                for region in regions:
                    for read in bam_in.fetch(
                        region["chrom"], region["start"], region["end"]
                    ):
                        if read.is_unmapped or read.query_name is None:
                            continue
                        if read.query_name not in read_ids_written:
                            bam_out.write(read)
                            read_ids_written.add(read.query_name)
                            total += 1

            pysam.sort("-o", out_bam, tmp_bam)
            pysam.index(out_bam)
            os.remove(tmp_bam)

            # Summary by region type according to mapping position
            gene_coords   = [(r["start"], r["end"]) for r in gene_regions]
            pseudo_coords = [(r["start"], r["end"]) for r in pseudo_regions]

            in_gene, in_pseudo, ambiguous = 0, 0, 0
            with pysam.AlignmentFile(out_bam, "rb") as bam_check:
                for read in bam_check.fetch():
                    pos = read.reference_start
                    in_g = any(s <= pos < e for s, e in gene_coords)
                    in_p = any(s <= pos < e for s, e in pseudo_coords)
                    if in_g and not in_p:
                        in_gene += 1
                    elif in_p and not in_g:
                        in_pseudo += 1
                    else:
                        ambiguous += 1

                print(f"[2_extract]   Total extracted: {total}")
                print(f"[2_extract]   Mapped in gene:      {in_gene}")
                print(f"[2_extract]   Mapped in pseudogene:{in_pseudo}")
                print(f"[2_extract]   Out/ambiguous:        {ambiguous}")
                print(f"[2_extract]   WARNING: pseudogene reads={in_pseudo} may be low")
                print(f"[2_extract]   if the aligner redirected reads from pseudogene to gene.")
                print(f"[2_extract]   Re-alignment (step 3) will determine the true origin.")
            print(f"[2_extract]   -> {out_bam}\n")

            print(f"[2_extract] Done. BAMs per pair in: {outdir}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract reads from gene+pseudogene regions by homologous pair"
    )
    parser.add_argument("--bam",    required=True, help="Original BAM aligned to the whole genome")
    parser.add_argument("--bed",    required=True, help="BED with 5 columns (chrom, start, end, name, pair_id)")
    parser.add_argument("--outdir", required=True, help="Output directory with one BAM per pair")
    args = parser.parse_args()

    extract_by_pairs(args.bam, args.bed, args.outdir)