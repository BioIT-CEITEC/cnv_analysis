#!/usr/bin/env python3


import pysam
import argparse
import subprocess
import os
import tempfile
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def _infer_region_type(name, extra_fields=()):
    """Detect gene vs pseudogene from name/annotations."""
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
    """Parse BED file
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
            strand = parts[5] if len(parts) > 5 else "+"
            extra  = parts[6:] if len(parts) > 6 else []
            rtype  = _infer_region_type(name, extra)
            region = {
                "chrom": chrom, "start": start, "end": end,
                "name": name, "strand": strand,
            }
            if rtype == "gene" and raw[pair_id]["gene"] is None:
                raw[pair_id]["gene"] = region
            elif rtype == "pseudogene" and raw[pair_id]["pseudogene"] is None:
                raw[pair_id]["pseudogene"] = region
            else:
                raw[pair_id]["unknown"].append(region)

    pairs = defaultdict(list)
    for pair_id, entry in raw.items():
        gene, pseudo, unknown = entry["gene"], entry["pseudogene"], list(entry["unknown"])
        if gene is None and unknown:
            gene = unknown.pop(0)
        if pseudo is None and unknown:
            pseudo = unknown.pop(0)
        if gene:
            pairs[pair_id].append(gene | {"type": "gene"})
        if pseudo:
            pairs[pair_id].append(pseudo | {"type": "pseudogene"})
    return pairs


def build_reference(ref_fa, regions, out_fasta, flank=300):
    """Build a FASTA for the target regions + flanks.
    """
    records = []
    coord_info = []
    with pysam.FastaFile(ref_fa) as fa:
        for region in regions:
            chrom_len = fa.get_reference_length(region["chrom"])
            start = max(0, region["start"] - flank)
            end   = min(chrom_len, region["end"] + flank)
            seq   = fa.fetch(region["chrom"], start, end)

            record = SeqRecord(
                Seq(seq),
                id=region["chrom"],
                description=(
                    f"region={region['chrom']}:{start}-{end} "
                    f"name={region['name']} type={region['type']} "
                    f"genomic_offset={start} flank={flank}"
                ),
            )
            records.append(record)
            coord_info.append({
                "chrom":  region["chrom"],
                "offset": start,
                "length": end - start,
                "strand": region.get("strand", "+"),
            })
    with open(out_fasta, "w") as f:
        SeqIO.write(records, f, "fasta")
    return coord_info


def liftover_bam_coords(in_bam, coord_info, genome_fai, out_bam):
    """Rewrite BAM coordinates from local (custom-ref) space to genomic space.

    For each read:
      - The contig name in the BAM header already matches the real chromosome
      - POS is local (0-based within the extracted region).
      - Add the genomic offset so POS becomes the true genomic position.
    """

    region_map = {}
    for ci in coord_info:
        if ci["chrom"] not in region_map:
            region_map[ci["chrom"]] = ci

    chrom_lengths = {}
    with open(genome_fai) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                chrom_lengths[parts[0]] = int(parts[1])

    with pysam.AlignmentFile(in_bam, "rb") as src:
        old_header = src.header.to_dict()

        new_sq = []
        for sq in old_header.get("SQ", []):
            name = sq["SN"]
            if name in chrom_lengths:
                sq = dict(sq)
                sq["LN"] = chrom_lengths[name]
            new_sq.append(sq)
        old_header["SQ"] = new_sq

        new_header = pysam.AlignmentHeader.from_dict(old_header)

        with pysam.AlignmentFile(out_bam, "wb", header=new_header) as dst:
            for read in src:
                if not read.is_unmapped:
                    ref_name = src.get_reference_name(read.reference_id)
                    ci       = region_map.get(ref_name, {"offset": 0, "strand": "+"})
                    offset   = ci["offset"]
                    neg      = ci["strand"] == "-"

                    read.reference_start += offset

                    if neg:
                        read.is_reverse       = not read.is_reverse
                        read.mate_is_reverse  = not read.mate_is_reverse
                        read.template_length  = -read.template_length

                    # Also fix mate coordinates if paired
                    if not read.mate_is_unmapped and read.next_reference_id >= 0:
                        mate_ref    = src.get_reference_name(read.next_reference_id)
                        mate_ci     = region_map.get(mate_ref, {"offset": 0})
                        read.next_reference_start += mate_ci["offset"]

                dst.write(read)


def count_reads_fastq(fq):
    if not os.path.exists(fq) or os.path.getsize(fq) == 0:
        return 0
    result = subprocess.run(["wc", "-l", fq], capture_output=True, text=True)
    return int(result.stdout.strip().split()[0]) // 4


def bam_to_fastq(bam_path, fq1, fq2, fq_singles):
    cmd = (
        f"samtools collate -u -O {bam_path} "
        f"| samtools fastq -1 {fq1} -2 {fq2} -s {fq_singles} -0 /dev/null -n"
    )
    subprocess.run(cmd, shell=True, check=True)
    r1 = count_reads_fastq(fq1)
    r2 = count_reads_fastq(fq2)
    s  = count_reads_fastq(fq_singles)
    print(f"[3_realign]   FASTQ: R1={r1}, R2={r2}, singletons={s}")
    return r1, r2, s


def align_to_ref(ref_fa, fq1, fq2, fq_singles, out_bam, threads=8):
    """Align reads against the custom reference and produce a sorted BAM."""
    subprocess.run(
        ["bwa", "index", ref_fa], check=True,
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )

    sam_paired  = out_bam.replace(".bam", ".paired.sam")
    sam_singles = out_bam.replace(".bam", ".singles.sam")

    with open(sam_paired, "w") as f:
        subprocess.run(
            ["bwa", "mem", "-t", str(threads), ref_fa, fq1, fq2],
            stdout=f, check=True, stderr=subprocess.DEVNULL,
        )

    singles_exist = os.path.exists(fq_singles) and os.path.getsize(fq_singles) > 0
    if singles_exist:
        with open(sam_singles, "w") as f:
            subprocess.run(
                ["bwa", "mem", "-t", str(threads), ref_fa, fq_singles],
                stdout=f, check=True, stderr=subprocess.DEVNULL,
            )

    bam_paired  = out_bam.replace(".bam", ".paired.bam")
    bam_singles = out_bam.replace(".bam", ".singles.bam")

    subprocess.run(["samtools", "sort", "-o", bam_paired, sam_paired], check=True)
    os.remove(sam_paired)

    if singles_exist:
        subprocess.run(["samtools", "sort", "-o", bam_singles, sam_singles], check=True)
        os.remove(sam_singles)
        subprocess.run(
            ["samtools", "merge", "-f", out_bam, bam_paired, bam_singles], check=True
        )
        os.remove(bam_paired)
        os.remove(bam_singles)
    else:
        os.rename(bam_paired, out_bam)


def realign_all_pairs(bam_dir, ref_fa, bed_path, outdir, threads=8, flank=300):
    """
    For each homologous pair:

    1. Build a local reference for the gene only
    2. Build a local reference for the pseudogene only
    3. Align the extracted reads against each custom reference
    4. Liftover BAM coordinates from local space to true genomic positions
    5. Sort, index, and write the final BAMs
    """
    os.makedirs(outdir, exist_ok=True)

    genome_fai = ref_fa + ".fai"
    if not os.path.exists(genome_fai):
        raise FileNotFoundError(
            f"Genome .fai not found: {genome_fai}\n"
            f"Run: samtools faidx {ref_fa}"
        )

    pairs = parse_bed(bed_path)

    for pair_id, regions in pairs.items():
        bam_path = os.path.join(bam_dir, f"{pair_id}.bam")
        if not os.path.exists(bam_path):
            print(f"[WARNING] {bam_path} not found, skipping {pair_id}")
            continue

        gene_regions   = [r for r in regions if r["type"] == "gene"]
        pseudo_regions = [r for r in regions if r["type"] == "pseudogene"]

        if not gene_regions or not pseudo_regions:
            print(f"[WARNING] {pair_id} does not have both gene and pseudogene defined")
            continue

        out_gene   = os.path.join(outdir, f"{pair_id}_gene.bam")
        out_pseudo = os.path.join(outdir, f"{pair_id}_pseudogene.bam")
        print(f"\n[3_realign] Processing pair {pair_id}...")
        print(f"  Gene:       {[(r['chrom'], r['start'], r['end']) for r in gene_regions]}")
        print(f"  Pseudogene: {[(r['chrom'], r['start'], r['end']) for r in pseudo_regions]}")

        with tempfile.TemporaryDirectory() as tmpdir:
            ref_gene   = os.path.join(tmpdir, "ref_gene.fa")
            ref_pseudo = os.path.join(tmpdir, "ref_pseudo.fa")
            fq1        = os.path.join(tmpdir, "R1.fastq")
            fq2        = os.path.join(tmpdir, "R2.fastq")
            fq_singles = os.path.join(tmpdir, "singles.fastq")

            gene_coord_info   = build_reference(ref_fa, gene_regions,   ref_gene,   flank=flank)
            pseudo_coord_info = build_reference(ref_fa, pseudo_regions, ref_pseudo, flank=flank)

            r1, r2, s = bam_to_fastq(bam_path, fq1, fq2, fq_singles)
            if r1 == 0 and r2 == 0 and s == 0:
                print(f"[WARNING] Empty FASTQ for {pair_id}, skipping")
                continue

            # Align to custom refs (produces local coordinates)
            bam_gene_local   = os.path.join(tmpdir, "gene_local.bam")
            bam_pseudo_local = os.path.join(tmpdir, "pseudo_local.bam")

            print(f"[3_realign] Aligning against gene reference...")
            align_to_ref(ref_gene,   fq1, fq2, fq_singles, bam_gene_local,   threads)

            print(f"[3_realign] Aligning against pseudogene reference...")
            align_to_ref(ref_pseudo, fq1, fq2, fq_singles, bam_pseudo_local, threads)

            # Liftover local coords → true genomic coords
            bam_gene_lifted   = os.path.join(tmpdir, "gene_lifted.bam")
            bam_pseudo_lifted = os.path.join(tmpdir, "pseudo_lifted.bam")

            print(f"[3_realign] Lifting over gene BAM coordinates...")
            liftover_bam_coords(bam_gene_local,   gene_coord_info,   genome_fai, bam_gene_lifted)

            print(f"[3_realign] Lifting over pseudogene BAM coordinates...")
            liftover_bam_coords(bam_pseudo_local, pseudo_coord_info, genome_fai, bam_pseudo_lifted)

            # Final sort + index (outside tmpdir so they survive)
            subprocess.run(["samtools", "sort", "-o", out_gene,   bam_gene_lifted],   check=True)
            subprocess.run(["samtools", "sort", "-o", out_pseudo, bam_pseudo_lifted], check=True)

        subprocess.run(["samtools", "index", out_gene],   check=True)
        subprocess.run(["samtools", "index", out_pseudo], check=True)

        with pysam.AlignmentFile(out_gene, "rb") as b:
            mapped_gene = b.mapped
        with pysam.AlignmentFile(out_pseudo, "rb") as b:
            mapped_pseudo = b.mapped

        print(f"[3_realign] Reads in gene BAM:       {mapped_gene}")
        print(f"[3_realign] Reads in pseudogene BAM: {mapped_pseudo}")

    print(f"\n[3_realign] Done. BAMs in: {outdir}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Re-align reads against specific gene/pseudogene references "
                    "and liftover coordinates to true genomic positions for IGV."
    )
    parser.add_argument("--bam_dir",  required=True)
    parser.add_argument("--ref",      required=True, help="Full genome FASTA (must have .fai)")
    parser.add_argument("--bed",      required=True)
    parser.add_argument("--outdir",   required=True)
    parser.add_argument("--threads",  type=int, default=8)
    parser.add_argument("--flank",    type=int, default=300)
    args = parser.parse_args()

    realign_all_pairs(
        args.bam_dir, args.ref, args.bed, args.outdir,
        threads=args.threads, flank=args.flank,
    )