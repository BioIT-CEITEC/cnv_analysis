#!/usr/bin/env python3
#
# Usage: python 3_realign_specific.py --bam_dir extracted/ --ref genome.fa
#                                    --bed regions.bed --outdir realigned/

import pysam
import argparse
import subprocess
import os
import tempfile
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


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
                "name": name, "type": region_type
            })
    return pairs


def build_reference(ref_fa, regions, out_fasta, flank=300):
    """Build a FASTA with the specified regions + flanks."""
    records = []
    with pysam.FastaFile(ref_fa) as fa:
        for region in regions:
            chrom_len = fa.get_reference_length(region["chrom"])
            start = max(0, region["start"] - flank)
            end   = min(chrom_len, region["end"] + flank)
            seq   = fa.fetch(region["chrom"], start, end)
            record = SeqRecord(
                Seq(seq),
                id=region["name"],
                description=f"{region['chrom']}:{start}-{end} type={region['type']} flank={flank}"
            )
            records.append(record)
    with open(out_fasta, "w") as f:
        SeqIO.write(records, f, "fasta")


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
    """
    Align the reads againsta specific reference (gene OR pseudogene) and produce 
    a BAM with the best alignment of each read against that reference.
    """
    subprocess.run(
        ["bwa", "index", ref_fa], check=True,
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )

    sam_paired  = out_bam.replace(".bam", ".paired.sam")
    sam_singles = out_bam.replace(".bam", ".singles.sam")

    with open(sam_paired, "w") as f:
        subprocess.run(
            ["bwa", "mem", "-t", str(threads), ref_fa, fq1, fq2],
            stdout=f, check=True, stderr=subprocess.DEVNULL
        )

    singles_exist = os.path.exists(fq_singles) and os.path.getsize(fq_singles) > 0
    if singles_exist:
        with open(sam_singles, "w") as f:
            subprocess.run(
                ["bwa", "mem", "-t", str(threads), ref_fa, fq_singles],
                stdout=f, check=True, stderr=subprocess.DEVNULL
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

    subprocess.run(["samtools", "index", out_bam], check=True)


def realign_all_pairs(bam_dir, ref_fa, bed_path, outdir, threads=8, flank=300):
    """
    For each homologous pair:

    1. Build a local reference ONLY for the gene
    2. Build a local reference ONLY for the pseudogene
    3. Align the extracted reads against each one separately
    4. Produce: pair1_gene.bam and pair1_pseudogene.bam

    Step 4 compares the AS between both BAMs to classify each read.
    """
    os.makedirs(outdir, exist_ok=True)
    pairs = parse_bed(bed_path)

    for pair_id, regions in pairs.items():
        bam_path = os.path.join(bam_dir, f"{pair_id}.bam")
        if not os.path.exists(bam_path):
            print(f"WARNING: {bam_path} not found, skipping {pair_id}")
            continue

        gene_regions   = [r for r in regions if r["type"] == "gene"]
        pseudo_regions = [r for r in regions if r["type"] == "pseudogene"]

        if not gene_regions or not pseudo_regions:
            print(f"[WARNING: {pair_id} does not have both gene and pseudogene defined")
            continue

        out_gene   = os.path.join(outdir, f"{pair_id}_gene.bam")
        out_pseudo = os.path.join(outdir, f"{pair_id}_pseudogene.bam")
        print(f"\n[3_realign] Processing pair {pair_id}...")

        with tempfile.TemporaryDirectory() as tmpdir:
            ref_gene   = os.path.join(tmpdir, "ref_gene.fa")
            ref_pseudo = os.path.join(tmpdir, "ref_pseudo.fa")
            fq1        = os.path.join(tmpdir, "R1.fastq")
            fq2        = os.path.join(tmpdir, "R2.fastq")
            fq_singles = os.path.join(tmpdir, "singles.fastq")

            build_reference(ref_fa, gene_regions,   ref_gene,   flank=flank)
            build_reference(ref_fa, pseudo_regions, ref_pseudo, flank=flank)
            print(f"Gene ref:   {[r['name'] for r in gene_regions]}")
            print(f"Pseudo ref: {[r['name'] for r in pseudo_regions]}")

            r1, r2, s = bam_to_fastq(bam_path, fq1, fq2, fq_singles)
            if r1 == 0 and r2 == 0 and s == 0:
                print(f"WARNING: Empty FASTQ for {pair_id}, skipping")
                continue

            print(f"Aligning against gene...")
            align_to_ref(ref_gene, fq1, fq2, fq_singles, out_gene, threads)

            print(f"Aligning against pseudogene...")
            align_to_ref(ref_pseudo, fq1, fq2, fq_singles, out_pseudo, threads)

        with pysam.AlignmentFile(out_gene, "rb") as b:
            mapped_gene = b.mapped
        with pysam.AlignmentFile(out_pseudo, "rb") as b:
            mapped_pseudo = b.mapped

        print(f"Reads aligned to gene:      {mapped_gene}")
        print(f"Reads aligned to pseudogene: {mapped_pseudo}")

    print(f"\n Done. BAMs in: {outdir}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Re-align each read separately against gene and pseudogene"
    )
    parser.add_argument("--bam_dir",  required=True)
    parser.add_argument("--ref",      required=True)
    parser.add_argument("--bed",      required=True)
    parser.add_argument("--outdir",   required=True)
    parser.add_argument("--threads",  type=int, default=8)
    parser.add_argument("--flank",    type=int, default=300)
    args = parser.parse_args()

    realign_all_pairs(
        args.bam_dir, args.ref, args.bed, args.outdir,
        threads=args.threads, flank=args.flank
    )
