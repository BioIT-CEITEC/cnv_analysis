#!/usr/bin/env python3
"""
Build a per-region GTF annotation file for a BED panel (run once per pipeline).

Reads the GTF and BED file, selects a canonical transcript per gene (MANE Select
if available, otherwise most exons), and writes one annotation row per BED region.

Output columns (tab-separated, no header, same row order as the BED):
  strand | gene_biotype | gtf_gene_name | gene_label

gene_label format:
  GENENAME[TRANSCRIPT_ID]:exon3      exon overlap
  GENENAME[TRANSCRIPT_ID]            gene overlaps but no exon (intronic/UTR)
  GENENAME[TRANSCRIPT_ID]:exon3,4    region spans multiple exons

BED col 4 may be comma-separated for multi-gene regions ("VPS13B,COH1").
Name matching is exact OR BED name as substring of GTF name (handles aliases
such as BED "AIP" matching GTF "AURKAIP1").
When no name matches, all coordinate-overlapping genes are reported so renamed
genes are never silently dropped.

Usage:
    build_bed_gtf_annotation.py <gtf_file> <bed_file> <out.tsv>
"""

import sys
import re
import gzip
from collections import defaultdict


def open_file(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def get_attr(attrs, key):
    m = re.search(r'{} "([^"]+)"'.format(re.escape(key)), attrs)
    return m.group(1) if m else ""


def norm_chrom(c):
    """Strip 'chr' prefix so BED and GTF chromosome names always match."""
    return c[3:] if c.startswith("chr") else c


def name_matches(bed_names, gtf_gene):
    """
    Exact match OR BED name is a substring of the GTF name.
    Coordinate overlap is the primary filter so substring matching is safe.
    """
    if not bed_names:
        return True
    for g in bed_names:
        if g == gtf_gene or g in gtf_gene:
            return True
    return False


def build_indexes(gtf_path):
    """
    Two-pass GTF read.  Returns:
      gene_canonical : {gene_name -> transcript_id}
      gene_by_chr    : {chrom -> [(start0, end1, strand, gene_name, gene_biotype)]}
      exon_by_chr    : {chrom -> [(start0, end1, gene_name, transcript_id, exon_number)]}

    Canonical transcript priority:
      1. MANE_Select tag
      2. Most exons
      3. First encountered (stable tiebreaker)
    """
    tx_meta  = {}               # transcript_id -> {gene_name, is_mane}
    tx_exons = defaultdict(int) # transcript_id -> exon count
    gene_by_chr = defaultdict(list)

    # Pass 1 — transcript metadata and exon counts
    with open_file(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            feature = p[2]
            if feature not in ("gene", "transcript", "exon"):
                continue

            chrom  = norm_chrom(p[0])
            start0 = int(p[3]) - 1
            end1   = int(p[4])
            strand = p[6]
            attrs  = p[8]
            gene_name = get_attr(attrs, "gene_name") or get_attr(attrs, "gene_id")

            if feature == "gene":
                gene_biotype = (get_attr(attrs, "gene_biotype")
                                or get_attr(attrs, "gene_type"))
                gene_by_chr[chrom].append(
                    (start0, end1, strand, gene_name, gene_biotype)
                )
            elif feature == "transcript":
                tx_id = get_attr(attrs, "transcript_id")
                tx_meta[tx_id] = {
                    "gene_name": gene_name,
                    "is_mane":   "MANE_Select" in attrs,
                }
            elif feature == "exon":
                tx_id = get_attr(attrs, "transcript_id")
                tx_exons[tx_id] += 1
                if tx_id not in tx_meta:
                    tx_meta[tx_id] = {"gene_name": gene_name, "is_mane": False}

    # Select canonical transcript per gene
    gene_tx_lists = defaultdict(list)
    for tx_id, meta in tx_meta.items():
        gene_tx_lists[meta["gene_name"]].append(
            (meta["is_mane"], tx_exons[tx_id], tx_id)
        )

    gene_canonical = {}
    for gene_name, tx_list in gene_tx_lists.items():
        tx_list.sort(key=lambda x: (-x[0], -x[1], x[2]))
        gene_canonical[gene_name] = tx_list[0][2]

    # Pass 2 — index only canonical-transcript exons
    exon_by_chr = defaultdict(list)
    seen = set()

    with open_file(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "exon":
                continue
            attrs     = p[8]
            tx_id     = get_attr(attrs, "transcript_id")
            gene_name = get_attr(attrs, "gene_name") or get_attr(attrs, "gene_id")

            if gene_canonical.get(gene_name) != tx_id:
                continue

            chrom  = norm_chrom(p[0])
            start0 = int(p[3]) - 1
            end1   = int(p[4])
            key    = (chrom, start0, end1, gene_name)
            if key in seen:
                continue
            seen.add(key)

            exon_by_chr[chrom].append(
                (start0, end1, gene_name, tx_id, get_attr(attrs, "exon_number"))
            )

    return gene_canonical, gene_by_chr, exon_by_chr


def annotate_region(chrom, start0, end1, bed_names,
                    gene_canonical, gene_by_chr, exon_by_chr):
    """Return (strand_str, biotype_str, gtf_names_str, gene_label) for one region."""

    def collect(name_filter):
        annots = {}
        for e_start, e_end, e_strand, e_gene, e_biotype in gene_by_chr.get(chrom, []):
            if e_end <= start0 or e_start >= end1:
                continue
            if name_filter and not name_matches(bed_names, e_gene):
                continue
            annots[e_gene] = {
                "strand":     e_strand,
                "biotype":    e_biotype,
                "transcript": gene_canonical.get(e_gene, ""),
                "exons":      set(),
            }
        for e_start, e_end, e_gene, e_tx, e_exon in exon_by_chr.get(chrom, []):
            if e_end <= start0 or e_start >= end1:
                continue
            if name_filter and not name_matches(bed_names, e_gene):
                continue
            if e_gene not in annots:
                annots[e_gene] = {
                    "strand":     "NA",
                    "biotype":    "NA",
                    "transcript": e_tx,
                    "exons":      set(),
                }
            if e_exon:
                annots[e_gene]["exons"].add(e_exon)
        return annots

    all_annots  = collect(name_filter=False)
    gene_annots = collect(name_filter=True) or all_annots

    gtf_names_str = ",".join(sorted(all_annots.keys())) if all_annots else "NA"

    # strand from all_annots so it reflects every gene overlapping the region,
    # consistent with gtf_gene_name (e.g. FANCC on - and AOPEP on + both shown)
    strands  = {v["strand"]  for v in all_annots.values()  if v["strand"]  != "NA"}
    # biotype from gene_annots (targeted gene only — what the capture was designed for)
    biotypes = {v["biotype"] for v in gene_annots.values()
                if v["biotype"] and v["biotype"] != "NA"}

    strand_str  = ";".join(sorted(strands))  if strands  else "NA"
    biotype_str = ";".join(sorted(biotypes)) if biotypes else "NA"

    label_parts = []
    for gname in sorted(gene_annots.keys()):
        info  = gene_annots[gname]
        tx    = info["transcript"]
        tx_tag = "[{}]".format(tx) if tx else ""
        exons  = info["exons"]
        if exons:
            try:
                exon_nums = sorted(exons, key=int)
            except ValueError:
                exon_nums = sorted(exons)
            label_parts.append(
                "{}{}:exon{}".format(gname, tx_tag, ",".join(str(x) for x in exon_nums))
            )
        else:
            label_parts.append("{}{}".format(gname, tx_tag))

    gene_label = ";".join(label_parts) if label_parts else "NA"

    return strand_str, biotype_str, gtf_names_str, gene_label


def main():
    gtf_path = sys.argv[1]
    bed_path = sys.argv[2]
    out_path = sys.argv[3]

    gene_canonical, gene_by_chr, exon_by_chr = build_indexes(gtf_path)

    with open(bed_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            chrom  = norm_chrom(p[0])
            start0 = int(p[1])
            end1   = int(p[2])

            bed_names = set()
            if len(p) > 3 and p[3]:
                bed_names = {g.strip() for g in p[3].split(",") if g.strip()}

            strand_str, biotype_str, gtf_names_str, gene_label = annotate_region(
                chrom, start0, end1, bed_names,
                gene_canonical, gene_by_chr, exon_by_chr
            )

            # Prefix with chr/start/end so the coverage process can join by coordinates
            # rather than relying on row-order alignment (paste).
            fout.write(
                p[0] + "\t" + p[1] + "\t" + p[2] + "\t"
                + strand_str + "\t" + biotype_str + "\t"
                + gtf_names_str + "\t" + gene_label + "\n"
            )


if __name__ == "__main__":
    main()
