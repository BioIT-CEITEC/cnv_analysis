#!/usr/bin/env python3
"""
Annotate a region_coverage.tsv with GTF-derived strand, gene_biotype, and gene_label.

Columns appended: strand | gene_biotype | gene_label

For each gene, a single canonical transcript is selected (MANE Select if present,
otherwise the transcript with the most exons) and used consistently across all
capture regions of that gene.  The transcript ID is embedded in the label so the
source of exon numbering is unambiguous.

gene_label examples:
  BRCA1[ENST00000357654.9]:exon5          single exon overlap
  BRCA1[ENST00000357654.9]:exon3,4        region spans two exons
  BRCA1[ENST00000357654.9]               gene overlaps but no exon (intronic/UTR)
  VPS13B[ENST00000261509.10]:exon3;COH1[ENST00000395306.5]:exon5   multi-gene region

BED col 4 may be comma-separated for multi-gene capture regions ("VPS13B,COH1").
Name matching is exact OR substring (BED name contained in GTF name) to handle
panel aliases (e.g. BED "AIP" matching GTF "AURKAIP1").

Usage:
    annotate_coverage_with_gtf.py <gtf_file> <region_coverage.tsv> <out.tsv>
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
    Return True when the GTF gene is considered a match for the BED name set.
    Accepts exact match or BED name as substring of GTF name (BED aliases are
    often shorter, e.g. BED "AIP" matches GTF "AURKAIP1").
    Coordinate overlap is the primary filter so substring matching is safe.
    """
    if not bed_names:
        return True
    for g in bed_names:
        if g == gtf_gene or g in gtf_gene:
            return True
    return False


def select_canonical_transcripts(gtf_path):
    """
    Two-pass GTF read that returns:
      gene_canonical : {gene_name -> transcript_id}
      gene_by_chr   : {chrom -> [(start0, end1, strand, gene_name, gene_biotype)]}

    Canonical transcript priority:
      1. MANE_Select tag on the transcript record
      2. Transcript with the greatest number of exons
      3. First encountered (stable sort tiebreaker)
    """
    # Pass 1 — collect transcript metadata and per-transcript exon counts
    tx_meta   = {}          # transcript_id -> {gene_name, is_mane}
    tx_exons  = defaultdict(int)   # transcript_id -> exon count
    gene_by_chr = defaultdict(list)

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
                tx_id   = get_attr(attrs, "transcript_id")
                is_mane = "MANE_Select" in attrs
                tx_meta[tx_id] = {"gene_name": gene_name, "is_mane": is_mane}

            elif feature == "exon":
                tx_id = get_attr(attrs, "transcript_id")
                tx_exons[tx_id] += 1
                # If transcript records are absent, register the transcript here
                if tx_id not in tx_meta:
                    tx_meta[tx_id] = {"gene_name": gene_name, "is_mane": False}

    # Group transcripts by gene and select canonical
    gene_tx_lists = defaultdict(list)
    for tx_id, meta in tx_meta.items():
        gene_tx_lists[meta["gene_name"]].append(
            (meta["is_mane"], tx_exons[tx_id], tx_id)
        )

    gene_canonical = {}
    for gene_name, tx_list in gene_tx_lists.items():
        # Sort: MANE_Select first (-is_mane), then most exons (-count), then ID (stable)
        tx_list.sort(key=lambda x: (-x[0], -x[1], x[2]))
        gene_canonical[gene_name] = tx_list[0][2]

    return gene_canonical, gene_by_chr


def build_exon_index(gtf_path, gene_canonical):
    """
    Return exon_by_chr indexing only the canonical transcript per gene.
    Each entry: (start0, end1, gene_name, transcript_id, exon_number)
    """
    exon_by_chr = defaultdict(list)
    seen        = set()   # (chrom, start0, end1, gene_name) — already rare with canonical tx

    with open_file(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "exon":
                continue
            attrs   = p[8]
            tx_id   = get_attr(attrs, "transcript_id")
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

            exon_number = get_attr(attrs, "exon_number")
            exon_by_chr[chrom].append((start0, end1, gene_name, tx_id, exon_number))

    return exon_by_chr


def main():
    gtf_path      = sys.argv[1]
    coverage_path = sys.argv[2]
    out_path      = sys.argv[3]

    gene_canonical, gene_by_chr = select_canonical_transcripts(gtf_path)
    exon_by_chr = build_exon_index(gtf_path, gene_canonical)

    with open(coverage_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            p = line.rstrip("\n").split("\t")
            chrom  = norm_chrom(p[0])
            start0 = int(p[1])
            end1   = int(p[2])

            # BED col 4 may have comma-separated gene names (e.g. "VPS13B,COH1")
            bed_names = set()
            if len(p) > 3 and p[3]:
                bed_names = {g.strip() for g in p[3].split(",") if g.strip()}

            def collect_gene_annots(name_filter):
                """
                Build per-gene annotation dict for genes overlapping this region.
                When name_filter is True, only include genes whose GTF name matches
                the BED names (exact or substring).  When False, include all overlapping
                genes regardless of name.
                """
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

            # Coordinate-only lookup → used for gtf_gene_name (always shows what is
            # physically at these coordinates, regardless of BED naming convention)
            all_annots = collect_gene_annots(name_filter=False)

            # Name-filtered lookup → used for strand / biotype / gene_label
            # Falls back to all overlapping genes when no name matches, so that
            # regions belonging to renamed genes are still fully annotated
            gene_annots = collect_gene_annots(name_filter=True)
            if not gene_annots:
                gene_annots = all_annots

            # GTF gene names: all coordinate-overlapping genes, NA only if truly none
            gtf_names_str = ",".join(sorted(all_annots.keys())) if all_annots else "NA"

            strands  = {v["strand"]  for v in gene_annots.values() if v["strand"]  != "NA"}
            biotypes = {v["biotype"] for v in gene_annots.values()
                        if v["biotype"] and v["biotype"] != "NA"}

            strand_str  = ";".join(sorted(strands))  if strands  else "NA"
            biotype_str = ";".join(sorted(biotypes)) if biotypes else "NA"

            # Build gene_label: GENENAME[TRANSCRIPT]:exon3,4  or  GENENAME[TRANSCRIPT]
            label_parts = []
            for gname in sorted(gene_annots.keys()):
                info = gene_annots[gname]
                tx   = info["transcript"]
                tx_tag = "[{}]".format(tx) if tx else ""

                exons = info["exons"]
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

            fout.write(
                "\t".join(p) + "\t" + strand_str + "\t" + biotype_str
                + "\t" + gtf_names_str + "\t" + gene_label + "\n"
            )


if __name__ == "__main__":
    main()
