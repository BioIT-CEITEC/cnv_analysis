import argparse
import os
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bed",    required=True)
    p.add_argument("--gtf",    required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--log",    default="conifer_prepare_probes.log")
    return p.parse_args()


def write_log(log_filename, message):
    with open(log_filename, "at") as f:
        f.write(message + "\n")


def normalise_chromosome(chrom):
    chrom = chrom.strip()
    return chrom[3:] if chrom.startswith("chr") else chrom


def parse_bed(path):
    intervals = []
    with open(path) as f:
        for line_number, raw_line in enumerate(f, start=1):
            line = raw_line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError("BED line %d has fewer than 3 columns" % line_number)
            chrom = fields[0].strip()
            start = int(fields[1])
            stop  = int(fields[2])
            if stop <= start:
                raise ValueError("BED line %d has invalid interval" % line_number)
            intervals.append({
                "chrom":      chrom,
                "start":      start,
                "stop":       stop,
                "chrom_norm": normalise_chromosome(chrom),
            })
    if not intervals:
        raise ValueError("No usable intervals found in %s" % path)
    return intervals


def parse_gtf_attributes(raw_attributes):
    attributes = {}
    for item in raw_attributes.strip().split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        attributes[key] = value.strip().strip('"')
    return attributes


def load_gtf_records(path, feature_name):
    records = defaultdict(list)
    with open(path) as f:
        for raw_line in f:
            if not raw_line or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != feature_name:
                continue
            attrs     = parse_gtf_attributes(fields[8])
            gene_name = attrs.get("gene_name") or attrs.get("gene_id")
            if not gene_name:
                continue
            gene_type = attrs.get("gene_type") or attrs.get("gene_biotype")
            if gene_type and gene_type != "protein_coding":
                continue
            chrom_norm = normalise_chromosome(fields[0])
            start      = int(fields[3]) - 1
            stop       = int(fields[4])
            records[chrom_norm].append((start, stop, gene_name))
    for chrom in records:
        records[chrom].sort(key=lambda item: (item[0], item[1], item[2]))
    return records


def load_gtf(path, log_filename):
    records = load_gtf_records(path, "gene")
    if records:
        return records
    write_log(log_filename, "No gene features found in GTF; falling back to exon features.")
    return load_gtf_records(path, "exon")


def annotate_intervals(intervals, gtf_records):
    by_chrom = defaultdict(list)
    for idx, interval in enumerate(intervals):
        by_chrom[interval["chrom_norm"]].append((idx, interval))

    annotations = ["." for _ in intervals]
    for chrom_norm, indexed_intervals in by_chrom.items():
        genes = gtf_records.get(chrom_norm, [])
        if not genes:
            continue
        indexed_intervals.sort(key=lambda item: (item[1]["start"], item[1]["stop"]))
        left = 0
        for interval_idx, interval in indexed_intervals:
            while left < len(genes) and genes[left][1] <= interval["start"]:
                left += 1
            j, seen, names = left, set(), []
            while j < len(genes) and genes[j][0] < interval["stop"]:
                g_start, g_stop, g_name = genes[j]
                if g_stop > interval["start"] and g_name not in seen:
                    seen.add(g_name)
                    names.append(g_name)
                j += 1
            if names:
                annotations[interval_idx] = ",".join(names)
    return annotations


def main():
    args         = parse_args()
    log_filename = args.log

    with open(log_filename, "wt") as f:
        f.write("\n##\n## conifer_prepare_probes\n##\n")

    write_log(log_filename, "Input BED  : %s" % args.bed)
    write_log(log_filename, "Input GTF  : %s" % args.gtf)
    write_log(log_filename, "Output TSV : %s" % args.output)

    intervals   = parse_bed(args.bed)
    gtf_records = load_gtf(args.gtf, log_filename)
    genes       = annotate_intervals(intervals, gtf_records)

    with open(args.output, "wt") as f:
        for interval, gene_name in zip(intervals, genes):
            f.write("%s\t%d\t%d\t%s\n" % (
                interval["chrom"], interval["start"], interval["stop"], gene_name,
            ))

    write_log(log_filename, "Prepared %d CONIFER probes." % len(intervals))


if __name__ == "__main__":
    main()
