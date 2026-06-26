import argparse
import os
import re


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--xcnv",        required=True)
    p.add_argument("--sample-name", required=True)
    p.add_argument("--output",      required=True)
    p.add_argument("--log",         default="xhmm_extract_sample.log")
    return p.parse_args()


def _write_empty(path):
    with open(path, "w") as f:
        f.write("chromosome\tstart\tend\tCN\ttype\n")


def main():
    args = parse_args()
    xcnv_file   = args.xcnv
    tsv_out     = args.output
    sample_name = args.sample_name
    log         = args.log

    with open(log, "wt") as f:
        f.write("\n##\n## XHMM extract sample\n##\n")

    if (not os.path.exists(xcnv_file)) or os.path.getsize(xcnv_file) == 0:
        _write_empty(tsv_out)
        with open(log, "at") as f:
            f.write(f"XHMM xcnv file missing or empty: {xcnv_file}\n")
        return

    with open(xcnv_file) as in_f:
        lines = [line.rstrip("\n") for line in in_f if line.strip()]

    if not lines:
        _write_empty(tsv_out)
        return

    header = re.split(r"\s+", lines[0].strip())
    idx    = {name: i for i, name in enumerate(header)}

    required = ["SAMPLE", "CNV", "INTERVAL"]
    if any(col not in idx for col in required):
        _write_empty(tsv_out)
        with open(log, "at") as f:
            f.write("XHMM xcnv header missing required columns: SAMPLE/CNV/INTERVAL\n")
        return

    records = []
    for line in lines[1:]:
        fields = re.split(r"\s+", line.strip())
        if len(fields) <= max(idx["SAMPLE"], idx["CNV"], idx["INTERVAL"]):
            continue
        if fields[idx["SAMPLE"]] != sample_name:
            continue

        cnv_type = fields[idx["CNV"]].upper()
        interval = fields[idx["INTERVAL"]]

        m = re.match(r"^([^:]+):(\d+)-(\d+)$", interval)
        if not m:
            continue

        chrom = m.group(1)
        start = int(m.group(2))
        end   = int(m.group(3))
        if end < start:
            continue

        if cnv_type == "DEL":
            cn, out_type = "1", "deletion"
        elif cnv_type == "DUP":
            cn, out_type = "3", "duplication"
        else:
            cn, out_type = "", cnv_type.lower()

        records.append((chrom, start, end, cn, out_type))

    with open(tsv_out, "w") as out_f:
        out_f.write("chromosome\tstart\tend\tCN\ttype\n")
        for rec in records:
            out_f.write("%s\t%d\t%d\t%s\t%s\n" % rec)

    with open(log, "at") as f:
        f.write(f"Wrote {len(records)} XHMM CNV calls for sample {sample_name}\n")


if __name__ == "__main__":
    main()
