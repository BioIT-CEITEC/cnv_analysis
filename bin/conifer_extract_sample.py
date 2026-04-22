import argparse
import os


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--calls",       required=True)
    p.add_argument("--sample-name", required=True)
    p.add_argument("--output",      required=True)
    p.add_argument("--log",         default="conifer_extract_sample.log")
    return p.parse_args()


def main():
    args        = parse_args()
    sample_name = args.sample_name
    log         = args.log

    with open(log, "wt") as f:
        f.write("\n##\n## conifer_extract_sample\n##\n")

    candidates = {sample_name, sample_name + ".rpkm"}
    if sample_name.endswith(".rpkm"):
        candidates.add(sample_name[:-5])

    rows_out = []
    if os.path.exists(args.calls):
        with open(args.calls) as fin:
            for line in fin:
                line = line.strip()
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) < 5:
                    continue
                if fields[0] in ("sampleID", "SampleID"):
                    continue
                if fields[0] not in candidates:
                    continue

                state = fields[4].strip().lower()
                if state in ("del", "deletion"):
                    cnv_type = "deletion"
                elif state in ("dup", "duplication"):
                    cnv_type = "duplication"
                else:
                    continue

                rows_out.append((fields[1], fields[2], fields[3], "", cnv_type))

    with open(args.output, "w") as fout:
        fout.write("chromosome\tstart\tend\tCN\ttype\n")
        for row in rows_out:
            fout.write("%s\t%s\t%s\t%s\t%s\n" % row)

    with open(log, "at") as f:
        f.write("sample: %s\n" % sample_name)
        f.write("input calls: %s\n" % args.calls)
        f.write("calls written: %d\n" % len(rows_out))


if __name__ == "__main__":
    main()
