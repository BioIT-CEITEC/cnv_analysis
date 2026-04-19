#!/usr/bin/env python

import os
import sys
import re
import math

# INPUT / OUTPUT
INPUT_DIR = sys.argv[1]
OUTPUT_DIR = sys.argv[2]

# raw_package_input = sys.argv[3] if len(sys.argv) > 3 else ""
# if raw_package_input and raw_package_input.lower() != "null":
#     PACKAGE_LIST = [p.strip().lower() for p in raw_package_input.split(",") if p.strip()]
# else:
#     PACKAGE_LIST = []
# 

os.makedirs(OUTPUT_DIR, exist_ok=True)

# CONSTANTS
KEEP_COLUMNS = ["CHROM", "START", "END", "CN_ESTIMATE", "TYPE"]

MISSING_HEADER = "CHROM\tSTART\tEND\tCOUNT\tTYPE"

MIN_KEYWORDS_REQUIRED = 2

HEADER_KEYWORDS = [
    "chromosome","start","end","length","id","type","nexons","BF",
    "reads.expected","reads.observed","reads.ratio","Gene","gene","Exon","Sample",
    "RC","medRC","RC.norm","medRC.norm","lowQual","CN","log2","cn","depth",
    "probes","weight","sample","region_id","chr","cn_pred","cn_id","cov",
    "norm_dist_mean","norm_dist_sd","pos","pop_HET_probability","alt_count","ref_count",
    "LOH","CNV_in_cohort","samples_in_cohort", "LINEAR_COPY_RATIO",
    "CHROM","POS","ID","REF","ALT","QUAL","FILTER","FOLD_CHANGE",
    "FOLD_CHANGE_LOG","CONTIG","IMPRECISE","SVLEN","SVTYPE",
    "WIDTH","STRAND","SAMPLENAME","MEDIAN","MEAN","COUNT"
]

REGEX = re.compile(
    r'^(GL\d+\.\d+|KI\d+\.\d+|Un_?GL\d+\.\d+|chrUn_GL\d+\.\d+|hs37d5|'
    r'chr[0-9XYM]+_.*_alt|GL\d+v\d+_alt)$'
    r'|SN:(GL\d+\.\d+|KI\d+\.\d+|Un_?GL\d+\.\d+|chrUn_GL\d+\.\d+|'
    r'hs37d5|chr[0-9XYM]+_.*_alt|GL\d+v\d+_alt)'
)

# HELPERS
def has_header(lines):
    for line in lines[:3]:
        line_lower = line.strip().lower()
        keyword_count = sum(1 for kw in HEADER_KEYWORDS if kw.lower() in line_lower)
        if keyword_count >= MIN_KEYWORDS_REQUIRED:
            return True
    return False


def normalize_header(header_line):
    cols = header_line.rstrip("\r\n").split("\t")
    # Pre-scan: only use POS as a START fallback when no explicit start column exists
    cols_upper = [c.strip().upper() for c in cols]
    has_explicit_start = "START" in cols_upper
    normalized = []
    for col in cols:
        col_clean = col.strip()
        col_upper = col_clean.upper()

        if col_upper in ["CHROMOSOME", "CONTIG", "CHR"]:
            normalized.append("CHROM")
        elif col_upper == "START" or (col_upper == "POS" and not has_explicit_start):
            normalized.append("START")
        elif col_upper == "END":
            normalized.append("END")
        elif col_upper in ["CN", "CN_PRED", "CN_PREDICTION", "CN_EST", "LINEAR_COPY_RATIO", "READS.RATIO", "FOLD_CHANGE"]:
            normalized.append("CN_ESTIMATE")
        elif col_upper in ["TYPE", "SVTYPE"]:
            normalized.append("TYPE")
        elif col_upper in ["SAMPLE", "SAMPLENAME", "SAMPLE_NAME"]:
            normalized.append("SAMPLE")
        else:
            normalized.append(col_clean)

    return "\t".join(normalized)


def round_half_up(value):
    v = value.strip()
    if v.upper().startswith("CN"):
        v = v[2:]
    return int(math.floor(float(v) + 0.5))


def normalize_type(value):
    v = value.strip().lower()
    if v in ["gain", "dup", "duplication"]:
        return "DUP"
    if v in ["loss", "del", "deletion"]:
        return "DEL"
    return value.strip()


def infer_type_from_cn(cn):
    if cn < 2:
        return "DEL"
    elif cn > 2:
        return "DUP"
    return "NEUTRAL"


# MAIN PROCESSING
abs_output_dir = os.path.abspath(OUTPUT_DIR)

for root, dirs, files in os.walk(INPUT_DIR, topdown=True):
    # Never recurse into the output directory — prevents re-processing normalized files
    dirs[:] = [
        d for d in dirs
        if os.path.abspath(os.path.join(root, d)) != abs_output_dir
    ]

    for filename in files:

        if not filename.lower().endswith((".tsv", ".cns")):
            continue

        in_path = os.path.join(root, filename)

        if not os.path.isfile(in_path):
            continue

        keep_cols = KEEP_COLUMNS + ["SAMPLE"] if re.search(r'jabcontool', root, re.IGNORECASE) else KEEP_COLUMNS

        # Output is always <original_basename>_normalized.tsv in a flat output dir
        name = os.path.splitext(filename)[0]
        out_filename = f"{name}_normalized.tsv"
        out_path = os.path.join(OUTPUT_DIR, out_filename)

        print("Processing:", in_path)

        with open(in_path) as infile:
            lines = infile.readlines()

        if not lines:
            continue

        lines = [line.rstrip("\n") for line in lines]

        if not has_header(lines):
            lines.insert(0, MISSING_HEADER)

        raw_header_cols = [c.strip().upper() for c in lines[0].rstrip("\r\n").split("\t")]
        skip_cn_rounding = "READS.RATIO" in raw_header_cols or "FOLD_CHANGE" in raw_header_cols
        header = normalize_header(lines[0])
        header_cols = header.split("\t")

        col_index = {col.upper(): i for i, col in enumerate(header_cols)}

        if not any(col in col_index for col in keep_cols):
            continue

        seen = set()

        with open(out_path, "w") as outfile:

            outfile.write("\t".join(keep_cols) + "\n")

            for line in lines[1:]:

                if not line.strip():
                    continue

                if line.startswith("@"):
                    continue

                cols = line.strip().split()

                if any(REGEX.match(col) for col in cols[:3]):
                    continue

                row = {}

                for col in keep_cols:
                    if col in col_index and col_index[col] < len(cols):
                        row[col] = cols[col_index[col]]
                    else:
                        row[col] = ""

                # CN rounding + remove neutral
                if row["CN_ESTIMATE"]:
                    try:
                        if skip_cn_rounding:
                            cn = float(row["CN_ESTIMATE"])
                        else:
                            cn = round_half_up(row["CN_ESTIMATE"])
                            if cn == 2:
                                continue
                            row["CN_ESTIMATE"] = str(cn)
                    except Exception:
                        continue

                # TYPE normalization / inference
                if row["TYPE"]:
                    row["TYPE"] = normalize_type(row["TYPE"])
                elif row["CN_ESTIMATE"]:
                    row["TYPE"] = infer_type_from_cn(float(row["CN_ESTIMATE"]))

                out_line = "\t".join(row[col] for col in keep_cols)
                if out_line in seen:
                    continue
                seen.add(out_line)
                outfile.write(out_line + "\n")

print("Variant normalization complete.")