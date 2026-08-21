"""
Produce one merged calls table across all BRONCO simulated samples, without
training or scoring anything -- just the candidate-region merge
(cnv_consensus_model.py's _merge_candidates) applied per sample, concatenated,
reshaped to match template.txt's column layout.

Reads directly from the 22 scattered run folders (no copying), reusing the
sample discovery from train_from_bronco.py.

Usage:
  python merge_from_bronco.py --out all_samples_merged.tsv
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "modules" / "consensus_model"))
from cnv_consensus_model import (  # noqa: E402
    ALL_CALLERS, MIN_OVERLAP_FRAC, norm_chr, load_sample_calls,
    _merge_candidates, overlaps,
)

sys.path.insert(0, str(Path(__file__).parent))
from train_from_bronco import (  # noqa: E402
    BRONCO_BASE, discover_samples, detect_callers_scattered,
)


def discover_ground_truth_with_names(base=BRONCO_BASE):
    """
    Same (run, sample) enumeration + duplicate-name disambiguation as
    discover_samples(), but keeps the 'name' (gene) column instead of running
    it through merge_gt_variants() -- the *_selected_regions.tsv files are
    already one row per simulated variant, not exon fragments needing merge.
    Returns {final_sample_name: [{chr, start, end, type, name}, ...]}.
    """
    runs = sorted(d.name for d in base.iterdir() if d.is_dir() and d.name.startswith("BRONCO_"))

    raw_pairs = []  # (run, sample, gt_tsv_path)
    for run in runs:
        mapped_dir = base / run / "simulated_data" / "mapped"
        sv_dir = base / run / "simulated_data" / "structural_varcalls"
        if not mapped_dir.is_dir() or not sv_dir.is_dir():
            continue
        for gt_tsv in mapped_dir.glob("*_selected_regions.tsv"):
            sample = gt_tsv.name[: -len("_selected_regions.tsv")]
            if (sv_dir / sample).is_dir():
                raw_pairs.append((run, sample, gt_tsv))

    name_to_runs = defaultdict(list)
    for run, sample, _ in raw_pairs:
        name_to_runs[sample].append(run)
    dupe_names = {n for n, rs in name_to_runs.items() if len(rs) > 1}

    gt_named = {}
    for run, sample, gt_tsv in raw_pairs:
        final_name = f"{run}_{sample}" if sample in dupe_names else sample
        df = pd.read_csv(gt_tsv, sep="\t", encoding="utf-8", encoding_errors="replace")
        df.columns = df.columns.str.strip()
        variants = []
        for _, row in df.iterrows():
            vtype = str(row["type"]).strip().upper()
            if vtype not in ("DEL", "DUP"):
                continue
            variants.append({
                "chr": norm_chr(row["chrom"]),
                "start": int(row["start"]),
                "end": int(row["end"]),
                "type": vtype,
                "name": str(row.get("name", "")).strip(),
            })
        gt_named[final_name] = variants
    return gt_named


def export_ground_truth(gt_named, out_path):
    """Flatten {sample: [variants]} to one TSV: sample, chr, start, end, type, name.

    Independent record of the true simulated events, since all_samples_merged.tsv
    only contains rows where >=1 caller fired -- a variant no caller detected at
    all has no row there. This file is the correct recall denominator.
    """
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(["sample", "chr", "start", "end", "type", "name"]) + "\n")
        for sample in sorted(gt_named):
            for v in gt_named[sample]:
                fh.write("\t".join(str(x) for x in
                          [sample, v["chr"], v["start"], v["end"], v["type"], v["name"]]) + "\n")


def best_gt_match(region, gt_variants):
    """(matches_gt, gt_gene) for one merged region against its sample's ground truth."""
    for g in gt_variants:
        if g["chr"] != region["chr"] or g["type"] != region["type"]:
            continue
        if overlaps(region["start"], region["end"], g["start"], g["end"], MIN_OVERLAP_FRAC):
            return 1, g["name"]
    return 0, ""


def merge_sample(sv_dir, callers, gt_variants):
    """Return the merged, template-shaped candidate-region table for one sample."""
    calls_by_caller = load_sample_calls(sv_dir, callers)
    candidates = []
    for caller, calls in calls_by_caller.items():
        for call in calls:
            candidates.append({"chr": call["chr"], "start": call["start"],
                               "end": call["end"], "type": call["type"],
                               "caller": caller, "score": call.get("score", 0.0),
                               "cn": call.get("cn")})
    if not candidates:
        return []

    merged = _merge_candidates(pd.DataFrame(candidates))
    if merged.empty:
        return []
    merged = merged.sort_values(["chr", "start"]).reset_index(drop=True)

    rows = []
    for _, reg in merged.iterrows():
        supporting = set(reg["callers_supporting"].split(","))
        matches_gt, gt_gene = best_gt_match(reg, gt_variants)
        row = {
            "chr": reg["chr"], "start": int(reg["start"]), "end": int(reg["end"]),
            "type": reg["type"], "callers_supporting": reg["callers_supporting"],
        }
        for c in ALL_CALLERS:
            row[c] = 1 if c in supporting else 0
        row["matches_gt"] = matches_gt
        row["gt_gene"] = gt_gene
        rows.append(row)
    return rows


HEADER = (["sample", "chr", "start", "end", "type", "callers_supporting"]
          + ALL_CALLERS + ["matches_gt", "gt_gene"])


def format_callers_supporting(value):
    """Quote only when it lists more than one caller, matching template.txt."""
    return f'"{value}"' if "," in value else value


def write_tsv(rows, out_path):
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(HEADER) + "\n")
        for row in rows:
            values = [str(row[col]) if col != "callers_supporting"
                      else format_callers_supporting(row[col])
                      for col in HEADER]
            fh.write("\t".join(values) + "\n")


def main():
    ap = argparse.ArgumentParser(description="Merge CNV calls across all BRONCO samples (no training/scoring)")
    ap.add_argument("--out", default="all_samples_merged.tsv")
    ap.add_argument("--gt-out", default=None,
                    help="Also export the flat ground-truth event list here (sample, chr, start, end, type, name)")
    ap.add_argument("--callers", nargs="+", default=None,
                    help="Callers to include (default: all detected)")
    args = ap.parse_args()

    print("Discovering samples across BRONCO run folders...")
    sample_dirs, _, all_samples = discover_samples()
    print(f"  {len(all_samples)} samples")

    print("Loading ground truth (with gene names)...")
    gt_named = discover_ground_truth_with_names()

    print("Detecting callers...")
    available = detect_callers_scattered(sample_dirs)
    callers = args.callers or available
    print(f"  Using callers: {', '.join(callers)}")

    all_rows = []
    for i, sample in enumerate(sorted(all_samples), 1):
        sv_dir = sample_dirs[sample]
        gt_variants = gt_named.get(sample, [])
        sample_rows = merge_sample(sv_dir, callers, gt_variants)
        for row in sample_rows:
            row["sample"] = sample
        all_rows.extend(sample_rows)
        if i % 50 == 0:
            print(f"  ... {i}/{len(all_samples)} samples merged")

    write_tsv(all_rows, args.out)
    print(f"\nDone. {len(all_rows)} merged candidate regions across {len(all_samples)} samples.")
    print(f"Written -> {args.out}")

    if args.gt_out:
        export_ground_truth(gt_named, args.gt_out)
        total_gt = sum(len(v) for v in gt_named.values())
        print(f"Ground truth: {total_gt} events across {len(gt_named)} samples -> {args.gt_out}")


if __name__ == "__main__":
    main()
