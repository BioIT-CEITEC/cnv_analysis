"""
Train the consensus model directly from Alejandro's BRONCO simulated cohort,
without staging/copying any files.

The 569 training samples live scattered across 22 run folders:
  mapped_final/BRONCO_<run>/simulated_data/structural_varcalls/<sample>/<caller>/...
  mapped_final/BRONCO_<run>/simulated_data/mapped/<sample>_selected_regions.tsv

train_consensus_model.py assumes a single flat --sim_dir/--gt_dir (one
subdirectory per sample). This script reuses its exact training logic
(BoundedConsensusModel, train_model, caller stats, feature matrix) but resolves
each sample's directory from a {sample_name: Path} map instead, built here from
enumerating the 22 run folders.

Usage:
  python train_from_bronco.py --out_model cnv_consensus_model.json --out_prefix cnv_train
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "modules" / "consensus_model"))
from cnv_consensus_model import (  # noqa: E402
    ALL_CALLERS, DEFAULT_CALLERS, CALLER_FILE_PATTERNS, MIN_OVERLAP_FRAC,
    norm_chr, merge_gt_variants, load_sample_calls, _best_hit, overlaps,
)
from train_consensus_model import (  # noqa: E402
    split_samples, train_model,
)

import numpy as np

BRONCO_BASE = Path(
    "/mnt/data/ceitec_cfg2/710000-CEITEC/713000-cmm/713016-bioit"
    "/base/workspace/alejandro/test_bronco/mapped_final"
)


# ---------------------------------------------------------------------------
# Enumeration: build {sample_name: structural_varcalls dir} + ground truth,
# reading directly from the 22 scattered run folders. No copying.
# ---------------------------------------------------------------------------

def discover_samples(base=BRONCO_BASE):
    runs = sorted(d.name for d in base.iterdir() if d.is_dir() and d.name.startswith("BRONCO_"))

    raw_pairs = []  # (run, sample, sv_dir, gt_tsv_path)
    for run in runs:
        mapped_dir = base / run / "simulated_data" / "mapped"
        sv_dir = base / run / "simulated_data" / "structural_varcalls"
        if not mapped_dir.is_dir() or not sv_dir.is_dir():
            continue
        for gt_tsv in mapped_dir.glob("*_selected_regions.tsv"):
            sample = gt_tsv.name[: -len("_selected_regions.tsv")]
            sample_sv_dir = sv_dir / sample
            if sample_sv_dir.is_dir():
                raw_pairs.append((run, sample, sample_sv_dir, gt_tsv))

    name_to_runs = defaultdict(list)
    for run, sample, _, _ in raw_pairs:
        name_to_runs[sample].append(run)
    dupe_names = {n for n, rs in name_to_runs.items() if len(rs) > 1}

    sample_dirs = {}
    gt = {}
    all_samples = set()
    for run, sample, sv_dir, gt_tsv in raw_pairs:
        final_name = f"{run}_{sample}" if sample in dupe_names else sample
        sample_dirs[final_name] = sv_dir

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
            })
        gt[final_name] = merge_gt_variants(variants)
        all_samples.add(final_name)

    print(f"Runs found: {len(runs)}")
    print(f"Total samples discovered: {len(all_samples)}")
    print(f"Duplicate names disambiguated with run prefix: {len(dupe_names)}")
    return sample_dirs, gt, all_samples


def detect_callers_scattered(sample_dirs):
    found = set()
    for sv_dir in sample_dirs.values():
        for sub in sv_dir.iterdir():
            if sub.is_dir() and sub.name in ALL_CALLERS:
                found.add(sub.name)
        for caller in ALL_CALLERS:
            if caller in found:
                continue
            pattern = CALLER_FILE_PATTERNS.get(caller, f"*_{caller}.tsv")
            if next(sv_dir.rglob(pattern), None) is not None:
                found.add(caller)
            elif pattern.endswith(".vcf.gz") and next(sv_dir.rglob(pattern[:-3]), None) is not None:
                found.add(caller)
        if len(found) == len(ALL_CALLERS):
            break
    return [c for c in ALL_CALLERS if c in found]


# ---------------------------------------------------------------------------
# Scattered-directory equivalents of compute_caller_stats / build_feature_matrix
# (identical logic to train_consensus_model.py, only the sample->dir lookup differs)
# ---------------------------------------------------------------------------

def compute_caller_stats_scattered(sample_dirs, gt, sample_set, callers):
    rec = {c: {"TP": 0, "FN": 0, "FP": 0} for c in callers}
    rec_del = {c: {"TP": 0, "FN": 0, "FP": 0} for c in callers}
    rec_dup = {c: {"TP": 0, "FN": 0, "FP": 0} for c in callers}

    for sample in sorted(sample_set):
        sv_dir = sample_dirs.get(sample)
        if sv_dir is None or not sv_dir.is_dir():
            continue
        sample_gt = gt.get(sample, [])
        sim_calls = load_sample_calls(sv_dir, callers)

        for caller in callers:
            calls = sim_calls.get(caller, [])
            for g in sample_gt:
                hit, _ = _best_hit(g, calls)
                bucket = rec_del[caller] if g["type"] == "DEL" else rec_dup[caller]
                if hit:
                    rec[caller]["TP"] += 1
                    bucket["TP"] += 1
                else:
                    rec[caller]["FN"] += 1
                    bucket["FN"] += 1
            for call in calls:
                is_tp = any(call["chr"] == g["chr"] and call.get("type") == g["type"]
                            and overlaps(call["start"], call["end"],
                                         g["start"], g["end"], MIN_OVERLAP_FRAC)
                            for g in sample_gt)
                if is_tp:
                    continue
                rec[caller]["FP"] += 1
                bucket = rec_del[caller] if call.get("type") == "DEL" else rec_dup[caller]
                bucket["FP"] += 1

    def _m(r):
        tp, fp, fn = r["TP"], r["FP"], r["FN"]
        prec = tp / (tp + fp) if (tp + fp) > 0 else float("nan")
        rec_ = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
        f1 = (2 * prec * rec_ / (prec + rec_)
              if not np.isnan(prec) and not np.isnan(rec_) and (prec + rec_) > 0
              else float("nan"))
        return tp, fp, fn, prec, rec_, f1

    rows = []
    for caller in callers:
        for vtype, table in (("ALL", rec), ("DEL", rec_del), ("DUP", rec_dup)):
            tp, fp, fn, prec, rec_, f1 = _m(table[caller])
            rows.append({"caller": caller, "variant_type": vtype,
                         "TP": tp, "FP": fp, "FN": fn,
                         "Precision": round(prec, 4) if not np.isnan(prec) else "NA",
                         "Recall": round(rec_, 4) if not np.isnan(rec_) else "NA",
                         "F1": round(f1, 4) if not np.isnan(f1) else "NA"})
    return pd.DataFrame(rows)


def build_feature_matrix_scattered(sample_dirs, gt, samples, callers):
    rows, labels, meta = [], [], []

    for sample in sorted(samples):
        sv_dir = sample_dirs.get(sample)
        if sv_dir is None or not sv_dir.is_dir():
            continue
        sample_gt = gt.get(sample, [])
        sim_calls = load_sample_calls(sv_dir, callers)

        for g in sample_gt:
            row = []
            for c in callers:
                hit, score = _best_hit(g, sim_calls.get(c, []))
                row.extend([hit, score])
            rows.append(row)
            labels.append(1)
            meta.append((sample, g["chr"], g["start"], g["end"], g["type"], "TP"))

        for caller in callers:
            for call in sim_calls.get(caller, []):
                is_tp = any(call["chr"] == g["chr"] and call["type"] == g["type"]
                            and overlaps(call["start"], call["end"],
                                         g["start"], g["end"], MIN_OVERLAP_FRAC)
                            for g in sample_gt)
                if is_tp:
                    continue
                row = []
                for c2 in callers:
                    if c2 == caller:
                        row.extend([1, call.get("score", 0.0)])
                    else:
                        hit2, score2 = _best_hit(call, sim_calls.get(c2, []))
                        row.extend([hit2, score2])
                rows.append(row)
                labels.append(0)
                meta.append((sample, call["chr"], call["start"], call["end"],
                             call["type"], "FP"))

    feature_names = []
    for c in callers:
        feature_names.extend([f"{c}_hit", f"{c}_score"])
    X = np.array(rows, dtype=float) if rows else np.zeros((0, 2 * len(callers)))
    y = np.array(labels, dtype=int)
    return X, y, feature_names, meta


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Train CNV consensus model directly from scattered BRONCO run folders")
    ap.add_argument("--out_model", default="cnv_consensus_model.json")
    ap.add_argument("--out_prefix", default="cnv_train")
    ap.add_argument("--callers", nargs="+", default=None)
    ap.add_argument("--weight_cap", type=float, default=1.0)
    ap.add_argument("--cap_mode", choices=("per_caller", "per_coef"), default="per_caller")
    ap.add_argument("--l2", type=float, default=0.01)
    ap.add_argument("--test_size", type=float, default=0.2)
    ap.add_argument("--random_seed", type=int, default=42)
    ap.add_argument("--score_cutoff", type=float, default=None)
    args = ap.parse_args()

    print("Discovering samples across BRONCO run folders...")
    sample_dirs, gt, all_samples = discover_samples()
    total_variants = sum(len(v) for v in gt.values())
    print(f"  {len(all_samples)} samples, {total_variants} true variants")

    print("Detecting callers...")
    available = detect_callers_scattered(sample_dirs)
    callers = args.callers or [c for c in DEFAULT_CALLERS if c in available]
    missing = [c for c in callers if c not in available]
    if missing:
        print(f"  WARNING: requested callers not found in data: {', '.join(missing)}")
    callers = [c for c in callers if c in available]
    print(f"  Using callers: {', '.join(callers)}")
    if not callers:
        sys.exit("ERROR: no requested callers present in the data.")

    print(f"\nSplitting samples ({int((1-args.test_size)*100)}% train / "
          f"{int(args.test_size*100)}% test, seed={args.random_seed})...")
    train_samples, test_samples = split_samples(gt, all_samples, args.test_size, args.random_seed)
    print(f"  Train: {len(train_samples)}  |  Test: {len(test_samples)}")
    print(f"  Train GT variants: {sum(len(gt.get(s, [])) for s in train_samples)}")

    split_rows = [{"sample": s, "split": "train"} for s in sorted(train_samples)]
    split_rows += [{"sample": s, "split": "test"} for s in sorted(test_samples)]
    split_path = f"{args.out_prefix}_train_test_split.tsv"
    pd.DataFrame(split_rows).to_csv(split_path, sep="\t", index=False)
    print(f"  Split saved -> {split_path}")

    print("\nComputing per-caller statistics on training samples (raw genome-wide FP)...")
    stats = compute_caller_stats_scattered(sample_dirs, gt, train_samples, callers)
    stats_all = stats[stats["variant_type"] == "ALL"].drop(columns="variant_type")
    stats_type = stats[stats["variant_type"] != "ALL"]
    stats_all.to_csv(f"{args.out_prefix}_caller_stats.tsv", sep="\t", index=False)
    stats_type.to_csv(f"{args.out_prefix}_caller_stats_by_type.tsv", sep="\t", index=False)
    print(f"\n  {'Caller':<14} {'TP':>5} {'FP':>6} {'FN':>5} {'Prec':>7} {'Rec':>7} {'F1':>7}")
    print("  " + "-" * 54)
    for _, r in stats_all.iterrows():
        print(f"  {r['caller']:<14} {r['TP']:>5} {r['FP']:>6} {r['FN']:>5} "
              f"{str(r['Precision']):>7} {str(r['Recall']):>7} {str(r['F1']):>7}")

    print("\nBuilding training feature matrix ([hit, score] per caller)...")
    X, y, feature_names, _ = build_feature_matrix_scattered(sample_dirs, gt, train_samples, callers)
    print(f"  Shape: {X.shape}  |  positives: {int(y.sum())}  negatives: {int((y == 0).sum())}")
    if y.sum() < 5:
        sys.exit("ERROR: too few positives to train.")

    cap_desc = ("combined hit+score" if args.cap_mode == "per_caller" else "each coefficient")
    print("\nTraining StandardScaler + logistic regression "
          f"({cap_desc} in [0, {args.weight_cap}], l2={args.l2})...")
    model, metrics, pr_df = train_model(X, y, feature_names, callers,
                                        args.weight_cap, args.l2, args.cap_mode)
    caller_weights = model.caller_weight_dict()
    cutoff = args.score_cutoff if args.score_cutoff is not None else metrics["suggested_cutoff"]

    pr_df.to_csv(f"{args.out_prefix}_pr_curve.tsv", sep="\t", index=False)

    print("\n=== Training metrics ===")
    for k, v in metrics.items():
        print(f"  {k}: {v}")
    cap_label = (f"combined <= {args.weight_cap}" if args.cap_mode == "per_caller"
                 else f"each coef <= {args.weight_cap}")
    print(f"\n=== Caller weights (non-negative, {cap_label}) ===")
    print(f"  {'Caller':<14} {'hit_coef':>9} {'score_coef':>11} {'combined':>10}")
    print("  " + "-" * 48)
    for c, w in sorted(caller_weights.items(), key=lambda x: -x[1]["combined"]):
        print(f"  {c:<14} {w['hit_coef']:>9.4f} {w['score_coef']:>11.4f} {w['combined']:>10.4f}")
    print(f"  {'(intercept)':<14} {model.intercept_:>9.4f}")
    print(f"\n  Suggested cutoff (max-F1 on train): {cutoff:.4f}")

    import json
    import pickle
    out = {
        "model_type": "BoundedConsensusModel",
        "strategy": "cnv_analysis (StandardScaler + LR on [hit,score]); "
                    "coefficients constrained non-negative and capped",
        "metrics": metrics,
        "caller_weights": caller_weights,
        "feature_weights": model.feature_weight_dict(),
        "intercept": round(float(model.intercept_), 6),
        "cutoff": cutoff,
        "callers": callers,
        "weight_cap": args.weight_cap,
        "cap_mode": args.cap_mode,
        "l2": args.l2,
        "min_overlap_frac": MIN_OVERLAP_FRAC,
        "train_samples": sorted(train_samples),
        "test_samples": sorted(test_samples),
    }
    with open(args.out_model, "w") as fh:
        json.dump(out, fh, indent=2)
    pkl_path = args.out_model.replace(".json", ".pkl")
    with open(pkl_path, "wb") as fh:
        pickle.dump(model, fh)
    print(f"\nModel saved -> {args.out_model}")
    print(f"Pipeline saved -> {pkl_path}")


if __name__ == "__main__":
    main()
