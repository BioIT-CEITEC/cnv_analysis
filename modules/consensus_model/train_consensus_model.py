"""
train_consensus_model.py  (simulated-variant benchmark)
=======================================================
Train the consensus model on the SIMULATED CNV cohort only (no baseline),
following the reference cnv_analysis strategy — StandardScaler + logistic
regression on [hit, score] features per caller — with one change: every
coefficient is NON-NEGATIVE and capped at `--weight_cap` (default 1.0).

Inputs
------
  --sim_dir   structural_varcalls_simulated_BR/   (caller calls WITH injected variants)
  --gt_dir    selected_regions_simulated/         (*_selected_regions.tsv ground truth)

Outputs (written next to --out_model / with --out_prefix)
--------
  <out_model>.json    caller weights + cutoff + callers + metrics + train/test split
  <out_model>.pkl     pickled BoundedConsensusModel (for scoring / evaluation)
  <out_prefix>_caller_stats.tsv          per-caller TP/FP/FN/Precision/Recall/F1 (train)
  <out_prefix>_caller_stats_by_type.tsv  same, split by DEL / DUP
  <out_prefix>_pr_curve.tsv              consensus precision/recall/F1 across thresholds
  <out_prefix>_train_test_split.tsv      which samples went to train vs test

Example
-------
  python train_consensus_model.py \
      --sim_dir  ~/Documents/model_training_data/structural_varcalls_simulated_BR \
      --gt_dir   ~/Documents/model_training_data/selected_regions_simulated \
      --out_model cnv_consensus_model.json \
      --out_prefix cnv_train \
      --test_size 0.2
"""

import argparse
import json
import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import (average_precision_score, precision_recall_curve,
                             roc_auc_score)
from sklearn.model_selection import StratifiedKFold, train_test_split

sys.path.insert(0, str(Path(__file__).parent))
from cnv_consensus_model import (  # noqa: E402
    DEFAULT_CALLERS, MIN_OVERLAP_FRAC, IGNORE_DIRS,
    load_ground_truth_dir, detect_callers, load_sample_calls,
    _best_hit, overlaps, build_feature_matrix, BoundedConsensusModel,
)

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Sample-level train/test split
# ---------------------------------------------------------------------------

def split_samples(gt, samples, test_size=0.2, random_seed=42):
    """Split at the sample level, stratified by whether the sample has any GT variant."""
    samples = sorted(samples)
    if test_size == 0.0:
        return set(samples), set()
    has_variant = [1 if gt.get(s) else 0 for s in samples]
    if sum(has_variant) < 2 or (len(has_variant) - sum(has_variant)) < 2:
        train, test = train_test_split(samples, test_size=test_size,
                                       random_state=random_seed)
    else:
        train, test = train_test_split(samples, test_size=test_size,
                                       stratify=has_variant, random_state=random_seed)
    return set(train), set(test)


# ---------------------------------------------------------------------------
# Per-caller statistics  (simulated calls vs ground truth, raw genome-wide FP)
# ---------------------------------------------------------------------------

def compute_caller_stats(sim_dir, gt, sample_set, callers):
    """
    For each caller: TP (GT hit), FN (GT missed), FP (call that overlaps no GT
    variant).  Computed directly from the caller's own calls — independent of
    any model.  Returns a long DataFrame with variant_type in {ALL, DEL, DUP}.
    """
    rec = {c: {"TP": 0, "FN": 0, "FP": 0} for c in callers}
    rec_del = {c: {"TP": 0, "FN": 0, "FP": 0} for c in callers}
    rec_dup = {c: {"TP": 0, "FN": 0, "FP": 0} for c in callers}

    for sample in sorted(sample_set):
        sim_sd = Path(sim_dir) / sample
        if not sim_sd.is_dir():
            continue
        sample_gt = gt.get(sample, [])
        sim_calls = load_sample_calls(sim_sd, callers)

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


# ---------------------------------------------------------------------------
# Model training
# ---------------------------------------------------------------------------

def train_model(X, y, feature_names, callers, weight_cap, l2, cap_mode):
    """Fit BoundedConsensusModel with 5-fold CV ROC-AUC; return (model, metrics, pr_df)."""
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_aucs = []
    for tr, va in cv.split(X, y):
        m = BoundedConsensusModel(feature_names, callers, weight_cap=weight_cap,
                                  l2=l2, cap_mode=cap_mode).fit(X[tr], y[tr])
        if y[va].sum() > 0 and (y[va] == 0).sum() > 0:
            cv_aucs.append(roc_auc_score(y[va], m.predict_proba(X[va])[:, 1]))
    cv_aucs = np.array(cv_aucs) if cv_aucs else np.array([float("nan")])

    model = BoundedConsensusModel(feature_names, callers, weight_cap=weight_cap,
                                  l2=l2, cap_mode=cap_mode).fit(X, y)
    y_prob = model.predict_proba(X)[:, 1]

    pr, rec, thr = precision_recall_curve(y, y_prob)
    f1 = 2 * pr * rec / np.where((pr + rec) > 0, pr + rec, 1)
    best_thr = float(thr[int(np.argmax(f1[:-1]))])

    metrics = {
        "cv_roc_auc_mean": round(float(np.nanmean(cv_aucs)), 4),
        "cv_roc_auc_std": round(float(np.nanstd(cv_aucs)), 4),
        "train_roc_auc": round(float(roc_auc_score(y, y_prob)), 4),
        "train_avg_precision": round(float(average_precision_score(y, y_prob)), 4),
        "n_train_positives": int(y.sum()),
        "n_train_negatives": int((y == 0).sum()),
        "weight_cap": weight_cap,
        "l2": l2,
        "suggested_cutoff": round(best_thr, 4),
    }
    pr_df = pd.DataFrame({"threshold": list(thr) + [1.0],
                          "precision": pr.tolist(), "recall": rec.tolist(),
                          "F1": f1.tolist()})
    return model, metrics, pr_df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Train non-negative capped CNV consensus model")
    ap.add_argument("--sim_dir", required=True, help="structural_varcalls_simulated_BR")
    ap.add_argument("--gt_dir", required=True, help="selected_regions_simulated")
    ap.add_argument("--out_model", default="cnv_consensus_model.json")
    ap.add_argument("--out_prefix", default="cnv_train")
    ap.add_argument("--callers", nargs="+", default=None,
                    help=f"Callers to use (default: {' '.join(DEFAULT_CALLERS)})")
    ap.add_argument("--weight_cap", type=float, default=1.0,
                    help="Upper bound for the (non-negative) caller weight (default 1.0)")
    ap.add_argument("--cap_mode", choices=("per_caller", "per_coef"),
                    default="per_caller",
                    help="per_caller: hit+score per caller <= weight_cap so NO weight "
                         "ever exceeds it (default). per_coef: each coefficient "
                         "<= weight_cap (combined can reach 2x).")
    ap.add_argument("--l2", type=float, default=0.01,
                    help="L2 penalty on coefficients (default 0.01)")
    ap.add_argument("--test_size", type=float, default=0.2)
    ap.add_argument("--random_seed", type=int, default=42)
    ap.add_argument("--score_cutoff", type=float, default=None,
                    help="Override the F1-optimal cutoff")
    args = ap.parse_args()

    print("Loading ground truth...")
    gt, gt_samples = load_ground_truth_dir(args.gt_dir)
    total_variants = sum(len(v) for v in gt.values())
    print(f"  {len(gt_samples)} labelled samples, {total_variants} true variants")

    sim_have = {d.name for d in Path(args.sim_dir).iterdir()
                if d.is_dir() and d.name not in IGNORE_DIRS}
    usable = sorted(gt_samples & sim_have)
    dropped = sorted(gt_samples - set(usable))
    print(f"  {len(usable)} samples usable (present in GT + simulated)")
    if dropped:
        print(f"  {len(dropped)} labelled samples dropped (no simulated dir): "
              f"{', '.join(dropped[:8])}{' ...' if len(dropped) > 8 else ''}")

    print("Detecting callers...")
    available = detect_callers(args.sim_dir)
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
    train_samples, test_samples = split_samples(gt, usable, args.test_size, args.random_seed)
    print(f"  Train: {len(train_samples)}  |  Test: {len(test_samples)}")
    print(f"  Train GT variants: {sum(len(gt.get(s, [])) for s in train_samples)}")

    split_rows = [{"sample": s, "split": "train"} for s in sorted(train_samples)]
    split_rows += [{"sample": s, "split": "test"} for s in sorted(test_samples)]
    split_path = f"{args.out_prefix}_train_test_split.tsv"
    pd.DataFrame(split_rows).to_csv(split_path, sep="\t", index=False)
    print(f"  Split saved -> {split_path}")

    print("\nComputing per-caller statistics on training samples "
          "(raw genome-wide FP)...")
    stats = compute_caller_stats(args.sim_dir, gt, train_samples, callers)
    stats_all = stats[stats["variant_type"] == "ALL"].drop(columns="variant_type")
    stats_type = stats[stats["variant_type"] != "ALL"]
    stats_all.to_csv(f"{args.out_prefix}_caller_stats.tsv", sep="\t", index=False)
    stats_type.to_csv(f"{args.out_prefix}_caller_stats_by_type.tsv", sep="\t", index=False)
    print(f"\n  {'Caller':<14} {'TP':>5} {'FP':>6} {'FN':>5} "
          f"{'Prec':>7} {'Rec':>7} {'F1':>7}")
    print("  " + "-" * 54)
    for _, r in stats_all.iterrows():
        print(f"  {r['caller']:<14} {r['TP']:>5} {r['FP']:>6} {r['FN']:>5} "
              f"{str(r['Precision']):>7} {str(r['Recall']):>7} {str(r['F1']):>7}")

    print("\nBuilding training feature matrix ([hit, score] per caller)...")
    X, y, feature_names, _ = build_feature_matrix(args.sim_dir, gt, train_samples, callers)
    print(f"  Shape: {X.shape}  |  positives: {int(y.sum())}  negatives: {int((y==0).sum())}")
    if y.sum() < 5:
        sys.exit("ERROR: too few positives to train.")

    cap_desc = ("combined hit+score" if args.cap_mode == "per_caller"
                else "each coefficient")
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
        print(f"  {c:<14} {w['hit_coef']:>9.4f} {w['score_coef']:>11.4f} "
              f"{w['combined']:>10.4f}")
    print(f"  {'(intercept)':<14} {model.intercept_:>9.4f}")
    print(f"\n  Suggested cutoff (max-F1 on train): {cutoff:.4f}")

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
    print("\nDone. Evaluate with evaluate_consensus_model.py.")


if __name__ == "__main__":
    main()
