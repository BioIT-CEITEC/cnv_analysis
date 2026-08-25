"""
Real-data evaluation: score the tuned XGBoost consensus model (no freec/
cnMOPS, modelling/models/xgboost_tuned_no_freec_cnmops.json) against the
Prague cohort's real, non-simulated ground truth.

INFERENCE ONLY. This script never calls .fit()/.train() on this data -- it
loads the already-trained model and scores it. Do not add training here:
this is real patient data, not simulated, and must not be used to fit
anything.

Source data (not committed to this repo -- real patient data):
    /mnt/data/ceitec_cfg2/710000-CEITEC/713000-cmm/713016-bioit/base/workspace/
    alejandro/test_bronco/test_prague/
        <sample>/<caller>/...                        caller output (8 callers)
        selected_regions/<sample>_selected_regions.tsv   ground truth

Outputs go to results/ next to this script, which is gitignored -- those
files contain real sample identifiers and are not meant to leave this
machine.

Usage (same conda env as every other modelling/src script):
    conda activate cnv_consensus_model
    python modelling/real_data_eval/evaluate_prague.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import (average_precision_score, fbeta_score,
                             precision_recall_curve, precision_score,
                             recall_score)

REPO_ROOT = Path(__file__).parent.parent.parent
PRAGUE_DIR = Path("/mnt/data/ceitec_cfg2/710000-CEITEC/713000-cmm/713016-bioit/base/"
                  "workspace/alejandro/test_bronco/test_prague")
GT_DIR = PRAGUE_DIR / "selected_regions"
MODEL_PATH = REPO_ROOT / "modelling" / "models" / "xgboost_tuned_no_freec_cnmops.json"
OUT_DIR = Path(__file__).parent / "results"

# Reuse cnv_consensus_model.py unchanged, same sys.path convention as
# modelling/merge_from_bronco.py.
sys.path.insert(0, str(REPO_ROOT / "modules" / "consensus_model"))
from cnv_consensus_model import (ALL_CALLERS, MIN_OVERLAP_FRAC,  # noqa: E402
                                 load_ground_truth_dir, load_sample_calls, overlaps)

# Reuse the exact candidate-building/feature-construction the deployed
# pipeline uses, so this evaluation matches production scoring exactly.
sys.path.insert(0, str(REPO_ROOT / "modelling" / "xgboost_pipeline" / "modules" / "xgboost_consensus_model"))
from score_samples_xgb import FEATURE_COLUMNS, build_candidate_table, build_features  # noqa: E402

sys.path.insert(0, str(REPO_ROOT / "modelling" / "src"))
from events import event_recall, match_rows_to_events  # noqa: E402

BETA = 2  # same recall-weighted F-beta used throughout this project
IGNORE_DIRS = {"selected_regions"}


def score_all_samples(prague_dir, model):
    frames = []
    sample_dirs = sorted(d for d in prague_dir.iterdir()
                         if d.is_dir() and d.name not in IGNORE_DIRS)
    for sample_dir in sample_dirs:
        calls_by_caller = load_sample_calls(sample_dir, ALL_CALLERS)
        merged = build_candidate_table(calls_by_caller)
        if merged.empty:
            continue
        features = build_features(merged)
        features["consensus_score"] = model.predict_proba(features[FEATURE_COLUMNS])[:, 1]
        # Ground truth uses the short sample id (e.g. "007_CZE3xBRCA1328_run137"),
        # not the full directory name (e.g. "..._run137.picard.sorted.RG.rmdup")
        # -- same short/long split as params.samples in the repo-root nextflow.config.
        features.insert(0, "sample", sample_dir.name.split(".")[0])
        frames.append(features)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def ground_truth_frame(gt_dict):
    rows = [{"sample": sample, "chr": e["chr"], "start": e["start"], "end": e["end"], "type": e["type"]}
            for sample, events in gt_dict.items() for e in events]
    return pd.DataFrame(rows)


def per_caller_event_recall(prague_dir, gt_dict):
    """For each true event, which of the 8 callers INDIVIDUALLY detected it
    (>=MIN_OVERLAP_FRAC overlap, matching type)? Diagnoses whether the
    real-data caller-coverage gap (Finding 1) is caller-specific or across
    the board -- separate from the model, which never sees an event a
    caller didn't propose a candidate for in the first place.
    """
    sample_dirs = {d.name.split(".")[0]: d for d in prague_dir.iterdir()
                   if d.is_dir() and d.name not in IGNORE_DIRS}
    hits = {c: 0 for c in ALL_CALLERS}
    event_rows = []
    total_events = 0

    for sample, events in gt_dict.items():
        sample_dir = sample_dirs.get(sample)
        if sample_dir is None:
            continue
        calls_by_caller = load_sample_calls(sample_dir, ALL_CALLERS)
        for e in events:
            total_events += 1
            detected_by = []
            for caller in ALL_CALLERS:
                hit = any(
                    call["chr"] == e["chr"] and call["type"] == e["type"]
                    and overlaps(call["start"], call["end"], e["start"], e["end"], MIN_OVERLAP_FRAC)
                    for call in calls_by_caller.get(caller, [])
                )
                if hit:
                    hits[caller] += 1
                    detected_by.append(caller)
            event_rows.append({
                "sample": sample, "chr": e["chr"], "type": e["type"],
                "span_width": e["end"] - e["start"],
                "n_callers_detecting": len(detected_by),
                "detected_by": ",".join(detected_by) if detected_by else "NONE",
            })

    per_caller = {c: {"n_events_detected": hits[c],
                      "recall": round(hits[c] / total_events, 4) if total_events else 0.0}
                  for c in ALL_CALLERS}
    return per_caller, pd.DataFrame(event_rows), total_events


def best_f2_threshold(y_true, y_prob):
    precision, recall, thresholds = precision_recall_curve(y_true, y_prob)
    beta2 = BETA ** 2
    denom = beta2 * precision + recall
    f_beta = np.where(denom > 0, (1 + beta2) * precision * recall / denom, 0.0)
    best_idx = int(np.argmax(f_beta[:-1]))  # last PR-curve point has no threshold
    return float(thresholds[best_idx])


def render_markdown(results):
    lines = [
        "# Real-data evaluation -- Prague cohort (XGBoost tuned, no freec/cnMOPS)",
        "",
        "**Inference only -- this model was NOT trained or fit on this data.** "
        "Source is real (non-simulated) patient data. Do not commit raw outputs "
        "from this folder.",
        "",
        f"- Samples scored: {results['n_samples']}",
        f"- Candidate rows: {results['n_candidate_rows']}",
        f"- True events (ground truth): {results['n_true_events']}",
        f"- Operating threshold (F{results['beta']}-optimal): {results['threshold']}",
        "",
        f"**Event-recall ceiling: {results['event_recall_ceiling']}** "
        f"({results['n_unreachable_events']} of {results['n_true_events']} true events have "
        f"zero candidate-row support from any of the 8 callers -- unreachable by any "
        f"row-based model, this one included). Compare `event recall` below against this "
        f"ceiling, not against 1.0, to separate a model problem from a caller-coverage problem.",
        "",
        "| metric | value |",
        "|---|---|",
        f"| PR-AUC | {results['pr_auc']} |",
        f"| row precision | {results['row_precision']} |",
        f"| row recall | {results['row_recall']} |",
        f"| row F{results['beta']} | {results['row_f2']} |",
        f"| event recall | {results['event_recall']} |",
        "",
    ]
    return "\n".join(lines)


def render_caller_coverage_markdown(per_caller, event_df, total_events):
    lines = [
        "# Real-data caller-coverage diagnostic -- Prague cohort", "",
        "**Inference/diagnostic only -- no training.** Investigates *why* the "
        "event-recall ceiling collapses on real data (Finding 1 in "
        "reports/discussion.md): is it caller-specific, or across the board?",
        "", f"Total true events checked: {total_events}", "",
        "## Per-caller recall against real ground truth", "",
        "| caller | events detected | recall |", "|---|---|---|",
    ]
    for caller, r in sorted(per_caller.items(), key=lambda kv: -kv[1]["recall"]):
        lines.append(f"| {caller} | {r['n_events_detected']} | {r['recall']} |")
    lines.append("")

    unreachable = event_df[event_df["n_callers_detecting"] == 0]
    reachable = event_df[event_df["n_callers_detecting"] > 0]
    lines += ["## Span width: reachable vs unreachable events", "",
             "Do callers miss real events because they're small, or is it independent of size?",
             "", "| group | n | min | median | max |", "|---|---|---|---|---|"]
    for label, sub in (("unreachable (0 callers)", unreachable), ("reachable (>=1 caller)", reachable)):
        if len(sub) == 0:
            lines.append(f"| {label} | 0 | - | - | - |")
        else:
            w = sub["span_width"]
            lines.append(f"| {label} | {len(sub)} | {int(w.min())} | {int(w.median())} | {int(w.max())} |")
    lines.append("")

    return "\n".join(lines)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Loading real ground truth from {GT_DIR}...")
    gt_dict, gt_samples = load_ground_truth_dir(GT_DIR)
    gt = ground_truth_frame(gt_dict)
    print(f"  {len(gt_samples)} samples, {len(gt)} true events")

    print(f"Loading model from {MODEL_PATH}...")
    model = xgb.XGBClassifier()
    model.load_model(str(MODEL_PATH))

    print(f"Scoring real Prague samples from {PRAGUE_DIR} (inference only, no training)...")
    df = score_all_samples(PRAGUE_DIR, model)
    print(f"  {len(df)} candidate rows across {df['sample'].nunique()} samples")

    print("Matching candidate rows to real ground-truth events...")
    event_keys = match_rows_to_events(df, gt)
    y_true = (event_keys.apply(len) > 0).astype(int)
    y_prob = df["consensus_score"].values

    # Event-recall ceiling (same check as Phase 3's blocking check on BRONCO,
    # baseline.py): events with zero candidate-row support from ANY caller
    # can never be recovered by any row-based model, regardless of quality.
    # Separates "the model missed it" from "no caller ever proposed it".
    all_true_mask = pd.Series(True, index=df.index)
    event_recall_ceiling = event_recall(all_true_mask, event_keys, gt)
    matched_events = set()
    for matches in event_keys:
        matched_events.update(matches)
    n_unreachable_events = len(gt) - len(matched_events)

    pr_auc = round(average_precision_score(y_true, y_prob), 4)
    threshold = best_f2_threshold(y_true, y_prob)
    predicted_mask = pd.Series(y_prob >= threshold, index=df.index)

    results = {
        "beta": BETA,
        "source": "real Prague cohort (non-simulated) -- NOT used for training",
        "n_samples": int(df["sample"].nunique()),
        "n_candidate_rows": int(len(df)),
        "n_true_events": int(len(gt)),
        "n_unreachable_events": int(n_unreachable_events),
        "event_recall_ceiling": round(event_recall_ceiling, 4),
        "threshold": round(threshold, 4),
        "pr_auc": pr_auc,
        "row_precision": round(precision_score(y_true, predicted_mask, zero_division=0), 4),
        "row_recall": round(recall_score(y_true, predicted_mask, zero_division=0), 4),
        "row_f2": round(fbeta_score(y_true, predicted_mask, beta=BETA, zero_division=0), 4),
        "event_recall": round(event_recall(predicted_mask, event_keys, gt), 4),
    }

    (OUT_DIR / "prague_eval.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    (OUT_DIR / "prague_eval.md").write_text(render_markdown(results), encoding="utf-8")

    df_out = df.copy()
    df_out["matches_gt"] = y_true.values
    df_out.to_csv(OUT_DIR / "prague_scored_candidates.tsv", sep="\t", index=False)

    print(f"\nWritten -> {OUT_DIR / 'prague_eval.json'}")
    print(f"Written -> {OUT_DIR / 'prague_eval.md'}")
    print(f"Written -> {OUT_DIR / 'prague_scored_candidates.tsv'}")
    print(f"\nRow precision={results['row_precision']}  recall={results['row_recall']}  "
          f"F{BETA}={results['row_f2']}  event_recall={results['event_recall']}  "
          f"PR-AUC={pr_auc}")

    print("\nPer-caller recall against real ground truth (diagnosing Finding 1)...")
    per_caller, event_df, total_events = per_caller_event_recall(PRAGUE_DIR, gt_dict)
    coverage_json = {"total_events": total_events, "per_caller": per_caller}
    (OUT_DIR / "prague_caller_coverage.json").write_text(
        json.dumps(coverage_json, indent=2), encoding="utf-8")
    (OUT_DIR / "prague_caller_coverage.md").write_text(
        render_caller_coverage_markdown(per_caller, event_df, total_events), encoding="utf-8")
    event_df.to_csv(OUT_DIR / "prague_event_detection_detail.tsv", sep="\t", index=False)

    print(f"Written -> {OUT_DIR / 'prague_caller_coverage.json'}")
    print(f"Written -> {OUT_DIR / 'prague_caller_coverage.md'}")
    print(f"Written -> {OUT_DIR / 'prague_event_detection_detail.tsv'}")
    print("\nPer-caller recall:")
    for caller, r in sorted(per_caller.items(), key=lambda kv: -kv[1]["recall"]):
        print(f"  {caller:<14} {r['n_events_detected']}/{total_events}  recall={r['recall']}")


if __name__ == "__main__":
    main()
