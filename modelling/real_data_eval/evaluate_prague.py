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
from cnv_consensus_model import ALL_CALLERS, load_ground_truth_dir, load_sample_calls  # noqa: E402

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
        features.insert(0, "sample", sample_dir.name)
        frames.append(features)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def ground_truth_frame(gt_dict):
    rows = [{"sample": sample, "chr": e["chr"], "start": e["start"], "end": e["end"], "type": e["type"]}
            for sample, events in gt_dict.items() for e in events]
    return pd.DataFrame(rows)


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
    y_true = event_keys.notna().astype(int)
    y_prob = df["consensus_score"].values

    pr_auc = round(average_precision_score(y_true, y_prob), 4)
    threshold = best_f2_threshold(y_true, y_prob)
    predicted_mask = pd.Series(y_prob >= threshold, index=df.index)

    results = {
        "beta": BETA,
        "source": "real Prague cohort (non-simulated) -- NOT used for training",
        "n_samples": int(df["sample"].nunique()),
        "n_candidate_rows": int(len(df)),
        "n_true_events": int(len(gt)),
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


if __name__ == "__main__":
    main()
