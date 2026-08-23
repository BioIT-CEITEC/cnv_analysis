"""
Phase 5 -- model selection & training (IMPLEMENTATION_PLAN.md).

Trains XGBoost (primary) and a plain L2 logistic regression (reference)
under the Phase 4 GroupKFold-by-sample splits, collects honest out-of-fold
(OOF) predictions, and gates the result against the Phase 3.5 baseline
(reports/metrics/baseline.json). Per the plan's pre-committed decision: if
the model doesn't beat the baseline materially, that's a legitimate finding
to report, not a failure to fix by tuning harder (that's what Phase 6 is
for, and only if this gate passes).
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, precision_recall_curve
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

sys.path.insert(0, str(Path(__file__).parent))
from baseline import BETA, evaluate_rule  # noqa: E402
from events import match_rows_to_events  # noqa: E402
from features import FEATURE_COLUMNS, OUT_PATH as FEATURES_PATH  # noqa: E402
from load import load_ground_truth  # noqa: E402

SEED = 42

REPORT_JSON = Path(__file__).parent.parent / "reports" / "metrics" / "training.json"
REPORT_MD = Path(__file__).parent.parent / "reports" / "metrics" / "training.md"
BASELINE_JSON = Path(__file__).parent.parent / "reports" / "metrics" / "baseline.json"

# Phase 3 Task 4 per-caller precision, findings.md -- the model's feature
# importance is checked against this ranking, not recomputed from raw data.
CALLER_PRECISION = {
    "XHMM": 0.8282, "gatk": 0.7842, "cnMOPS": 0.4618, "exomeDepth": 0.314,
    "cnvkit": 0.2929, "conifer": 0.2366, "panelcnMOPS": 0.1611, "freec": 0.0479,
}


def load_features():
    if not FEATURES_PATH.exists():
        raise FileNotFoundError(
            f"{FEATURES_PATH} not found. Run: python modelling/src/features.py"
        )
    df = pd.read_csv(FEATURES_PATH, sep="\t",
                     dtype={"sample": str, "chr": str, "type": str})
    df["type_conflict"] = df["type_conflict"].astype(bool)
    return df


def best_f2_threshold(y_true, y_prob):
    """Threshold maximizing F-beta (beta=BETA) from the PR curve."""
    precision, recall, thresholds = precision_recall_curve(y_true, y_prob)
    beta2 = BETA ** 2
    denom = beta2 * precision + recall
    f_beta = np.where(denom > 0, (1 + beta2) * precision * recall / denom, 0.0)
    best_idx = int(np.argmax(f_beta[:-1]))  # last PR-curve point has no threshold
    return float(thresholds[best_idx])


def _oof_predict(df, fit_predict_fold):
    """Run fit_predict_fold(train_df, test_df) once per fold, assemble one
    out-of-fold probability Series covering every row exactly once."""
    prob = pd.Series(np.nan, index=df.index)
    for fold in sorted(df["fold"].unique()):
        train = df[df["fold"] != fold]
        test = df[df["fold"] == fold]
        prob.loc[test.index] = fit_predict_fold(train, test)
    assert prob.notna().all(), "every row must get an out-of-fold prediction"
    return prob


def oof_predict_xgb(df):
    def fit_predict_fold(train, test):
        n_pos = (train["matches_gt"] == 1).sum()
        n_neg = (train["matches_gt"] == 0).sum()
        model = XGBClassifier(scale_pos_weight=n_neg / n_pos, random_state=SEED)
        model.fit(train[FEATURE_COLUMNS], train["matches_gt"])
        return model.predict_proba(test[FEATURE_COLUMNS])[:, 1]
    return _oof_predict(df, fit_predict_fold)


def oof_predict_logreg(df):
    def fit_predict_fold(train, test):
        model = Pipeline([
            ("scale", StandardScaler()),
            ("clf", LogisticRegression(penalty="l2", class_weight="balanced",
                                       max_iter=1000, random_state=SEED)),
        ])
        model.fit(train[FEATURE_COLUMNS], train["matches_gt"])
        return model.predict_proba(test[FEATURE_COLUMNS])[:, 1]
    return _oof_predict(df, fit_predict_fold)


def score_model(name, y_true, oof_prob, event_keys, gt):
    pr_auc = round(average_precision_score(y_true, oof_prob), 4)
    threshold = best_f2_threshold(y_true, oof_prob)
    result = evaluate_rule(y_true, oof_prob >= threshold, event_keys, gt,
                           f"{name} (OOF, threshold={threshold:.4f})")
    result["pr_auc"] = pr_auc
    result["threshold"] = round(threshold, 4)
    return result


def feature_importance_check(df):
    """Fit one XGBoost model on the full frame (importance only, never used
    for scoring) and check its caller-flag ranking against CALLER_PRECISION."""
    model = XGBClassifier(
        scale_pos_weight=(df["matches_gt"] == 0).sum() / (df["matches_gt"] == 1).sum(),
        random_state=SEED,
    )
    model.fit(df[FEATURE_COLUMNS], df["matches_gt"])
    importances = dict(zip(FEATURE_COLUMNS, model.feature_importances_.tolist()))

    callers = list(CALLER_PRECISION)
    rho, pvalue = spearmanr([CALLER_PRECISION[c] for c in callers],
                            [importances[c] for c in callers])
    return {
        "importances": {k: round(v, 4) for k, v in importances.items()},
        "caller_ranking_spearman_rho": round(float(rho), 4),
        "caller_ranking_spearman_pvalue": round(float(pvalue), 4),
    }


def render_markdown(results):
    lines = ["# Phase 5 -- Model Selection & Training", "",
             f"F-beta = F{results['beta']} (same metric as the Phase 3.5 baseline gate).", ""]

    lines += ["## Out-of-fold model comparison", ""]
    lines.append("| model | threshold | pr_auc | row_precision | row_recall | row_F2 | event_recall |")
    lines.append("|---|---|---|---|---|---|---|")
    for r in (results["xgboost"], results["logreg"], results["baseline_best"]):
        label = r.get("description", r.get("rule", "baseline"))
        threshold = r.get("threshold", "-")
        pr_auc = r.get("pr_auc") if r.get("pr_auc") is not None else "-"
        lines.append(f"| {label} | {threshold} | {pr_auc} | "
                     f"{r['row_precision']} | {r['row_recall']} | {r['row_f2']} | "
                     f"{r['event_recall']} |")
    lines.append("")

    lines += ["## Gate", ""]
    verdict = "PASSED" if results["gate_passed"] else "FAILED"
    lines.append(f"**{verdict}** -- XGBoost OOF row_F2={results['xgboost']['row_f2']} vs "
                 f"baseline row_F2={results['baseline_best']['row_f2']} "
                 f"({results['baseline_best']['description']}).")
    if not results["gate_passed"]:
        lines.append("\nPer the plan's pre-committed decision: the model did not beat the "
                     "baseline materially. Do not proceed to Phase 6 tuning -- report the "
                     "baseline rule as sufficient.")
    lines.append("")

    lines += ["## Feature importance (XGBoost, full-data fit, gain-based) vs "
              "Phase 3 caller precision", ""]
    lines.append(f"Spearman rho (caller precision vs importance) = "
                 f"{results['feature_importance']['caller_ranking_spearman_rho']} "
                 f"(p={results['feature_importance']['caller_ranking_spearman_pvalue']})")
    lines.append("")
    lines.append("| feature | importance |")
    lines.append("|---|---|")
    for feat, imp in sorted(results["feature_importance"]["importances"].items(),
                            key=lambda kv: -kv[1]):
        lines.append(f"| {feat} | {imp} |")
    lines.append("")

    return "\n".join(lines)


def main():
    df = load_features()
    gt = load_ground_truth()
    y_true = df["matches_gt"]

    print("Matching rows to ground-truth events (reused from events.py)...")
    event_keys = match_rows_to_events(df, gt)

    print("XGBoost OOF training (5 folds)...")
    xgb_prob = oof_predict_xgb(df)
    xgb_result = score_model("xgboost", y_true, xgb_prob, event_keys, gt)

    print("Logistic regression OOF training (5 folds)...")
    logreg_prob = oof_predict_logreg(df)
    logreg_result = score_model("logreg", y_true, logreg_prob, event_keys, gt)

    print("Feature importance vs Phase 3 caller precision...")
    importance = feature_importance_check(df)

    baseline_best = json.loads(BASELINE_JSON.read_text(encoding="utf-8"))["best_overall"]
    gate_passed = xgb_result["row_f2"] > baseline_best["row_f2"]

    results = {
        "beta": BETA,
        "xgboost": xgb_result,
        "logreg": logreg_result,
        "feature_importance": importance,
        "baseline_best": baseline_best,
        "gate_passed": gate_passed,
    }

    REPORT_JSON.parent.mkdir(parents=True, exist_ok=True)
    REPORT_JSON.write_text(json.dumps(results, indent=2), encoding="utf-8")
    REPORT_MD.write_text(render_markdown(results), encoding="utf-8")

    print(f"\nWritten -> {REPORT_JSON}")
    print(f"Written -> {REPORT_MD}")
    print(f"\nXGBoost OOF:  row_F2={xgb_result['row_f2']}, event_recall={xgb_result['event_recall']}")
    print(f"Logreg OOF:   row_F2={logreg_result['row_f2']}, event_recall={logreg_result['event_recall']}")
    print(f"Baseline:     row_F2={baseline_best['row_f2']}, event_recall={baseline_best['event_recall']}")
    print(f"\nGate {'PASSED' if gate_passed else 'FAILED'}")


if __name__ == "__main__":
    main()
