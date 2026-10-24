"""
Phase 6 -- hyperparameter tuning & optimization (IMPLEMENTATION_PLAN.md).

Tunes XGBoost's hyperparameters via a grouped-CV random search that is
deliberately independent of the Phase 4 persisted `fold` column (different
seed and fold count), so hyperparameter selection can't leak into the
honest out-of-fold metric computed afterwards with those persisted folds.
Then selects a fixed operating threshold against a confirmed precision
target (>=0.90, maximize recall -- confirmed with the user this session)
and persists the tuned model.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve
from sklearn.model_selection import RandomizedSearchCV, StratifiedGroupKFold
from xgboost import XGBClassifier

sys.path.insert(0, str(Path(__file__).parent))
from baseline import BETA  # noqa: E402
from events import match_rows_to_events  # noqa: E402
from features import FEATURE_COLUMNS  # noqa: E402
from load import load_ground_truth  # noqa: E402
from train import BASELINE_JSON  # noqa: E402
from train import REPORT_JSON as TRAINING_JSON  # noqa: E402
from train import _oof_predict, load_features, score_model  # noqa: E402

SEED = 42
SEARCH_SEED = 123
N_ITER = 25
PRECISION_TARGET = 0.90

REPORT_JSON = Path(__file__).parent.parent / "reports" / "metrics" / "tuning.json"
REPORT_MD = Path(__file__).parent.parent / "reports" / "metrics" / "tuning.md"
MODEL_PATH = Path(__file__).parent.parent / "models" / "xgboost_tuned.json"

PARAM_DISTRIBUTIONS = {
    "max_depth": [3, 4, 5, 6, 8],
    "n_estimators": [100, 200, 300, 500],
    "learning_rate": [0.01, 0.03, 0.05, 0.1, 0.2],
    "min_child_weight": [1, 3, 5, 10],
    "subsample": [0.6, 0.8, 1.0],
    "colsample_bytree": [0.6, 0.8, 1.0],
}
# scale_pos_weight is searched as a multiplier on the dataset's own
# negative/positive ratio, not an absolute value -- an absolute value
# wouldn't mean the same thing on a different class balance.
SPW_MULTIPLIERS = [0.5, 0.75, 1.0, 1.5, 2.0]


def search_hyperparameters(df):
    """Randomized search over PARAM_DISTRIBUTIONS x SPW_MULTIPLIERS, scored
    on PR-AUC via a 4-fold grouped-by-sample CV independent of the Phase 4
    persisted fold column (different seed, different fold count) -- so the
    search can't leak into the honest OOF metric computed afterwards.
    Returns the fitted RandomizedSearchCV (best_params_/best_estimator_
    already available) and the top-5 configs for auditability.
    """
    X = df[FEATURE_COLUMNS]
    y = df["matches_gt"]
    base_ratio = (y == 0).sum() / (y == 1).sum()

    param_distributions = dict(PARAM_DISTRIBUTIONS)
    param_distributions["scale_pos_weight"] = [base_ratio * m for m in SPW_MULTIPLIERS]

    cv_splits = list(
        StratifiedGroupKFold(n_splits=4, shuffle=True, random_state=SEARCH_SEED)
        .split(X, y, groups=df["sample"])
    )

    search = RandomizedSearchCV(
        XGBClassifier(random_state=SEED),
        param_distributions=param_distributions,
        n_iter=N_ITER,
        cv=cv_splits,
        scoring="average_precision",
        random_state=SEED,
        n_jobs=-1,
    )
    search.fit(X, y)

    cv_df = pd.DataFrame(search.cv_results_).sort_values("mean_test_score", ascending=False).head(5)
    top5 = [
        {**{k.replace("param_", ""): v for k, v in row.items() if k.startswith("param_")},
         "mean_pr_auc": round(float(row["mean_test_score"]), 4)}
        for _, row in cv_df.iterrows()
    ]
    return search, top5


def select_threshold_for_precision(y_true, y_prob, target=PRECISION_TARGET):
    """Highest-recall point on the PR curve with precision >= target. If no
    point clears the target, report the max precision actually achieved
    instead of failing -- this is a real possibility worth surfacing, not
    an assumption that the target is always reachable."""
    precision, recall, thresholds = precision_recall_curve(y_true, y_prob)
    precision, recall = precision[:-1], recall[:-1]  # last point has no threshold
    meets_target = precision >= target
    if meets_target.any():
        idx = np.where(meets_target)[0]
        best_idx = idx[np.argmax(recall[idx])]
        target_met = True
    else:
        best_idx = int(np.argmax(precision))
        target_met = False
    return {
        "threshold": float(thresholds[best_idx]),
        "precision": float(precision[best_idx]),
        "recall": float(recall[best_idx]),
        "target_met": target_met,
    }


def render_markdown(results):
    lines = ["# Phase 6 -- Hyperparameter Tuning & Optimization", "",
             f"F-beta = F{results['beta']} (comparison metric only -- the operating "
             f"threshold below uses the confirmed precision target instead).", ""]

    lines += ["## Best hyperparameters", "",
             "(4-fold grouped-by-sample CV, seed=123, scored on PR-AUC -- "
             "independent of the Phase 4 reporting folds, seed=42)", ""]
    lines.append("```json")
    lines.append(json.dumps(results["best_params"], indent=2))
    lines.append("```")
    lines.append("")

    lines += ["## Top 5 configs from the search (auditability)", ""]
    keys = list(results["cv_results_top5"][0].keys())
    lines.append("| " + " | ".join(keys) + " |")
    lines.append("|" + "---|" * len(keys))
    for row in results["cv_results_top5"]:
        lines.append("| " + " | ".join(str(row[k]) for k in keys) + " |")
    lines.append("")

    lines += ["## Honest OOF comparison (Phase 4 persisted folds -- same procedure as Phase 5)", ""]
    lines.append("| model | row_precision | row_recall | row_F2 | event_recall | pr_auc |")
    lines.append("|---|---|---|---|---|---|")
    for label, r in (("xgboost_tuned", results["xgboost_tuned"]),
                     ("xgboost_untuned (Phase 5)", results["xgboost_untuned"]),
                     ("baseline (gatk fired)", results["baseline_best"])):
        pr_auc = r.get("pr_auc") if r.get("pr_auc") is not None else "-"
        lines.append(f"| {label} | {r['row_precision']} | {r['row_recall']} | "
                     f"{r['row_f2']} | {r['event_recall']} | {pr_auc} |")
    lines.append("")

    lines += [f"## Operating threshold (target: precision >= {results['precision_target']})", ""]
    op = results["operating_point"]
    if op["target_met"]:
        lines.append(f"**threshold = {op['threshold']:.4f}** -- precision={op['precision']:.4f}, "
                     f"recall={op['recall']:.4f} at this point on the tuned model's OOF PR curve.")
    else:
        lines.append(f"**Target not met.** Best achievable precision on the OOF PR curve is "
                     f"{op['precision']:.4f} (recall={op['recall']:.4f}) at threshold="
                     f"{op['threshold']:.4f} -- reported as the operating point since "
                     f"{results['precision_target']} isn't reachable on this data.")
    lines.append("")
    lines.append(f"Model persisted -> `{results['model_path']}`")
    lines.append("")

    return "\n".join(lines)


def main():
    df = load_features()
    gt = load_ground_truth()
    y_true = df["matches_gt"]

    print("Matching rows to ground-truth events (reused from events.py)...")
    event_keys = match_rows_to_events(df, gt)

    print(f"Randomized search ({N_ITER} configs x 4-fold grouped CV, scoring=PR-AUC)...")
    search, cv_results_top5 = search_hyperparameters(df)
    best_params = search.best_params_

    print("Honest OOF evaluation of the tuned hyperparameters (Phase 4 persisted folds)...")
    def fit_predict_fold(train, test):
        model = XGBClassifier(random_state=SEED, **best_params)
        model.fit(train[FEATURE_COLUMNS], train["matches_gt"])
        return model.predict_proba(test[FEATURE_COLUMNS])[:, 1]
    tuned_oof_prob = _oof_predict(df, fit_predict_fold)
    tuned_result = score_model("xgboost_tuned", y_true, tuned_oof_prob, event_keys, gt)

    print(f"Selecting operating threshold for precision >= {PRECISION_TARGET}...")
    operating_point = select_threshold_for_precision(y_true, tuned_oof_prob)

    print("Persisting tuned model (refit on full data by the search)...")
    MODEL_PATH.parent.mkdir(parents=True, exist_ok=True)
    search.best_estimator_.save_model(str(MODEL_PATH))

    training = json.loads(TRAINING_JSON.read_text(encoding="utf-8"))
    baseline_best = json.loads(BASELINE_JSON.read_text(encoding="utf-8"))["best_overall"]

    results = {
        "beta": BETA,
        "precision_target": PRECISION_TARGET,
        "best_params": {k: (float(v) if isinstance(v, np.floating) else v)
                        for k, v in best_params.items()},
        "cv_results_top5": cv_results_top5,
        "xgboost_tuned": tuned_result,
        "xgboost_untuned": training["xgboost"],
        "baseline_best": baseline_best,
        "operating_point": operating_point,
        "model_path": str(MODEL_PATH),
    }

    REPORT_JSON.parent.mkdir(parents=True, exist_ok=True)
    REPORT_JSON.write_text(json.dumps(results, indent=2), encoding="utf-8")
    REPORT_MD.write_text(render_markdown(results), encoding="utf-8")

    print(f"\nWritten -> {REPORT_JSON}")
    print(f"Written -> {REPORT_MD}")
    print(f"Model saved -> {MODEL_PATH}")
    print(f"\nTuned OOF:   row_F2={tuned_result['row_f2']}, event_recall={tuned_result['event_recall']}")
    print(f"Untuned OOF: row_F2={training['xgboost']['row_f2']}, event_recall={training['xgboost']['event_recall']}")
    print(f"Baseline:    row_F2={baseline_best['row_f2']}, event_recall={baseline_best['event_recall']}")
    met = "" if operating_point["target_met"] else "  [TARGET NOT MET -- best achievable shown]"
    print(f"\nOperating point (target precision >= {PRECISION_TARGET}): "
          f"threshold={operating_point['threshold']:.4f}, "
          f"precision={operating_point['precision']:.4f}, "
          f"recall={operating_point['recall']:.4f}{met}")


if __name__ == "__main__":
    main()
