"""
Feature-selection experiment: XGBoost with `freec` and `cnMOPS` dropped from
the feature set -- flagged by the user and a colleague, from reviewing the
Phase 5/6 results, as not contributing meaningfully. Retrains, re-tunes, and
compares against the existing all-8-caller results.

This file does NOT modify train.py, tune.py, or features.py -- it only
imports their generic (feature-list-agnostic) helpers and defines its own
reduced feature list, so the original all-feature implementation and its
results stay intact in case that's what ends up deployed.
"""

import json
import sys
from pathlib import Path

import numpy as np
from sklearn.model_selection import RandomizedSearchCV, StratifiedGroupKFold
from xgboost import XGBClassifier

sys.path.insert(0, str(Path(__file__).parent))
from baseline import BETA  # noqa: E402
from events import match_rows_to_events  # noqa: E402
from features import FEATURE_COLUMNS as ALL_FEATURE_COLUMNS  # noqa: E402
from load import load_ground_truth  # noqa: E402
from train import BASELINE_JSON  # noqa: E402
from train import REPORT_JSON as TRAINING_JSON  # noqa: E402
from train import _oof_predict, load_features, score_model  # noqa: E402
from tune import PARAM_DISTRIBUTIONS, SPW_MULTIPLIERS  # noqa: E402
from tune import REPORT_JSON as TUNING_JSON  # noqa: E402
from tune import (N_ITER, PRECISION_TARGET, SEARCH_SEED, SEED,
                  select_threshold_for_precision)  # noqa: E402

DROPPED_CALLERS = ["freec", "cnMOPS"]
FEATURE_COLUMNS = [c for c in ALL_FEATURE_COLUMNS if c not in DROPPED_CALLERS]
assert not (set(DROPPED_CALLERS) & set(FEATURE_COLUMNS)), \
    "dropped callers leaked back into the reduced feature list"

REPORT_JSON = Path(__file__).parent.parent / "reports" / "metrics" / "no_freec_cnmops.json"
REPORT_MD = Path(__file__).parent.parent / "reports" / "metrics" / "no_freec_cnmops.md"
MODEL_PATH = Path(__file__).parent.parent / "models" / "xgboost_tuned_no_freec_cnmops.json"


def oof_predict_xgb_reduced(df):
    """Same shape as train.oof_predict_xgb, closed over the reduced
    FEATURE_COLUMNS instead of the full 8-caller list."""
    def fit_predict_fold(train, test):
        n_pos = (train["matches_gt"] == 1).sum()
        n_neg = (train["matches_gt"] == 0).sum()
        model = XGBClassifier(scale_pos_weight=n_neg / n_pos, random_state=SEED)
        model.fit(train[FEATURE_COLUMNS], train["matches_gt"])
        return model.predict_proba(test[FEATURE_COLUMNS])[:, 1]
    return _oof_predict(df, fit_predict_fold)


def search_hyperparameters_reduced(df):
    """Same shape and search config as tune.search_hyperparameters, closed
    over the reduced FEATURE_COLUMNS."""
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
    return search


def render_markdown(results):
    lines = ["# Feature-selection experiment -- drop freec + cnMOPS", "",
             f"F-beta = F{results['beta']}. Dropped features: "
             f"{', '.join(results['dropped_features'])}.", ""]

    lines += ["## Five-way comparison (all OOF, Phase 4 persisted folds)", ""]
    lines.append("| model | row_precision | row_recall | row_F2 | event_recall | pr_auc |")
    lines.append("|---|---|---|---|---|---|")
    rows = [
        ("baseline (gatk fired)", results["baseline_best"]),
        ("xgboost, all features (Phase 5)", results["xgboost_full_untuned"]),
        ("xgboost tuned, all features (Phase 6)", results["xgboost_full_tuned"]),
        ("xgboost, reduced features (no freec/cnMOPS)", results["xgboost_reduced_untuned"]),
        ("xgboost tuned, reduced features", results["xgboost_reduced_tuned"]),
    ]
    for label, r in rows:
        pr_auc = r.get("pr_auc") if r.get("pr_auc") is not None else "-"
        lines.append(f"| {label} | {r['row_precision']} | {r['row_recall']} | "
                     f"{r['row_f2']} | {r['event_recall']} | {pr_auc} |")
    lines.append("")

    full_f2 = results["xgboost_full_tuned"]["row_f2"]
    reduced_f2 = results["xgboost_reduced_tuned"]["row_f2"]
    delta = reduced_f2 - full_f2
    verdict = ("no material difference" if abs(delta) < 0.002 else
              ("an improvement" if delta > 0 else "a regression"))
    lines += ["## Verdict", ""]
    lines.append(f"Dropping {' and '.join(results['dropped_features'])} changes tuned OOF "
                 f"row_F2 from {full_f2:.4f} to {reduced_f2:.4f} ({delta:+.4f}) -- "
                 f"**{verdict}** compared to keeping all 8 callers.")
    lines.append("")

    lines += ["## Best hyperparameters, reduced feature set", "",
             "(same search config as Phase 6: 4-fold grouped-by-sample CV, seed=123, "
             "scored on PR-AUC)", ""]
    lines.append("```json")
    lines.append(json.dumps(results["best_params_reduced"], indent=2))
    lines.append("```")
    lines.append("")

    op = results["operating_point_reduced"]
    lines += [f"## Operating threshold, reduced model (target: precision >= {PRECISION_TARGET})", ""]
    if op["target_met"]:
        lines.append(f"**threshold = {op['threshold']:.4f}** -- precision={op['precision']:.4f}, "
                     f"recall={op['recall']:.4f}.")
    else:
        lines.append(f"**Target not met.** Best achievable precision is {op['precision']:.4f} "
                     f"(recall={op['recall']:.4f}) at threshold={op['threshold']:.4f}.")
    lines.append("")
    lines.append(f"Model persisted -> `{results['model_path']}`")
    lines.append("")

    return "\n".join(lines)


def main():
    df = load_features()
    gt = load_ground_truth()
    y_true = df["matches_gt"]

    print(f"Reduced feature set ({len(FEATURE_COLUMNS)} columns, dropped "
          f"{DROPPED_CALLERS}): {FEATURE_COLUMNS}")
    print("Matching rows to ground-truth events (reused from events.py)...")
    event_keys = match_rows_to_events(df, gt)

    print("XGBoost OOF training on reduced feature set (retrain step)...")
    reduced_prob = oof_predict_xgb_reduced(df)
    reduced_result = score_model("xgboost_reduced", y_true, reduced_prob, event_keys, gt)

    print(f"Randomized search on reduced feature set ({N_ITER} configs x 4-fold "
          f"grouped CV, scoring=PR-AUC)...")
    search = search_hyperparameters_reduced(df)
    best_params = search.best_params_

    print("Honest OOF evaluation of the tuned reduced-feature model...")
    def fit_predict_fold(train, test):
        model = XGBClassifier(random_state=SEED, **best_params)
        model.fit(train[FEATURE_COLUMNS], train["matches_gt"])
        return model.predict_proba(test[FEATURE_COLUMNS])[:, 1]
    tuned_prob = _oof_predict(df, fit_predict_fold)
    tuned_result = score_model("xgboost_reduced_tuned", y_true, tuned_prob, event_keys, gt)

    print(f"Selecting operating threshold for precision >= {PRECISION_TARGET}...")
    operating_point = select_threshold_for_precision(y_true, tuned_prob, PRECISION_TARGET)

    print("Persisting tuned reduced-feature model (refit on full data by the search)...")
    MODEL_PATH.parent.mkdir(parents=True, exist_ok=True)
    search.best_estimator_.save_model(str(MODEL_PATH))

    existing_training = json.loads(TRAINING_JSON.read_text(encoding="utf-8"))
    existing_tuning = json.loads(TUNING_JSON.read_text(encoding="utf-8"))
    baseline_best = json.loads(BASELINE_JSON.read_text(encoding="utf-8"))["best_overall"]

    results = {
        "beta": BETA,
        "dropped_features": DROPPED_CALLERS,
        "feature_columns": FEATURE_COLUMNS,
        "xgboost_reduced_untuned": reduced_result,
        "best_params_reduced": {k: (float(v) if isinstance(v, np.floating) else v)
                                for k, v in best_params.items()},
        "xgboost_reduced_tuned": tuned_result,
        "operating_point_reduced": operating_point,
        "xgboost_full_untuned": existing_training["xgboost"],
        "xgboost_full_tuned": existing_tuning["xgboost_tuned"],
        "baseline_best": baseline_best,
        "model_path": str(MODEL_PATH),
    }

    REPORT_JSON.parent.mkdir(parents=True, exist_ok=True)
    REPORT_JSON.write_text(json.dumps(results, indent=2), encoding="utf-8")
    REPORT_MD.write_text(render_markdown(results), encoding="utf-8")

    print(f"\nWritten -> {REPORT_JSON}")
    print(f"Written -> {REPORT_MD}")
    print(f"Model saved -> {MODEL_PATH}")
    print(f"\nXGBoost (all features, untuned):     row_F2={existing_training['xgboost']['row_f2']}, "
          f"event_recall={existing_training['xgboost']['event_recall']}")
    print(f"XGBoost (all features, tuned):        row_F2={existing_tuning['xgboost_tuned']['row_f2']}, "
          f"event_recall={existing_tuning['xgboost_tuned']['event_recall']}")
    print(f"XGBoost (reduced features, untuned):  row_F2={reduced_result['row_f2']}, "
          f"event_recall={reduced_result['event_recall']}")
    print(f"XGBoost (reduced features, tuned):    row_F2={tuned_result['row_f2']}, "
          f"event_recall={tuned_result['event_recall']}")
    print(f"Baseline:                              row_F2={baseline_best['row_f2']}, "
          f"event_recall={baseline_best['event_recall']}")


if __name__ == "__main__":
    main()
