"""
Phase 4 -- feature engineering (IMPLEMENTATION_PLAN.md).

Drops the leakage/redundant columns identified in Phase 3
(findings.md "Leakage / redundancy"), engineers the model's input features,
and assigns GroupKFold-by-sample folds. Output is the same row-per-candidate
table plus new columns -- `chr`/`start`/`end`/`type`/`sample` stay in the
frame so events.py can still reconcile rows to ground-truth events downstream.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

sys.path.insert(0, str(Path(__file__).parent))
from load import CALLER_COLUMNS, load_merged  # noqa: E402

DROP_COLUMNS = ["callers_supporting", "gt_gene"]

# Default modeling feature set. `chr` is deliberately left out here (kept in
# the persisted frame, not in this list) -- Phase 3/4 expect it to be
# low-value and the with/without ablation is a Phase 5/6 concern.
FEATURE_COLUMNS = CALLER_COLUMNS + ["agreement_count", "type_DUP",
                                     "log_span_width", "type_conflict"]

SEED = 42
N_FOLDS = 5

OUT_PATH = Path(__file__).parent.parent / "data" / "processed" / "features.tsv"


def build_features(df):
    """Drop leakage/redundant columns and add the engineered feature columns."""
    df = df.drop(columns=DROP_COLUMNS)
    assert not (set(DROP_COLUMNS) & set(df.columns)), "leakage columns survived the drop"

    df["agreement_count"] = df[CALLER_COLUMNS].sum(axis=1)
    df["type_DUP"] = (df["type"] == "DUP").astype(int)
    df["log_span_width"] = np.log1p(df["end"] - df["start"])
    df["type_conflict"] = (
        df.groupby(["sample", "chr", "start", "end"])["type"].transform("nunique") > 1
    )
    return df


def add_folds(df, n_folds=N_FOLDS, seed=SEED):
    """Assign a `fold` column via GroupKFold-by-sample, stratified on matches_gt."""
    splitter = StratifiedGroupKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    df = df.copy()
    df["fold"] = -1
    for fold, (_, holdout_idx) in enumerate(
        splitter.split(df, df["matches_gt"], groups=df["sample"])
    ):
        df.iloc[holdout_idx, df.columns.get_loc("fold")] = fold
    assert (df["fold"] >= 0).all(), "every row must get a fold"
    return df


if __name__ == "__main__":
    merged = load_merged()
    features = build_features(merged)

    expected_agreement = {1: 46767, 2: 4580, 3: 1363, 4: 1164,
                           5: 793, 6: 410, 7: 120, 8: 24}
    actual_agreement = features["agreement_count"].value_counts().to_dict()
    assert actual_agreement == expected_agreement, (
        f"agreement_count distribution drifted from findings.md Task 5: "
        f"expected {expected_agreement}, got {actual_agreement}"
    )

    features = add_folds(features)
    folds_per_sample = features.groupby("sample")["fold"].nunique()
    assert (folds_per_sample == 1).all(), "a sample was split across folds"

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    features.to_csv(OUT_PATH, sep="\t", index=False)

    print(f"Written -> {OUT_PATH} ({len(features)} rows, {features['sample'].nunique()} samples)")
    print("\nFold sizes and positive rate:")
    print(features.groupby("fold").agg(n_rows=("matches_gt", "size"),
                                        n_samples=("sample", "nunique"),
                                        positive_rate=("matches_gt", "mean")))
