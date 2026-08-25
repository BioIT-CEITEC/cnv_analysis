"""
score_samples_xgb.py -- score one sample's CNV caller calls with the tuned
XGBoost consensus model (no freec/cnMOPS in FEATURE_COLUMNS). The model was
produced by modelling/src/train_no_freec_cnmops.py; this script must stay in
sync with that file's feature construction.

Unlike score_samples.py (which loops over an sv_dir of many sample
subdirectories), this scores exactly ONE sample per invocation -- matching
the one-Nextflow-task-per-sample granularity XGBOOST_CONSENSUS_SCORE runs at.

Expected --sample_dir layout: a directory holding that one sample's caller
output files (matched by CALLER_FILE_PATTERNS in cnv_consensus_model.py).

IMPORTANT: agreement_count is computed from ALL 8 callers in ALL_CALLERS,
even though only 6 of them (excluding freec/cnMOPS) are fed to the model
directly. This matches exactly how the training data was built --
modelling/src/features.py computes agreement_count BEFORE
train_no_freec_cnmops.py drops freec/cnMOPS from FEATURE_COLUMNS. Do not
"simplify" this to 6 callers -- it would silently mismatch what the model
was trained on.

Output TSV columns (tab-separated, one row per merged candidate region) --
identical shape to score_samples.py's output, so main.nf's awk reshape can
be reused unchanged:
    sample | chr | start | end | type | n_callers | callers_supporting | consensus_score | cn_label
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb

# Reuse the shared cnv_consensus_model.py (repo root modules/consensus_model/)
# unchanged -- same sys.path convention already used by modelling/merge_from_bronco.py.
# This file lives at modelling/xgboost_pipeline/modules/xgboost_consensus_model/,
# four levels under the repo root.
sys.path.insert(0, str(Path(__file__).parents[4] / "modules" / "consensus_model"))
from cnv_consensus_model import ALL_CALLERS, _merge_candidates, load_sample_calls  # noqa: E402

EMPTY_HEADER = "sample\tchr\tstart\tend\ttype\tn_callers\tcallers_supporting\tconsensus_score\tcn_label\n"

# Must match modelling/src/train_no_freec_cnmops.py's FEATURE_COLUMNS exactly,
# same order: CALLER_COLUMNS (modelling/src/load.py) minus freec/cnMOPS, plus
# the four engineered features from modelling/src/features.py.
FEATURE_COLUMNS = ["cnvkit", "exomeDepth", "panelcnMOPS", "XHMM", "conifer", "gatk",
                   "agreement_count", "type_DUP", "log_span_width", "type_conflict"]


def build_candidate_table(calls_by_caller):
    """Same candidate-building step as cnv_consensus_model.score_calls(), up
    to _merge_candidates() -- our own feature vector is built from there
    instead of the [hit, score] pairs BoundedConsensusModel expects."""
    candidates = []
    for caller, calls in calls_by_caller.items():
        for call in calls:
            candidates.append({"chr": call["chr"], "start": call["start"],
                               "end": call["end"], "type": call["type"],
                               "caller": caller})
    if not candidates:
        return pd.DataFrame()
    return _merge_candidates(pd.DataFrame(candidates))


def build_features(merged):
    """Reproduce modelling/src/features.py's build_features() from the
    production candidate table: explode callers_supporting into the 8 named
    binary columns, then the 4 engineered features."""
    df = merged.copy()
    for caller in ALL_CALLERS:
        df[caller] = df["callers_supporting"].apply(
            lambda s, c=caller: int(c in s.split(","))
        )
    df["agreement_count"] = df[ALL_CALLERS].sum(axis=1)
    df["type_DUP"] = (df["type"] == "DUP").astype(int)
    df["log_span_width"] = np.log1p(df["end"] - df["start"])
    df["type_conflict"] = (
        df.groupby(["chr", "start", "end"])["type"].transform("nunique") > 1
    )
    return df


def main():
    ap = argparse.ArgumentParser(description="Score one sample's CNV calls with the tuned XGBoost consensus model")
    ap.add_argument("--sample_dir", required=True, help="Directory of this sample's caller output files")
    ap.add_argument("--sample_name", required=True)
    ap.add_argument("--model", required=True, help="Path to the XGBoost native-format model JSON")
    ap.add_argument("--out", required=True, help="Output TSV path")
    ap.add_argument("--cutoff", type=float, default=None,
                    help="If set, also write a filtered TSV with score >= cutoff")
    args = ap.parse_args()

    calls_by_caller = load_sample_calls(Path(args.sample_dir), ALL_CALLERS)
    merged = build_candidate_table(calls_by_caller)

    if merged.empty:
        with open(args.out, "w") as fh:
            fh.write(EMPTY_HEADER)
        return

    features = build_features(merged)

    model = xgb.XGBClassifier()
    model.load_model(args.model)
    features["consensus_score"] = model.predict_proba(features[FEATURE_COLUMNS])[:, 1]

    out = features[["chr", "start", "end", "type", "n_callers",
                    "callers_supporting", "consensus_score", "cn_label"]].copy()
    out.insert(0, "sample", args.sample_name)
    out = out.sort_values("consensus_score", ascending=False)
    out.to_csv(args.out, sep="\t", index=False)

    if args.cutoff is not None:
        filtered = out[out["consensus_score"] >= args.cutoff]
        filtered_path = args.out.replace(".tsv", f"_filtered_{args.cutoff:.4f}.tsv")
        filtered.to_csv(filtered_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
