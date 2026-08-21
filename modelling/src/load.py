"""
Phase 3 data loading + schema validation for the consensus-model project.

Validates the hazards called out in IMPLEMENTATION_PLAN.md's hazard table
before any EDA or feature work touches the data:
  - exactly 8 caller columns, strictly {0,1}
  - matches_gt present and non-null on every row
  - ground truth loaded separately (all_samples_merged.tsv is candidate-only:
    it can only contain rows where >=1 caller fired, so it cannot serve as the
    recall denominator on its own -- see eda.py task 6/7).
"""

from pathlib import Path

import pandas as pd

CALLER_COLUMNS = ["cnvkit", "freec", "exomeDepth", "panelcnMOPS",
                   "cnMOPS", "XHMM", "conifer", "gatk"]

MERGED_COLUMNS = (["sample", "chr", "start", "end", "type", "callers_supporting"]
                   + CALLER_COLUMNS + ["matches_gt", "gt_gene"])

DEFAULT_MERGED_PATH = Path(__file__).parent.parent / "data" / "raw" / "all_samples_merged.tsv"
DEFAULT_GT_PATH = Path(__file__).parent.parent / "data" / "raw" / "ground_truth_events.tsv"


def load_merged(path=DEFAULT_MERGED_PATH):
    """Load and validate all_samples_merged.tsv. Raises on any hazard-table violation."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Regenerate with:\n"
            f"  python modelling/merge_from_bronco.py --out {path} "
            f"--gt-out {DEFAULT_GT_PATH}"
        )

    df = pd.read_csv(path, sep="\t", dtype={"sample": str, "chr": str, "type": str})

    missing = set(MERGED_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(f"Missing expected column(s): {sorted(missing)}")

    present_callers = [c for c in CALLER_COLUMNS if c in df.columns]
    if len(present_callers) != 8:
        raise ValueError(
            f"Expected exactly 8 caller columns, found {len(present_callers)}: {present_callers}"
        )

    for c in CALLER_COLUMNS:
        bad = set(df[c].unique()) - {0, 1}
        if bad:
            raise ValueError(f"Caller column '{c}' has non-binary values: {sorted(bad)}")

    if df["matches_gt"].isna().any():
        n = int(df["matches_gt"].isna().sum())
        raise ValueError(f"matches_gt is null on {n} row(s) -- must be populated on every row")
    bad_labels = set(df["matches_gt"].unique()) - {0, 1}
    if bad_labels:
        raise ValueError(f"matches_gt has non-binary values: {sorted(bad_labels)}")

    return df


def load_ground_truth(path=DEFAULT_GT_PATH):
    """Load the flat ground-truth event list (sample, chr, start, end, type, name).

    This is the independent source of truth for event counts and recall
    denominators -- all_samples_merged.tsv cannot be used for this because it
    only contains rows where at least one caller fired (see module docstring).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Regenerate with:\n"
            f"  python modelling/merge_from_bronco.py --out {DEFAULT_MERGED_PATH} --gt-out {path}"
        )
    df = pd.read_csv(path, sep="\t", dtype={"sample": str, "chr": str, "type": str, "name": str})
    expected = {"sample", "chr", "start", "end", "type", "name"}
    missing = expected - set(df.columns)
    if missing:
        raise ValueError(f"Missing expected column(s) in ground truth file: {sorted(missing)}")
    bad_types = set(df["type"].unique()) - {"DEL", "DUP"}
    if bad_types:
        raise ValueError(f"Ground truth has non-DEL/DUP type value(s): {sorted(bad_types)}")
    return df


if __name__ == "__main__":
    merged = load_merged()
    gt = load_ground_truth()
    print(f"all_samples_merged.tsv: {len(merged)} rows, {merged['sample'].nunique()} samples -- OK")
    print(f"ground_truth_events.tsv: {len(gt)} events, {gt['sample'].nunique()} samples -- OK")
