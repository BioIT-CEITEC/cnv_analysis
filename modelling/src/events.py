"""
Phase 4 module (row<->event reconciliation), promoted early in Phase 3.5
because the baseline gate already needs event-level recall.

A "row" is one candidate CNV region in all_samples_merged.tsv. An "event" is
one true simulated variant in ground_truth_events.tsv. Multiple rows can match
the same event (Phase 3 Task 7: usually 1, sometimes 2-3); some events match
zero rows (55 of 3,414 -- no caller fired on them at all, so they can never be
recovered by any row-based rule or model -- see findings.md Task 7).
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "modules" / "consensus_model"))
from cnv_consensus_model import MIN_OVERLAP_FRAC, overlaps  # noqa: E402


def match_rows_to_events(df, gt):
    """Per-row ground-truth event key, or None if the row matches no event.

    Event key = (sample, chr, start, end, type) of the matched ground-truth
    row -- unique within a sample since duplicate sample names were already
    disambiguated with a run prefix upstream (merge_from_bronco.py).
    """
    gt_by_sample = {
        sample: sub[["chr", "start", "end", "type"]].to_dict("records")
        for sample, sub in gt.groupby("sample")
    }

    keys = []
    for _, row in df.iterrows():
        match = None
        for g in gt_by_sample.get(row["sample"], []):
            if g["chr"] != row["chr"] or g["type"] != row["type"]:
                continue
            if overlaps(row["start"], row["end"], g["start"], g["end"], MIN_OVERLAP_FRAC):
                match = (row["sample"], g["chr"], g["start"], g["end"], g["type"])
                break
        keys.append(match)
    return pd.Series(keys, index=df.index)


def event_recall(predicted_mask, event_keys, gt):
    """Fraction of ALL ground-truth events (denominator = len(gt), per Phase 3
    Task 6's decision) that have >=1 matched row where predicted_mask is True.

    event_keys: the Series returned by match_rows_to_events(df, gt) -- passed
    in rather than recomputed, since it's identical for every rule evaluated
    against the same df/gt and is not cheap to build per-row.
    """
    matched_and_predicted = set(event_keys[predicted_mask].dropna())
    return len(matched_and_predicted) / len(gt)


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).parent))
    from load import load_ground_truth, load_merged  # noqa: E402

    df = load_merged()
    gt = load_ground_truth()
    keys = match_rows_to_events(df, gt)
    n_matched_events = keys.dropna().nunique()
    print(f"Rows matched to an event: {keys.notna().sum()}/{len(df)}")
    print(f"Distinct events matched by >=1 row: {n_matched_events}/{len(gt)}")
    print(f"Oracle event recall (predict every row true): {event_recall(pd.Series(True, index=df.index), keys, gt):.4f}")
