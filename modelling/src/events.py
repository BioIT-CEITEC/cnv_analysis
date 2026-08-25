"""
Phase 4 module (row<->event reconciliation), promoted early in Phase 3.5
because the baseline gate already needs event-level recall.

A "row" is one candidate CNV region in all_samples_merged.tsv. An "event" is
one true simulated variant in ground_truth_events.tsv. Multiple rows can match
the same event (Phase 3 Task 7: usually 1, sometimes 2-3); some events match
zero rows (55 of 3,414 -- no caller fired on them at all, so they can never be
recovered by any row-based rule or model -- see findings.md Task 7).

A row can also legitimately match MORE THAN ONE event: one large candidate
call can span several ground-truth events that are themselves >10kb apart
(merge_gt_variants()'s merge-gap threshold splits them into separate events).
Confirmed real case: reports/discussion.md's real-data section, one exomeDepth
call spanning 4 separate EPCAM/MSH2 sub-events. match_rows_to_events/
event_recall must not assume at most one match per row -- an earlier version
did, and silently undercounted event recall/ceiling whenever this happened.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "modules" / "consensus_model"))
from cnv_consensus_model import MIN_OVERLAP_FRAC, overlaps  # noqa: E402


def match_rows_to_events(df, gt):
    """Per-row list of ALL matching ground-truth event keys (possibly empty,
    possibly more than one -- see module docstring).

    Event key = (sample, chr, start, end, type) of a matched ground-truth
    row -- unique within a sample since duplicate sample names were already
    disambiguated with a run prefix upstream (merge_from_bronco.py).
    """
    gt_by_sample = {
        sample: sub[["chr", "start", "end", "type"]].to_dict("records")
        for sample, sub in gt.groupby("sample")
    }

    all_matches = []
    for _, row in df.iterrows():
        matches = []
        for g in gt_by_sample.get(row["sample"], []):
            if g["chr"] != row["chr"] or g["type"] != row["type"]:
                continue
            if overlaps(row["start"], row["end"], g["start"], g["end"], MIN_OVERLAP_FRAC):
                matches.append((row["sample"], g["chr"], g["start"], g["end"], g["type"]))
        all_matches.append(matches)
    return pd.Series(all_matches, index=df.index)


def event_recall(predicted_mask, event_keys, gt):
    """Fraction of ALL ground-truth events (denominator = len(gt), per Phase 3
    Task 6's decision) that have >=1 matched row where predicted_mask is True.
    A row can match multiple events and an event can be matched by multiple
    rows (many-to-many) -- union every predicted-True row's match list.

    event_keys: the Series returned by match_rows_to_events(df, gt) -- passed
    in rather than recomputed, since it's identical for every rule evaluated
    against the same df/gt and is not cheap to build per-row.
    """
    matched_and_predicted = set()
    for matches in event_keys[predicted_mask]:
        matched_and_predicted.update(matches)
    return len(matched_and_predicted) / len(gt)


def _self_test():
    """Hand-checked sample: exact match, boundary-overlap match (frac==MIN_OVERLAP_FRAC
    exactly), a too-small overlap, a type mismatch, a gt event with zero
    supporting rows, and one row spanning TWO separate gt events (the real
    case this module's docstring describes) -- covers every branch
    match_rows_to_events/event_recall take, including the multi-event-per-row
    case the original one-match-per-row implementation undercounted.
    """
    gt = pd.DataFrame([
        {"sample": "S1", "chr": "chr1", "start": 1000, "end": 2000, "type": "DEL"},  # G1
        {"sample": "S1", "chr": "chr1", "start": 5000, "end": 6000, "type": "DEL"},  # G3
        {"sample": "S2", "chr": "chr2", "start": 500, "end": 600, "type": "DUP"},    # G2, never matched
        {"sample": "S3", "chr": "chr3", "start": 1000, "end": 1200, "type": "DEL"},  # G4
        {"sample": "S3", "chr": "chr3", "start": 5000, "end": 5200, "type": "DEL"},  # G5
    ])
    df = pd.DataFrame([
        {"sample": "S1", "chr": "chr1", "start": 1000, "end": 2000, "type": "DEL"},  # R1: exact -> G1
        {"sample": "S1", "chr": "chr1", "start": 1500, "end": 2500, "type": "DEL"},  # R3: frac==0.5 -> G1
        {"sample": "S1", "chr": "chr1", "start": 1900, "end": 2900, "type": "DEL"},  # R4: frac==0.1 -> no match
        {"sample": "S1", "chr": "chr1", "start": 1000, "end": 2000, "type": "DUP"},  # R5: type mismatch -> no match
        {"sample": "S1", "chr": "chr1", "start": 5000, "end": 6000, "type": "DEL"},  # R6: exact -> G3
        {"sample": "S3", "chr": "chr3", "start": 900,  "end": 5300, "type": "DEL"},  # R7: spans BOTH G4 and G5
    ])
    g1_key = ("S1", "chr1", 1000, 2000, "DEL")
    g3_key = ("S1", "chr1", 5000, 6000, "DEL")
    g4_key = ("S3", "chr3", 1000, 1200, "DEL")
    g5_key = ("S3", "chr3", 5000, 5200, "DEL")

    keys = match_rows_to_events(df, gt)
    expected = [[g1_key], [g1_key], [], [], [g3_key], [g4_key, g5_key]]
    assert list(keys) == expected, f"expected {expected}, got {list(keys)}"

    oracle_recall = event_recall(pd.Series(True, index=df.index), keys, gt)
    # 4 of 5 gt events reachable (G1, G3, G4, G5); G2 (S2) unreachable -> 4/5
    assert oracle_recall == 4 / 5, f"expected 4/5, got {oracle_recall}"

    # The regression check: R7 alone must be credited for BOTH G4 and G5,
    # not just the first one it happens to overlap.
    r7_only_mask = pd.Series([False, False, False, False, False, True], index=df.index)
    r7_recall = event_recall(r7_only_mask, keys, gt)
    assert r7_recall == 2 / 5, f"expected 2/5 (G4 and G5 both credited to R7), got {r7_recall}"

    partial_mask = pd.Series([True, True, False, False, False, False], index=df.index)
    partial_recall = event_recall(partial_mask, keys, gt)
    assert partial_recall == 1 / 5, f"expected 1/5 (only G1 predicted), got {partial_recall}"

    print("_self_test: OK")


if __name__ == "__main__":
    _self_test()

    sys.path.insert(0, str(Path(__file__).parent))
    from load import load_ground_truth, load_merged  # noqa: E402

    df = load_merged()
    gt = load_ground_truth()
    keys = match_rows_to_events(df, gt)
    matched_mask = keys.apply(len) > 0
    matched_events = set()
    for matches in keys:
        matched_events.update(matches)
    print(f"Rows matched to an event: {matched_mask.sum()}/{len(df)}")
    print(f"Distinct events matched by >=1 row: {len(matched_events)}/{len(gt)}")
    print(f"Oracle event recall (predict every row true): {event_recall(pd.Series(True, index=df.index), keys, gt):.4f}")
