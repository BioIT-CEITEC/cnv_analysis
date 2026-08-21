"""
Phase 3 -- EDA & Data Understanding (IMPLEMENTATION_PLAN.md).

Runs the 8 Phase 3 tasks against the validated data from load.py and writes
concrete numeric answers to reports/eda/findings.md. No plots -- the
acceptance criterion is numbers, not generic charts.

Usage:
  python modelling/src/eda.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from load import CALLER_COLUMNS, load_ground_truth, load_merged  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "modules" / "consensus_model"))
from cnv_consensus_model import MIN_OVERLAP_FRAC, overlaps  # noqa: E402

REPORT_PATH = Path(__file__).parent.parent / "reports" / "eda" / "findings.md"


def task1_schema(df):
    """Schema confirmed by load_merged()'s asserts (8 caller columns, present).
    Reported here rather than re-validated."""
    return {
        "n_rows": len(df),
        "n_samples": df["sample"].nunique(),
        "caller_columns": CALLER_COLUMNS,
    }


def task2_label_coverage(df):
    """matches_gt populated on every row (also asserted by load_merged())."""
    return {
        "rows_with_label": int(df["matches_gt"].notna().sum()),
        "total_rows": len(df),
        "fully_populated": bool(df["matches_gt"].notna().all()),
    }


def task3_class_balance(df):
    overall = df["matches_gt"].value_counts().to_dict()
    per_sample_pos = df.groupby("sample")["matches_gt"].sum()
    return {
        "overall_counts": {int(k): int(v) for k, v in overall.items()},
        "overall_ratio_pos": round(df["matches_gt"].mean(), 4),
        "per_sample_positive_min": int(per_sample_pos.min()),
        "per_sample_positive_median": float(per_sample_pos.median()),
        "per_sample_positive_max": int(per_sample_pos.max()),
    }


def task4_per_caller_precision(df):
    rows = []
    for c in CALLER_COLUMNS:
        fired = df[df[c] == 1]
        n_fired = len(fired)
        precision = fired["matches_gt"].mean() if n_fired else float("nan")
        rows.append({"caller": c, "n_fired": n_fired,
                     "precision": round(precision, 4) if n_fired else None})
    return sorted(rows, key=lambda r: (r["precision"] is None, -(r["precision"] or 0)))


def task5_precision_by_agreement(df):
    agreement = df[CALLER_COLUMNS].sum(axis=1)
    out = []
    for k in sorted(agreement.unique()):
        subset = df[agreement == k]
        out.append({
            "agreement_count": int(k),
            "n_rows": len(subset),
            "precision": round(subset["matches_gt"].mean(), 4),
        })
    return out


def task6_blocking_check():
    """Documented, not computed: all_samples_merged.tsv is built purely from
    _merge_candidates() over caller calls (merge_from_bronco.py), so it can only
    contain rows where >=1 caller fired. A true event with zero caller support
    has no row here at all. Therefore the event-level recall denominator MUST
    come from ground_truth_events.tsv, never from counting matches_gt==1 rows.
    """
    return (
        "all_samples_merged.tsv is candidate-only: every row originates from "
        "_merge_candidates() over caller calls, so a true event with zero "
        "caller support cannot appear in this table. Decision: the recall "
        "denominator for all event-level metrics is the row count of "
        "ground_truth_events.tsv, per sample -- never a count derived from "
        "matches_gt rows in the merged table."
    )


def _match_row_to_event(row, gt_events):
    """Return the matching ground-truth event's (sample, chr, start, end, type)
    key, or None. Same matching rule as merge_from_bronco.py's best_gt_match."""
    for g in gt_events:
        if g["chr"] != row["chr"] or g["type"] != row["type"]:
            continue
        if overlaps(row["start"], row["end"], g["start"], g["end"], MIN_OVERLAP_FRAC):
            return (g["chr"], g["start"], g["end"], g["type"])
    return None


def task7_event_structure(df, gt):
    gt_by_sample = {
        sample: sub[["chr", "start", "end", "type"]].to_dict("records")
        for sample, sub in gt.groupby("sample")
    }

    n_gt_events = len(gt)
    events_per_sample = gt.groupby("sample").size()

    positive_rows = df[df["matches_gt"] == 1]
    rows_per_event = {}  # (sample, event_key) -> row count
    for _, row in positive_rows.iterrows():
        events = gt_by_sample.get(row["sample"], [])
        key = _match_row_to_event(row, events)
        if key is not None:
            rows_per_event[(row["sample"], key)] = rows_per_event.get((row["sample"], key), 0) + 1

    matched_events = set(rows_per_event.keys())
    all_event_keys = set()
    for sample, events in gt_by_sample.items():
        for g in events:
            all_event_keys.add((sample, (g["chr"], g["start"], g["end"], g["type"])))
    missed_events = all_event_keys - matched_events

    rows_per_event_counts = pd.Series(list(rows_per_event.values()))
    return {
        "n_gt_events": n_gt_events,
        "n_gt_samples": gt["sample"].nunique(),
        "events_per_sample_mean": round(events_per_sample.mean(), 4),
        "events_per_sample_min": int(events_per_sample.min()),
        "events_per_sample_max": int(events_per_sample.max()),
        "n_events_matched_by_at_least_one_row": len(matched_events),
        "n_events_missed_by_all_callers": len(missed_events),
        "rows_per_event_distribution": (
            rows_per_event_counts.value_counts().sort_index().to_dict()
            if not rows_per_event_counts.empty else {}
        ),
        "rows_per_event_mean": round(rows_per_event_counts.mean(), 4) if not rows_per_event_counts.empty else None,
    }


def task8_span_width(df):
    span = df["end"] - df["start"]
    out = {}
    for label, subset in [("matches_gt=0", span[df["matches_gt"] == 0]),
                           ("matches_gt=1", span[df["matches_gt"] == 1]),
                           ("type=DEL", span[df["type"] == "DEL"]),
                           ("type=DUP", span[df["type"] == "DUP"])]:
        out[label] = {
            "min": int(subset.min()), "p25": float(subset.quantile(0.25)),
            "median": float(subset.median()), "p75": float(subset.quantile(0.75)),
            "max": int(subset.max()), "mean": round(float(subset.mean()), 1),
        }
    return out


def render_markdown(results):
    lines = ["# Phase 3 -- EDA Findings", ""]

    lines += ["## Task 1 -- Schema", ""]
    t1 = results["task1"]
    lines.append(f"- Rows: {t1['n_rows']}, samples: {t1['n_samples']}")
    lines.append(f"- Caller columns (8/8 confirmed present): {', '.join(t1['caller_columns'])}")
    lines.append("")

    lines += ["## Task 2 -- Label coverage", ""]
    t2 = results["task2"]
    lines.append(f"- matches_gt populated on {t2['rows_with_label']}/{t2['total_rows']} rows "
                 f"(fully_populated={t2['fully_populated']})")
    lines.append("")

    lines += ["## Task 3 -- Class balance", ""]
    t3 = results["task3"]
    lines.append(f"- Overall: {t3['overall_counts']} (positive ratio = {t3['overall_ratio_pos']})")
    lines.append(f"- Per-sample positive-row count: min={t3['per_sample_positive_min']}, "
                 f"median={t3['per_sample_positive_median']}, max={t3['per_sample_positive_max']}")
    lines.append("")

    lines += ["## Task 4 -- Per-caller precision `P(matches_gt=1 | caller fired)`", ""]
    lines.append("| caller | n_fired | precision |")
    lines.append("|---|---|---|")
    for r in results["task4"]:
        lines.append(f"| {r['caller']} | {r['n_fired']} | {r['precision']} |")
    lines.append("")

    lines += ["## Task 5 -- Precision by agreement count", ""]
    lines.append("| agreement_count | n_rows | precision |")
    lines.append("|---|---|---|")
    for r in results["task5"]:
        lines.append(f"| {r['agreement_count']} | {r['n_rows']} | {r['precision']} |")
    lines.append("")

    lines += ["## Task 6 -- Blocking check: recall denominator", ""]
    lines.append(results["task6"])
    lines.append("")

    lines += ["## Task 7 -- Event structure", ""]
    t7 = results["task7"]
    lines.append(f"- Ground-truth events: {t7['n_gt_events']} across {t7['n_gt_samples']} samples "
                 f"(mean {t7['events_per_sample_mean']}/sample, range "
                 f"{t7['events_per_sample_min']}-{t7['events_per_sample_max']})")
    lines.append(f"- Events matched by >=1 candidate row: {t7['n_events_matched_by_at_least_one_row']}")
    lines.append(f"- **Events missed by every caller (true false negatives, invisible in "
                 f"all_samples_merged.tsv): {t7['n_events_missed_by_all_callers']}**")
    lines.append(f"- Rows-per-matched-event distribution: {t7['rows_per_event_distribution']} "
                 f"(mean {t7['rows_per_event_mean']})")
    lines.append("")

    lines += ["## Task 8 -- Span width (end - start)", ""]
    lines.append("| group | min | p25 | median | p75 | max | mean |")
    lines.append("|---|---|---|---|---|---|---|")
    for label, s in results["task8"].items():
        lines.append(f"| {label} | {s['min']} | {s['p25']} | {s['median']} | {s['p75']} | {s['max']} | {s['mean']} |")
    lines.append("")

    lines += ["## Leakage / redundancy (for Phase 4)", ""]
    lines.append("- `gt_gene` is populated only when `matches_gt=1` -- confirmed leakage, drop in Phase 4.")
    lines.append("- `callers_supporting` (string) is redundant with the 8 binary flags -- drop the string, "
                 "recompute agreement_count from the flags (already done in Task 5 above).")
    lines.append("")

    return "\n".join(lines)


def main():
    df = load_merged()
    gt = load_ground_truth()

    results = {
        "task1": task1_schema(df),
        "task2": task2_label_coverage(df),
        "task3": task3_class_balance(df),
        "task4": task4_per_caller_precision(df),
        "task5": task5_precision_by_agreement(df),
        "task6": task6_blocking_check(),
        "task7": task7_event_structure(df, gt),
        "task8": task8_span_width(df),
    }

    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(render_markdown(results), encoding="utf-8")
    print(f"Findings written -> {REPORT_PATH}")
    print(f"Rows: {results['task1']['n_rows']}, samples: {results['task1']['n_samples']}")
    print(f"Class balance: {results['task3']['overall_counts']}")
    print(f"GT events missed by every caller: {results['task7']['n_events_missed_by_all_callers']}")


if __name__ == "__main__":
    main()
