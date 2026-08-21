"""
Phase 3.5 -- Baseline gate (IMPLEMENTATION_PLAN.md). Must run before any ML.

Establishes the bar Phase 5's model has to beat "materially" -- if it can't,
the pre-committed decision is to stop and report that a rule suffices.

Usage:
  python modelling/src/baseline.py
"""

import json
import sys
from pathlib import Path

import pandas as pd
from sklearn.metrics import (average_precision_score, fbeta_score,
                             precision_score, recall_score)

sys.path.insert(0, str(Path(__file__).parent))
from events import event_recall, match_rows_to_events  # noqa: E402
from load import CALLER_COLUMNS, load_ground_truth, load_merged  # noqa: E402

BETA = 2  # confirmed with user: favor recall (missing a real CNV costs more
          # than an extra false candidate for manual review)

REPORT_JSON = Path(__file__).parent.parent / "reports" / "metrics" / "baseline.json"
REPORT_MD = Path(__file__).parent.parent / "reports" / "metrics" / "baseline.md"


def evaluate_rule(y_true, predicted_mask, event_keys, gt, description):
    predicted_mask = predicted_mask.astype(bool)
    return {
        "description": description,
        "n_predicted_positive": int(predicted_mask.sum()),
        "row_precision": round(precision_score(y_true, predicted_mask, zero_division=0), 4),
        "row_recall": round(recall_score(y_true, predicted_mask, zero_division=0), 4),
        "row_f2": round(fbeta_score(y_true, predicted_mask, beta=BETA, zero_division=0), 4),
        "event_recall": round(event_recall(predicted_mask, event_keys, gt), 4),
    }


def rule_a_agreement_sweep(df, y_true, agreement, event_keys, gt):
    pr_auc = round(average_precision_score(y_true, agreement), 4)
    sweep = []
    for k in range(1, 9):
        result = evaluate_rule(y_true, agreement >= k, event_keys, gt,
                               f"agreement_count >= {k}")
        result["k"] = k
        sweep.append(result)
    return {"pr_auc": pr_auc, "sweep": sweep}


def rule_b_single_caller(df, y_true, event_keys, gt):
    out = {}
    for caller in ("XHMM", "gatk"):  # XHMM is the named Rule B (highest Task-4
                                     # precision); gatk kept for comparison only
        out[caller] = evaluate_rule(y_true, df[caller] == 1, event_keys, gt,
                                    f"{caller} fired")
        out[caller]["pr_auc"] = None  # single binary decision, no curve
    return out


def rule_c_combination(df, y_true, agreement, event_keys, gt):
    candidates = []
    for k in range(2, 6):
        mask = (agreement >= k) | (df["XHMM"] == 1) | (df["gatk"] == 1)
        result = evaluate_rule(y_true, mask, event_keys, gt,
                               f"(agreement_count >= {k}) OR (XHMM fired) OR (gatk fired)")
        result["k"] = k
        result["pr_auc"] = None
        candidates.append(result)
    best = max(candidates, key=lambda r: r["row_f2"])
    return {"candidates": candidates, "best": best}


def render_markdown(results):
    lines = ["# Phase 3.5 -- Baseline Gate", "", f"F-beta = F{BETA} (recall-weighted, confirmed with user).", ""]

    lines += ["## Event-recall ceiling", ""]
    lines.append(f"No row-based rule (or later model) can exceed **{results['event_recall_ceiling']:.4f}** "
                 f"event-level recall -- {results['n_unreachable_events']} of "
                 f"{results['n_total_events']} true events have zero candidate-row "
                 f"support from any caller (Phase 3 Task 7).")
    lines.append("")

    lines += ["## Rule A -- agreement_count >= k", "", f"PR-AUC (row-level): {results['rule_a']['pr_auc']}", ""]
    lines.append("| k | n_predicted | row_precision | row_recall | row_F2 | event_recall |")
    lines.append("|---|---|---|---|---|---|")
    for r in results["rule_a"]["sweep"]:
        lines.append(f"| {r['k']} | {r['n_predicted_positive']} | {r['row_precision']} | "
                     f"{r['row_recall']} | {r['row_f2']} | {r['event_recall']} |")
    lines.append("")

    lines += ["## Rule B -- single high-precision caller", ""]
    lines.append("| caller | n_predicted | row_precision | row_recall | row_F2 | event_recall |")
    lines.append("|---|---|---|---|---|---|")
    for caller, r in results["rule_b"].items():
        lines.append(f"| {caller} | {r['n_predicted_positive']} | {r['row_precision']} | "
                     f"{r['row_recall']} | {r['row_f2']} | {r['event_recall']} |")
    lines.append("")

    lines += ["## Rule C -- agreement_count >= k OR XHMM OR gatk", ""]
    lines.append("| k | n_predicted | row_precision | row_recall | row_F2 | event_recall |")
    lines.append("|---|---|---|---|---|---|")
    for r in results["rule_c"]["candidates"]:
        lines.append(f"| {r['k']} | {r['n_predicted_positive']} | {r['row_precision']} | "
                     f"{r['row_recall']} | {r['row_f2']} | {r['event_recall']} |")
    lines.append(f"\n**Best Rule C:** {results['rule_c']['best']['description']} "
                 f"(row_F2={results['rule_c']['best']['row_f2']})")
    lines.append("")

    lines += ["## Overall best baseline", ""]
    best = results["best_overall"]
    lines.append(f"**{best['rule']}** -- {best['description']}")
    lines.append(f"- row_precision={best['row_precision']}, row_recall={best['row_recall']}, "
                 f"row_F2={best['row_f2']}, event_recall={best['event_recall']}, "
                 f"pr_auc={best['pr_auc']}")
    lines.append("")
    lines.append("**Gate:** Phase 5's model must beat this row_F2 (and ideally the "
                 "event_recall) materially, or the pre-committed decision is to stop "
                 "and report that this rule suffices.")

    return "\n".join(lines)


def main():
    df = load_merged()
    gt = load_ground_truth()
    y_true = df["matches_gt"]
    agreement = df[CALLER_COLUMNS].sum(axis=1)

    print("Matching rows to ground-truth events (reused across all rules)...")
    event_keys = match_rows_to_events(df, gt)
    all_true_mask = pd.Series(True, index=df.index)
    ceiling = event_recall(all_true_mask, event_keys, gt)
    n_unreachable = len(gt) - int(len(set(event_keys.dropna())))

    print("Rule A: agreement_count sweep...")
    rule_a = rule_a_agreement_sweep(df, y_true, agreement, event_keys, gt)

    print("Rule B: single high-precision caller...")
    rule_b = rule_b_single_caller(df, y_true, event_keys, gt)

    print("Rule C: bounded combination search...")
    rule_c = rule_c_combination(df, y_true, agreement, event_keys, gt)

    candidates_for_best = (
        [{**r, "rule": "A"} for r in rule_a["sweep"]]
        + [{**rule_b["XHMM"], "rule": "B (XHMM)"}]
        + [{**rule_b["gatk"], "rule": "B (gatk)"}]
        + [{**r, "rule": "C"} for r in rule_c["candidates"]]
    )
    best_overall = max(candidates_for_best, key=lambda r: r["row_f2"])

    results = {
        "beta": BETA,
        "n_total_events": len(gt),
        "n_unreachable_events": n_unreachable,
        "event_recall_ceiling": round(ceiling, 4),
        "rule_a": rule_a,
        "rule_b": rule_b,
        "rule_c": rule_c,
        "best_overall": best_overall,
    }

    REPORT_JSON.parent.mkdir(parents=True, exist_ok=True)
    REPORT_JSON.write_text(json.dumps(results, indent=2), encoding="utf-8")
    REPORT_MD.write_text(render_markdown(results), encoding="utf-8")

    print(f"\nWritten -> {REPORT_JSON}")
    print(f"Written -> {REPORT_MD}")
    print(f"\nEvent-recall ceiling: {ceiling:.4f} ({n_unreachable} events unreachable)")
    print(f"Best baseline: {best_overall['description']} "
          f"(row_F2={best_overall['row_f2']}, event_recall={best_overall['event_recall']})")


if __name__ == "__main__":
    main()
