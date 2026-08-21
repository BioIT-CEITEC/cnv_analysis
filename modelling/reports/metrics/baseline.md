# Phase 3.5 -- Baseline Gate

F-beta = F2 (recall-weighted, confirmed with user).

## Event-recall ceiling

No row-based rule (or later model) can exceed **0.9839** event-level recall -- 55 of 3414 true events have zero candidate-row support from any caller (Phase 3 Task 7).

## Rule A -- agreement_count >= k

PR-AUC (row-level): 0.7659

| k | n_predicted | row_precision | row_recall | row_F2 | event_recall |
|---|---|---|---|---|---|
| 1 | 55221 | 0.062 | 1.0 | 0.2484 | 0.9839 |
| 2 | 8454 | 0.3908 | 0.965 | 0.7458 | 0.9619 |
| 3 | 3874 | 0.7447 | 0.8426 | 0.821 | 0.845 |
| 4 | 2511 | 0.8467 | 0.6209 | 0.6559 | 0.6227 |
| 5 | 1347 | 0.8931 | 0.3513 | 0.3999 | 0.3524 |
| 6 | 554 | 0.9458 | 0.153 | 0.1839 | 0.1535 |
| 7 | 144 | 0.9653 | 0.0406 | 0.0502 | 0.0407 |
| 8 | 24 | 0.9583 | 0.0067 | 0.0084 | 0.0067 |

## Rule B -- single high-precision caller

| caller | n_predicted | row_precision | row_recall | row_F2 | event_recall |
|---|---|---|---|---|---|
| XHMM | 1682 | 0.8282 | 0.4068 | 0.4529 | 0.408 |
| gatk | 4027 | 0.7842 | 0.9223 | 0.8909 | 0.9227 |

## Rule C -- agreement_count >= k OR XHMM OR gatk

| k | n_predicted | row_precision | row_recall | row_F2 | event_recall |
|---|---|---|---|---|---|
| 2 | 8693 | 0.3808 | 0.9667 | 0.7392 | 0.9625 |
| 3 | 4730 | 0.6808 | 0.9404 | 0.8738 | 0.9402 |
| 4 | 4327 | 0.7338 | 0.9273 | 0.8808 | 0.9274 |
| 5 | 4272 | 0.7416 | 0.9252 | 0.8816 | 0.9253 |

**Best Rule C:** (agreement_count >= 5) OR (XHMM fired) OR (gatk fired) (row_F2=0.8816)

## Overall best baseline

**B (gatk)** -- gatk fired
- row_precision=0.7842, row_recall=0.9223, row_F2=0.8909, event_recall=0.9227, pr_auc=None

**Gate:** Phase 5's model must beat this row_F2 (and ideally the event_recall) materially, or the pre-committed decision is to stop and report that this rule suffices.