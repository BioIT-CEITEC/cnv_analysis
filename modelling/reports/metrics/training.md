# Phase 5 -- Model Selection & Training

F-beta = F2 (same metric as the Phase 3.5 baseline gate).

## Out-of-fold model comparison

| model | threshold | pr_auc | row_precision | row_recall | row_F2 | event_recall |
|---|---|---|---|---|---|---|
| xgboost (OOF, threshold=0.6992) | 0.6992 | 0.9394 | 0.8715 | 0.9232 | 0.9124 | 0.9233 |
| logreg (OOF, threshold=0.7608) | 0.7608 | 0.9217 | 0.8226 | 0.9279 | 0.9047 | 0.9285 |
| gatk fired | - | - | 0.7842 | 0.9223 | 0.8909 | 0.9227 |

## Gate

**PASSED** -- XGBoost OOF row_F2=0.9124 vs baseline row_F2=0.8909 (gatk fired).

## Feature importance (XGBoost, full-data fit, gain-based) vs Phase 3 caller precision

Spearman rho (caller precision vs importance) = 0.0476 (p=0.9108)

| feature | importance |
|---|---|
| gatk | 0.9227 |
| conifer | 0.0203 |
| agreement_count | 0.0133 |
| panelcnMOPS | 0.0086 |
| cnvkit | 0.0082 |
| XHMM | 0.0053 |
| type_DUP | 0.0048 |
| exomeDepth | 0.0045 |
| log_span_width | 0.0045 |
| freec | 0.0039 |
| cnMOPS | 0.0028 |
| type_conflict | 0.001 |
