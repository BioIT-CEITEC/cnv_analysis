# Phase 5 -- Model Selection & Training

F-beta = F2 (same metric as the Phase 3.5 baseline gate).

## Out-of-fold model comparison

| model | threshold | pr_auc | row_precision | row_recall | row_F2 | event_recall |
|---|---|---|---|---|---|---|
| xgboost (OOF, threshold=0.6992) | 0.6992 | 0.9394 | 0.8715 | 0.9232 | 0.9124 | 0.9233 |
| logreg (OOF, threshold=0.7608) | 0.7608 | 0.9217 | 0.8226 | 0.9279 | 0.9047 | 0.9285 |
| xgboost_no_gatk (OOF, threshold=0.6964) | 0.6964 | 0.9424 | 0.8671 | 0.9258 | 0.9134 | 0.9259 |
| gatk fired | - | - | 0.7842 | 0.9223 | 0.8909 | 0.9227 |

## Gate

**PASSED** -- XGBoost OOF row_F2=0.9124 vs baseline row_F2=0.8909 (gatk fired).

## Feature importance (XGBoost, full-data fit) vs Phase 3 caller precision

Spearman rho, gain-based importance = 0.0476 (p=0.9108)
Spearman rho, permutation importance = -0.2619 (p=0.5309)

| feature | gain importance | permutation importance |
|---|---|---|
| gatk | 0.9227 | 0.4588 |
| conifer | 0.0203 | 0.0212 |
| agreement_count | 0.0133 | 0.0396 |
| panelcnMOPS | 0.0086 | 0.1443 |
| cnvkit | 0.0082 | 0.0128 |
| XHMM | 0.0053 | 0.002 |
| type_DUP | 0.0048 | 0.0161 |
| exomeDepth | 0.0045 | 0.0017 |
| log_span_width | 0.0045 | 0.099 |
| freec | 0.0039 | 0.0097 |
| cnMOPS | 0.0028 | 0.0027 |
| type_conflict | 0.001 | 0.0 |

## Caller redundancy with gatk (on true positives only)

Of the true positives each caller fires on, the fraction gatk also fires on -- high redundancy means the caller adds few TPs gatk doesn't already cover, which is why it gets little marginal importance above.

| caller | n_tp_fired | n_also_gatk | n_unique_tp | redundancy_frac |
|---|---|---|---|---|
| XHMM | 1393 | 1384 | 9 | 0.9935 |
| cnMOPS | 399 | 394 | 5 | 0.9875 |
| panelcnMOPS | 3250 | 3143 | 107 | 0.9671 |
| exomeDepth | 774 | 735 | 39 | 0.9496 |
| cnvkit | 1042 | 988 | 54 | 0.9482 |
| freec | 1449 | 1341 | 108 | 0.9255 |
| conifer | 2163 | 1995 | 168 | 0.9223 |
