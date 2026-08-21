# Phase 3 -- EDA Findings

## Task 1 -- Schema

- Rows: 55221, samples: 569
- Caller columns (8/8 confirmed present): cnvkit, freec, exomeDepth, panelcnMOPS, cnMOPS, XHMM, conifer, gatk

## Task 2 -- Label coverage

- matches_gt populated on 55221/55221 rows (fully_populated=True)

## Task 3 -- Class balance

- Overall: {0: 51797, 1: 3424} (positive ratio = 0.062)
- Per-sample positive-row count: min=5, median=6.0, max=8

## Task 4 -- Per-caller precision `P(matches_gt=1 | caller fired)`

| caller | n_fired | precision |
|---|---|---|
| XHMM | 1682 | 0.8282 |
| gatk | 4027 | 0.7842 |
| cnMOPS | 864 | 0.4618 |
| exomeDepth | 2465 | 0.314 |
| cnvkit | 3557 | 0.2929 |
| conifer | 9142 | 0.2366 |
| panelcnMOPS | 20169 | 0.1611 |
| freec | 30223 | 0.0479 |

## Task 5 -- Precision by agreement count

| agreement_count | n_rows | precision |
|---|---|---|
| 1 | 46767 | 0.0026 |
| 2 | 4580 | 0.0915 |
| 3 | 1363 | 0.5569 |
| 4 | 1164 | 0.793 |
| 5 | 793 | 0.8562 |
| 6 | 410 | 0.939 |
| 7 | 120 | 0.9667 |
| 8 | 24 | 0.9583 |

## Task 6 -- Blocking check: recall denominator

all_samples_merged.tsv is candidate-only: every row originates from _merge_candidates() over caller calls, so a true event with zero caller support cannot appear in this table. Decision: the recall denominator for all event-level metrics is the row count of ground_truth_events.tsv, per sample -- never a count derived from matches_gt rows in the merged table.

## Task 7 -- Event structure

- Ground-truth events: 3414 across 569 samples (mean 6.0/sample, range 6-6)
- Events matched by >=1 candidate row: 3359
- **Events missed by every caller (true false negatives, invisible in all_samples_merged.tsv): 55**
- Rows-per-matched-event distribution: {1: 3296, 2: 61, 3: 2} (mean 1.0194)

## Task 8 -- Span width (end - start)

| group | min | p25 | median | p75 | max | mean |
|---|---|---|---|---|---|---|
| matches_gt=0 | 13 | 377.0 | 10407.0 | 81516.0 | 155767482 | 4139779.2 |
| matches_gt=1 | 13 | 2880.0 | 11047.0 | 38298.5 | 102558249 | 2059377.2 |
| type=DEL | 13 | 293.0 | 8977.5 | 35795.0 | 152306184 | 2567455.5 |
| type=DUP | 13 | 1126.0 | 15712.0 | 1880057.0 | 155767482 | 6773143.4 |

## Leakage / redundancy (for Phase 4)

- `gt_gene` is populated only when `matches_gt=1` -- confirmed leakage, drop in Phase 4.
- `callers_supporting` (string) is redundant with the 8 binary flags -- drop the string, recompute agreement_count from the flags (already done in Task 5 above).
