# Phase 6 -- Hyperparameter Tuning & Optimization

F-beta = F2 (comparison metric only -- the operating threshold below uses the confirmed precision target instead).

## Best hyperparameters

(4-fold grouped-by-sample CV, seed=123, scored on PR-AUC -- independent of the Phase 4 reporting folds, seed=42)

```json
{
  "subsample": 1.0,
  "scale_pos_weight": 7.563814252336448,
  "n_estimators": 500,
  "min_child_weight": 3,
  "max_depth": 4,
  "learning_rate": 0.05,
  "colsample_bytree": 0.6
}
```

## Top 5 configs from the search (auditability)

| subsample | scale_pos_weight | n_estimators | min_child_weight | max_depth | learning_rate | colsample_bytree | mean_pr_auc |
|---|---|---|---|---|---|---|---|
| 1.0 | 7.563814252336448 | 500 | 3 | 4 | 0.05 | 0.6 | 0.9474 |
| 0.8 | 11.345721378504672 | 500 | 1 | 6 | 0.01 | 0.6 | 0.9471 |
| 0.8 | 11.345721378504672 | 100 | 1 | 5 | 0.2 | 0.8 | 0.9467 |
| 0.6 | 30.255257009345794 | 300 | 5 | 8 | 0.03 | 0.8 | 0.9452 |
| 1.0 | 11.345721378504672 | 100 | 5 | 5 | 0.05 | 0.6 | 0.9449 |

## Honest OOF comparison (Phase 4 persisted folds -- same procedure as Phase 5)

| model | row_precision | row_recall | row_F2 | event_recall | pr_auc |
|---|---|---|---|---|---|
| xgboost_tuned | 0.8813 | 0.9299 | 0.9198 | 0.9291 | 0.9466 |
| xgboost_untuned (Phase 5) | 0.8715 | 0.9232 | 0.9124 | 0.9233 | 0.9394 |
| baseline (gatk fired) | 0.7842 | 0.9223 | 0.8909 | 0.9227 | - |

## Operating threshold (target: precision >= 0.9)

**threshold = 0.7161** -- precision=0.9006, recall=0.9238 at this point on the tuned model's OOF PR curve.

Model persisted -> `/mnt/data/ceitec_cfg2/710000-CEITEC/713000-cmm/713016-bioit/base/workspace/aman/cnv/modelling/models/xgboost_tuned.json`
