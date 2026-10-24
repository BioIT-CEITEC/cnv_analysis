# Feature-selection experiment -- drop freec + cnMOPS

F-beta = F2. Dropped features: freec, cnMOPS.

## Five-way comparison (all OOF, Phase 4 persisted folds)

| model | row_precision | row_recall | row_F2 | event_recall | pr_auc |
|---|---|---|---|---|---|
| baseline (gatk fired) | 0.7842 | 0.9223 | 0.8909 | 0.9227 | - |
| xgboost, all features (Phase 5) | 0.8715 | 0.9232 | 0.9124 | 0.9233 | 0.9394 |
| xgboost tuned, all features (Phase 6) | 0.8813 | 0.9299 | 0.9198 | 0.9291 | 0.9466 |
| xgboost, reduced features (no freec/cnMOPS) | 0.8819 | 0.9203 | 0.9123 | 0.9206 | 0.9388 |
| xgboost tuned, reduced features | 0.8835 | 0.9279 | 0.9186 | 0.9277 | 0.9449 |

## Verdict

Dropping freec and cnMOPS changes tuned OOF row_F2 from 0.9198 to 0.9186 (-0.0012) -- **no material difference** compared to keeping all 8 callers.

## Best hyperparameters, reduced feature set

(same search config as Phase 6: 4-fold grouped-by-sample CV, seed=123, scored on PR-AUC)

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

## Operating threshold, reduced model (target: precision >= 0.9)

**threshold = 0.7499** -- precision=0.9002, recall=0.9191.

Model persisted -> `/mnt/data/ceitec_cfg2/710000-CEITEC/713000-cmm/713016-bioit/base/workspace/aman/cnv/modelling/models/xgboost_tuned_no_freec_cnmops.json`
