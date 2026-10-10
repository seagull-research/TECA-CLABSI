# RELEASE_NOTES

Scope, thresholds, environment and usage notes for this repository.

## Scope

Published:

- `src/01_cv_tuning_and_fivefold.R` — fold construction, per-fold hyperparameter
  optimization, five-fold cross-validation training and evaluation;
- `src/02_test_eval_ensemble.R` — soft-voting ensemble evaluation on the hold-out test
  set, statistical comparisons and SHAP analysis;
- `functions/teca_feature_extraction.R` — standalone reference implementation of the TECA method;
- `demo/demo_synthetic.R` — end-to-end demonstration on synthetic data.

Not published:

- individual patient data;
- the institutional data-cleaning and database-construction scripts. The expected input
  schema is documented in the README instead.

## TECA thresholds

The following normal ranges were fixed a priori from clinical reference ranges, applied
to all patients and all landmarks, and used for the TECA features reported in the paper.
Areas are computed over the 72-hour feature window preceding each landmark (or over the
shortened window in the sensitivity analysis).

| Biomarker | Normal range | Unit |
|---|---|---|
| Body temperature (`temp`) | 36.0 – 37.3 | °C |
| Respiration rate (`RR`) | 12 – 20 | breaths/min |
| Heart rate (`HR`) | 60 – 100 | beats/min |
| Mean arterial pressure (`MAP`) | 70 – 105 | mmHg |
| White blood cell count (`WBC`) | 4 – 10 | ×10⁹/L |
| C-reactive protein (`CRP`) | 0 – 10 | mg/L |
| Procalcitonin (`PCT`) | 0 – 0.05 | ng/mL |

## Environment

The reported analysis used R 4.5.1. Package versions matter for exact reproduction:
`mlrMBO`, `mlr`, `smoof` and `xgboost` in particular should be pinned to the versions
used for the reported results if bit-for-bit reproduction is required.

## Known limitations

- The pipeline expects a pre-computed landmark table; all aggregation columns
  (`<metric>_mean_general`, `<metric>_std_general`, `<metric>_last_general`,
  `temp_max_general`) and the TECA columns (`<metric>_area`) must already be present.
  `functions/teca_feature_extraction.R` computes the TECA columns from raw time series.
- The scripts call `rm(list = ls())` at the start of their main analysis blocks, because
  each block is meant to be run in a fresh session. If you `source()` them into a session
  that holds other objects, run them in their own environment.
- `src/01_cv_tuning_and_fivefold.R` is parameterised at the top of the script:
  `feature_type` (`"AREA"` or `"GENERAL"`) selects the feature set that is tuned, and
  `outcome_col` selects the prediction window. `outcome_col` is selected again at the top
  of the evaluation block, so the two selections must be kept consistent.
