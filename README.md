# TECA-CLABSI: Threshold-Exceeding Cumulative Area feature extraction and dynamic CLABSI prediction

Code accompanying the study *"XGBoost-Based Dynamic Prediction of Central Line-Associated Bloodstream Infection Using the Threshold-Exceeding Cumulative Area Feature Extraction Method"*.

This repository contains the **model training and evaluation pipeline** and a self-contained **reference implementation of the TECA feature extraction method**.

It deliberately excludes:

- any patient data;
- the institutional data-cleaning and database-construction scripts.

The input to the training pipeline is a landmark-level analytical table (schema below). TECA feature extraction is provided as a standalone implementation so that it can be applied to any comparable time-series dataset.

## Repository structure

```
├── functions/
│   └── teca_feature_extraction.R      # Reference implementation of the TECA method
├── src/
│   ├── 01_cv_tuning_and_fivefold.R    # Fold construction, Bayesian tuning, five-fold CV
│   └── 02_test_eval_ensemble.R        # Soft-voting ensemble, statistical tests, SHAP
├── demo/
│   └── demo_synthetic.R               # End-to-end demo on synthetic data
└── docs/
    └── RELEASE_NOTES.md               # Scope, thresholds, environment and limitations
```

## The TECA method

For a biomarker with a clinically defined normal range `[lower, upper]`, TECA integrates the **excess amplitude** of measurements beyond the normal range over the feature window (72 h) preceding each landmark:

- inside the range the excess amplitude is 0;
- beyond a boundary it is the distance to that boundary;
- adjacent measurements are integrated with the trapezoidal rule, with linear interpolation at threshold crossings.

Thresholds are fixed a priori from clinical reference ranges and are applied identically to all patients and all landmarks. See `functions/teca_feature_extraction.R` for the exact implementation and `docs/RELEASE_NOTES.md` for the threshold table.

## Requirements

- R >= 4.2
- Modeling: `xgboost`, `mlrMBO`, `mlr`, `smoof`, `ParamHelpers`, `lhs`, `dplyr`, `tidyr`, `purrr`
- Evaluation and reporting: `pROC`, `PRROC`, `readxl`, `openxlsx`, `ggplot2`, `gridExtra`, `patchwork`, `ggbeeswarm`, `viridis`
- Demo (optional): `xgboost`, `pROC`

## Running the pipeline

The scripts read their input and output locations from environment variables, so no path needs to be edited:

```r
Sys.setenv(CLABSI_DATA_PATH   = "path/to/landmark_table/")   # default: data/
Sys.setenv(CLABSI_OUTPUT_PATH = "path/to/output/")           # default: output/
```

`CLABSI_DATA_PATH` must contain the landmark-level table named
`LM_train_t_vs_lab_clean.xlsx` (column schema below). Run the scripts from the
repository root so that the default relative paths (`data/`, `output/`) resolve
correctly.

Then, in order:

```r
source("src/01_cv_tuning_and_fivefold.R")   # step 1: folds, tuning, five-fold CV models
source("src/02_test_eval_ensemble.R")       # step 2: ensemble evaluation on the hold-out set
```

Step 1 must be run before step 2: step 2 loads the fold models written by step 1. Both scripts write every model, parameter set and report under `CLABSI_OUTPUT_PATH`.

## Input data schema

One row per landmark:

- **Outcome columns**: `clabsi_7d`, `clabsi_3d`, `clabsi_1d` — CLABSI within the corresponding prediction window after the landmark;
- **Static features** (9): `lm`, `gender_numeric`, `age`, `admission_los`, `los_cvc`, `ventilator`, `u_catheter`, `OTI_lm`, `cvc_location_lm`;
- **Aggregated dynamic features** (21) for the seven biomarkers `temp`, `RR`, `HR`, `MAP`, `WBC`, `CRP`, `PCT`: `<metric>_mean_general`, `<metric>_std_general` for all seven, `<metric>_last_general` for the six biomarkers other than body temperature, and `temp_max_general` for body temperature;
- **TECA features** (7): one column per biomarker, `<metric>_area`, produced by `teca_features()`.

This gives a baseline feature set of 9 + 21 = 30 features and a TECA feature set of 30 + 7 = 37 features.

Individual patient data are **not** included. De-identified data are available from the corresponding author upon reasonable request and with approval of the institutional ethics committee.

## Modeling protocol

- **Landmarks and prediction windows**: catheter episodes are divided into 24-hour landmarks; the primary outcome is CLABSI within 7 days, with 3-day and 1-day windows as secondary outcomes.
- **Temporal split**: 80% of encounters (by admission, no patient overlap) form the training/validation set and 20% form the hold-out test set. The test set is never used for tuning, threshold selection or feature-set selection.
- **Fold construction**: patient-level stratified five-fold split with a forced positive distribution, so that every fold contains at least one positive sample (`seed = 2025`).
- **Hyperparameter optimization**: model-based optimization (Gaussian-process surrogate with a Matern 5/2 kernel and standard-error prediction, expected-improvement infill) run independently within each fold, on training/validation data only. Seven hyperparameters are tuned (`eta`, `max_depth`, `min_child_weight`, `gamma`, `reg_lambda`, `subsample`, `colsample_bytree`). Budget: a 10-point maximin Latin hypercube initial design followed by 20 expected-improvement iterations (`set.seed(2024)`). The objective is log-loss on the validation fold, with early stopping.
- **Class imbalance**: handled by a fixed positive-class weight (`scale_pos_weight`), selected in preliminary experiments on the training set and then held constant across folds.
- **Boosting**: `binary:logistic` objective with log-loss evaluation, `nrounds = 500`, `early_stopping_rounds = 25`.
- **Hold-out predictions**: soft voting across the five fold models.
- **Statistics**: AUROC differences by the DeLong test; other metrics by bootstrap with 1000 replicates (`set.seed(42)`).

## Reproducing the reported tables

| Reported result | How to reproduce |
|---|---|
| Table 2, cross-validation columns | `src/01_cv_tuning_and_fivefold.R` |
| Table 2, hold-out columns and Figure 2 (SHAP) | `src/02_test_eval_ensemble.R` |
| Table 3 (48 h and 24 h feature windows) | Re-run both scripts with `CLABSI_DATA_PATH` pointing at the landmark table rebuilt with the shortened feature window; no separate script is required. |

## Demo

Run from the repository root:

```r
source("demo/demo_synthetic.R")
```

This generates synthetic time series, extracts TECA features, and (if `xgboost` is installed) fits a minimal model. No real data are used at any point.

## Citation

If you use this code, please cite the manuscript (citation to be added upon publication).

## License

Released under the MIT License — see [`LICENSE`](LICENSE).
