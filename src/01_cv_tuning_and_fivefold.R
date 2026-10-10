# Positive rate is too low to justify nesting (the nested-CV script is kept separately); five-fold CV instead
# add regularization
# add the evaluation part
# biomarkers: temp HR RR MAP WBC CRP PCT
# general features: maximum for temp, last available value for the others, mean and SD for all; area = general + _area


# ============================================================================
# Script 1: five-fold cross-validation - split with forced positive distribution (8:2, single round)
# Strategy: distribute the positive patients evenly, assign the negative patients at random, then verify that every fold has at least one positive sample
# keep every patient (subject_id) inside a single fold, and make sure each fold contains positive cases
# ============================================================================


# ===== 1. Library loading =====
library(readxl)
library(dplyr)

# ===== 2. Path configuration =====
DATA_PATH <- Sys.getenv("CLABSI_DATA_PATH", "data/")
SAVE_PATH <- Sys.getenv("CLABSI_OUTPUT_PATH", "output/")
dir.create(SAVE_PATH, showWarnings = FALSE, recursive = TRUE)

# ===== 3. Data loading =====
df <- read_excel(file.path(DATA_PATH, "LM_train_t_vs_lab_clean.xlsx"))
df$clabsi_7d <- as.numeric(df$clabsi_7d)

# ===== 4. Five-fold split with forced positive distribution =====
create_stratified_folds <- function(df, k = 5, seed = 2025) {
  set.seed(seed)
  
  # identify the patient-level label (a patient counts as positive if any landmark is positive)
  patient_labels <- df %>%
    group_by(subject_id) %>%
    summarise(is_positive = max(clabsi_7d, na.rm = TRUE), .groups = "drop")
  
  pos_patients <- patient_labels$subject_id[patient_labels$is_positive == 1]
  neg_patients <- patient_labels$subject_id[patient_labels$is_positive == 0]
  
  cat(sprintf("Total patients: %d (positive %d, negative %d)\n",
              nrow(patient_labels), length(pos_patients), length(neg_patients)))
  
  # positive patients distributed evenly, so that every fold gets some
  pos_folds <- split(sample(pos_patients), rep(1:k, length.out = length(pos_patients)))
  
  # negative patients assigned at random
  neg_folds <- split(sample(neg_patients), sample(1:k, length(neg_patients), replace = TRUE))
  
  # merge the patient IDs belonging to each fold
  all_folds <- lapply(1:k, function(i) c(pos_folds[[i]], neg_folds[[i]]))
  
  # check that every fold has at least one positive sample
  cat("\n===== Forced positive-distribution check =====\n")
  for(i in 1:k) {
    fold_data <- df %>% filter(subject_id %in% all_folds[[i]])
    pos_count <- sum(fold_data$clabsi_7d, na.rm = TRUE)
    cat(sprintf("Fold %d: patients %d | samples %d | positives %d | positive rate %.2f%%\n",
                i, length(all_folds[[i]]), nrow(fold_data), pos_count,
                100*pos_count/nrow(fold_data)))
    
    if(pos_count == 0) stop(sprintf("Fold %d: no positive sample in the validation set; split failed", i))
  }
  
  cat("\nSplit complete: every fold contains at least one positive sample\n")
  return(all_folds)
}

# run the split (8:2, i.e. the validation set holds 20% of the patients)
# note: the p argument of create_folds controls the validation fraction; with 5 folds p = 0.2 gives the 8:2 split
all_folds <- create_stratified_folds(df, k = 5)

# build the training / validation ID lists (8:2)
folds_out_subject_ids <- all_folds  # validation set (20% of the patients)
folds_in_subject_ids <- lapply(1:5, function(i) {
  # training set = the other four folds combined (80% of the patients)
  unique(unlist(all_folds[-i]))
})

# ===== 5. Saving the results =====
save(folds_in_subject_ids, 
     file = file.path(SAVE_PATH, "folds_in_subject_ids_5CV.rds"))
save(folds_out_subject_ids, 
     file = file.path(SAVE_PATH, "folds_out_subject_ids_5CV.rds"))

cat(sprintf("\nFive-fold split complete (80:20); saved to: %s\n", SAVE_PATH))












# XGBoost
# ============================================================================
# Script 2: five-fold cross-validation
# ============================================================================
# use the clean version of the table, with the data types already processed
# switch between the general and area (general + area) models by commenting in the block below; the feature engineering and the saved parameters follow automatically
# SBP and DBP are left out for now, only MAP is used.
# missingness indicators and SpO2 are left out for now. Mean and SD are added for every feature that enters the model
# Bayesian tuning, carried out separately for each fold
# uses the regularization parameters obtained in the preliminary experiments (three sets for general, three for area)

# ===== 1. Library loading =====
library(xgboost)
library(readxl)
library(dplyr)
library(mlrMBO)
library(smoof)
library(lhs)

# ===== 2. Path configuration =====
DATA_PATH <- Sys.getenv("CLABSI_DATA_PATH", "data/")
SAVE_PATH <- Sys.getenv("CLABSI_OUTPUT_PATH", "output/")
dir.create(SAVE_PATH, showWarnings = FALSE, recursive = TRUE)



# ===== 3. Model configuration =====
feature_type <- "AREA"
# feature_type <- "GENERAL"
outcome_col <- "clabsi_7d"
# outcome_col <- "clabsi_3d"
# outcome_col <- "clabsi_1d"



# ===== 4. Data loading =====
df <- read_excel(file.path(DATA_PATH, "LM_train_t_vs_lab_clean.xlsx"))
df$clabsi_7d <- as.numeric(df$clabsi_7d)

load(file.path(SAVE_PATH, "folds_in_subject_ids_5CV.rds"))
load(file.path(SAVE_PATH, "folds_out_subject_ids_5CV.rds"))

# ===== 5. Feature definition =====
BASE_FEATS <- c(
  "lm","gender_numeric", "age", "admission_los", "los_cvc",
  "ventilator", "u_catheter", "OTI_lm", "cvc_location_lm"
)
METRICS <- c("temp", "RR", "HR", "MAP", "WBC", "CRP", "PCT") 
general_value_features <- c(
  "temp_max_general",
  paste0(METRICS[-1], "_last_general")
)
general_stat_features <- unlist(lapply(METRICS, function(m) {
  c(paste0(m, "_mean_general"), paste0(m, "_std_general"))
}))
GENERAL_FEATS <- c(BASE_FEATS, general_value_features, general_stat_features)

if (feature_type == "AREA") {
  area_value_features <- paste0(METRICS, "_area")
  feature_cols <- unique(c(GENERAL_FEATS, area_value_features))
} else {
  feature_cols <- GENERAL_FEATS
}
feature_cols <- intersect(feature_cols, colnames(df))
df[feature_cols] <- df[feature_cols] %>% mutate_all(~replace(., is.na(.), 0))

# ===== Outcome-specific regularization parameters =====
# lookup table: feature_type x outcome_col
reg_params <- list(
  AREA = list(
    "clabsi_1d" = list(lambda = 1.0, gamma = 0.3),
    "clabsi_3d" = list(lambda = 10.5, gamma = 0.2),
    "clabsi_7d" = list(lambda = 5.0, gamma = 0.1)
  ),
  GENERAL = list(
    "clabsi_1d" = list(lambda = 10.0, gamma = 1.0),
    "clabsi_3d" = list(lambda = 5.0, gamma = 0.3),
    "clabsi_7d" = list(lambda = 1.5, gamma = 0.8)
  )
)

# fetch the parameters for the configuration currently selected
current_params <- reg_params[[feature_type]][[outcome_col]]
if (is.null(current_params)) {
  stop(sprintf("Unsupported configuration: feature_type='%s', outcome_col='%s'",
               feature_type, outcome_col))
}

# search range derived from the recommended values of the current configuration
# note: the min(0.4, ...) cap used earlier is dropped so that the upper bound can adapt
lambda_range <- c(max(0.1, current_params$lambda - 2), current_params$lambda + 2)
gamma_range <- c(max(0, current_params$gamma - 0.1), current_params$gamma + 0.1)

# safety check: the lower bound must not be greater than the upper bound
if (lambda_range[1] > lambda_range[2]) lambda_range <- sort(lambda_range)
if (gamma_range[1] > gamma_range[2]) gamma_range <- sort(gamma_range)

cat(sprintf("Recommended parameters for %s-%s: reg_lambda=%.2f, gamma=%.2f\n",
            feature_type, outcome_col, current_params$lambda, current_params$gamma))
cat(sprintf("Search ranges: reg_lambda in [%.2f, %.2f], gamma in [%.2f, %.2f]\n",
            lambda_range[1], lambda_range[2], gamma_range[1], gamma_range[2]))

# ===== 6. Bayesian tuning configuration (regularization parameters included) =====
param.set <- makeParamSet(
  makeNumericParam("eta", lower = 0.01, upper = 0.2, default = 0.1),
  makeIntegerParam("max_depth", lower = 3, upper = 10, default = 6),
  makeIntegerParam("min_child_weight", lower = 5, upper = 80, default = 10),
  makeNumericParam("gamma", lower = gamma_range[1], upper = gamma_range[2], default = current_params$gamma),
  makeNumericParam("reg_lambda", lower = lambda_range[1], upper = lambda_range[2], default = current_params$lambda),
  makeNumericParam("subsample", lower = 0.7, upper = 1.0, default = 0.8),
  makeNumericParam("colsample_bytree", lower = 0.7, upper = 1.0, default = 0.8)
)

set.seed(2024)
design <- generateDesign(n = 10, par.set = param.set, fun = lhs::maximinLHS)

ctrl <- makeMBOControl()
ctrl <- setMBOControlTermination(ctrl, iters = 20)
ctrl <- setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI())

# ===== key fix: predict.type must be set to "se" =====
# wrong: predict.type defaults to "response"
# right: predict.type has to be "se" for the expected-improvement criterion to work
gp.learner <- makeLearner(
  "regr.km", 
  predict.type = "se",  # required by the expected-improvement criterion
  covtype = "matern5_2",
  nugget = 1e-4,
  nugget.stability = 1e-2,
  optim.method = "BFGS"
)

# ===== 7. Objective function factory =====
create_objective_function <- function(fold_idx) {
  train_data <- df %>% filter(subject_id %in% folds_in_subject_ids[[fold_idx]])
  valid_data <- df %>% filter(subject_id %in% folds_out_subject_ids[[fold_idx]])
  
  function(params) {
    # return a penalty value when the validation fold has no positive sample
    if(sum(valid_data[[outcome_col]], na.rm = TRUE) < 1) {
      return(10)
    }
    
    # build the DMatrix
    dtrain <- xgb.DMatrix(
      data = as.matrix(train_data %>% select(all_of(feature_cols))),
      label = train_data[[outcome_col]]
    )
    dvalid <- xgb.DMatrix(
      data = as.matrix(valid_data %>% select(all_of(feature_cols))),
      label = valid_data[[outcome_col]]
    )
    
    # train the model (the regularization parameters are included automatically)
    model <- xgb.train(
      params = c(as.list(params), list(objective = "binary:logistic", eval_metric = "logloss")),
      data = dtrain,
      nrounds = 500,
      watchlist = list(val = dvalid),
      early_stopping_rounds = 25,
      verbose = 0
    )
    
    # predict and compute the log-loss
    preds <- predict(model, dvalid)
    preds <- pmax(pmin(preds, 1 - 1e-6), 1e-6)
    logloss <- mean(-(valid_data[[outcome_col]] * log(preds) + 
                        (1 - valid_data[[outcome_col]]) * log(1 - preds)), 
                    na.rm = TRUE)
    
    # free memory
    rm(model, dtrain, dvalid)
    
    return(ifelse(is.na(logloss) | logloss > 10, 10, logloss))
  }
}

# ===== 8. Running the tuning =====
best_params_list <- list()

for(fold_i in 1:5) {
  cat(sprintf("\n========== Fold %d: Bayesian tuning [%s-%s] ==========",
              fold_i, feature_type, outcome_col))
  
  cat(sprintf("\nTraining: %d samples | %d positives | Validation: %d samples | %d positives\n",
              nrow(df %>% filter(subject_id %in% folds_in_subject_ids[[fold_i]])),
              sum(df$clabsi_7d[df$subject_id %in% folds_in_subject_ids[[fold_i]]]),
              nrow(df %>% filter(subject_id %in% folds_out_subject_ids[[fold_i]])),
              sum(df$clabsi_7d[df$subject_id %in% folds_out_subject_ids[[fold_i]]])))
  
  objective_fn <- create_objective_function(fold_i)
  
  smoof_fn <- makeSingleObjectiveFunction(
    name = sprintf("fold_%d_objective", fold_i),
    fn = objective_fn,
    par.set = param.set,
    minimize = TRUE,
    has.simple.signature = TRUE
  )
  
  # run the Bayesian optimization
  run <- tryCatch({
    mbo(
      fun = smoof_fn,
      control = ctrl,
      show.info = TRUE,
      learner = gp.learner,
      design = design
    )
  }, error = function(e) {
    cat(sprintf("\nFold %d tuning failed: %s\n", fold_i, as.character(e)))
    return(NULL)
  })
  
  # store the best parameters
  if (!is.null(run)) {
    best_params_list[[fold_i]] <- as.list(run$x)
    # print the information that includes the regularization parameters
    cat(sprintf("\nBest parameters: eta=%.4f, max_depth=%d, gamma=%.4f, reg_lambda=%.4f, logloss=%.4f\n",
                run$x$eta, run$x$max_depth, run$x$gamma, run$x$reg_lambda, run$y))
  } else {
    best_params_list[[fold_i]] <- list(
      eta = 0.1,
      max_depth = 6,
      min_child_weight = 10,
      gamma = 0.1,
      reg_lambda = 5.0,  # fallback value
      subsample = 0.8,
      colsample_bytree = 0.8
    )
  }
}

# ===== 9. Saving the results =====
save_file <- file.path(SAVE_PATH, sprintf("best_params_%s_%s.rds", feature_type, outcome_col))
save(best_params_list, file = save_file)

cat(sprintf("\nFive-fold Bayesian tuning complete; parameters saved to: %s\n", save_file))











# ============================================================================
# Script 3: five-fold cross-validation - cost-sensitive evaluation and comparison
# Purpose: evaluate the GENERAL and AREA models with the tuned parameters (AUC, AUPRC, SEN, SPE, PPV, NPV and F1) and bootstrap CIs
# between-model differences for AUROC, AUPRC and Brier score
# saves the models
# ============================================================================

# ===== 1. Library loading =====
library(xgboost)
library(dplyr)
library(tidyr)
library(pROC)
library(PRROC)
library(readxl)
library(openxlsx)
library(purrr)
library(grDevices)
library(ggplot2)
library(gridExtra)

set.seed(42)
rm(list = ls()); gc()

safe_div <- function(num, den) if(den > 0) num/den else 0
outcome_col <- "clabsi_7d"
# outcome_col <- "clabsi_3d"
# outcome_col <- "clabsi_1d"



# ===== 2. Path configuration =====
DATA_PATH <- Sys.getenv("CLABSI_DATA_PATH", "data/")
SAVE_PATH <- Sys.getenv("CLABSI_OUTPUT_PATH", "output/")

# create the directory that will hold the models
MODEL_SAVE_PATH <- file.path(SAVE_PATH, "saved_models", outcome_col)
dir.create(MODEL_SAVE_PATH, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Models will be saved to: %s\n", MODEL_SAVE_PATH))

# ===== 3. Data loading =====
df <- read_excel(file.path(DATA_PATH, "LM_train_t_vs_lab_clean.xlsx"))
df[[outcome_col]] <- as.numeric(df[[outcome_col]])
load(file.path(SAVE_PATH, "folds_out_subject_ids_5CV.rds"))

# ===== 4. Imbalance diagnostics =====
total_pos <- sum(df[[outcome_col]], na.rm = TRUE)
total_neg <- sum(1 - df[[outcome_col]], na.rm = TRUE)
scale_pos_weight <- total_neg / total_pos
custom_weights <- c("clabsi_7d" = 12.45, "clabsi_3d" = 15, "clabsi_1d" = 10)
scale_pos_weight <- custom_weights[outcome_col]

cat("\n===== Class imbalance diagnostics =====\n")
cat(sprintf("Total: %d | positive: %d | negative: %d | positive rate: %.2f%% | class weight: %.1f\n",
            nrow(df), total_pos, total_neg, 100*total_pos/nrow(df), scale_pos_weight))

# ===== 5. Feature definition =====
BASE_FEATS <- c("lm", "gender_numeric", "age", "admission_los", "los_cvc",
                "ventilator", "u_catheter", "OTI_lm", "cvc_location_lm")
METRICS <- c("temp", "RR", "HR", "MAP", "WBC", "CRP", "PCT") 

get_features <- function(type) {
  # mean and SD for every biomarker, plus one value-type summary per biomarker
  # (72-h maximum for body temperature, last available value for the other six)
  stat_feats <- unlist(lapply(METRICS, function(m) {
    c(paste0(m, "_mean_general"), paste0(m, "_std_general"))
  }))
  value_feats <- c("temp_max_general", paste0(METRICS[-1], "_last_general"))
  feats <- c(BASE_FEATS, value_feats, stat_feats)
  if (toupper(type) == "AREA") {
    feats <- c(feats, paste0(METRICS, "_area"))
  }
  selected <- intersect(feats, colnames(df))
  if (length(selected) == 0) stop(sprintf("Feature set %s matched no column", type))
  return(selected)
}

# ===== 6. Loading the tuned parameters =====
load_params_with_weight <- function(param_file, scale_pos_weight) {
  tmp_env <- new.env()
  load(param_file, envir = tmp_env)
  params_list <- tmp_env$best_params_list
  
  lapply(1:5, function(fold_i) {
    fold_params <- as.list(params_list[[fold_i]])
    fold_params$scale_pos_weight <- scale_pos_weight
    fold_params$objective <- "binary:logistic"
    fold_params$eval_metric <- "logloss"
    fold_params
  })
}

# ===== Saving the model =====
save_model_metadata <- function(model, feat_cols, threshold, model_type, fold_idx, 
                                train_data_info, save_path) {
  # save the XGBoost model
  model_file <- file.path(save_path, sprintf("xgb_model_%s_fold%d.json", model_type, fold_idx))
  xgb.save(model, fname = model_file)
  
  # save the model metadata (RDS, contains everything)
  metadata <- list(
    model_type = model_type,
    fold_idx = fold_idx,
    outcome_col = outcome_col,
    feature_names = feat_cols,
    n_features = length(feat_cols),
    optimal_threshold = threshold,
    scale_pos_weight = scale_pos_weight,
    train_data_info = train_data_info,
    saved_time = Sys.time(),
    xgb_model_file = model_file
  )
  
  metadata_file <- file.path(save_path, sprintf("metadata_%s_fold%d.rds", model_type, fold_idx))
  saveRDS(metadata, metadata_file)
  
  # also write the feature list to a plain text file, so that it can be inspected
  feat_file <- file.path(save_path, sprintf("features_%s_fold%d.txt", model_type, fold_idx))
  writeLines(feat_cols, feat_file)
  
  cat(sprintf("  Model saved: %s\n", model_file))
  return(metadata)
}

# ===== 7. Evaluating one fold (model saving added) =====
evaluate_fold <- function(fold_idx, feature_type, params_list, save_model = TRUE) {
  .fold <- fold_idx
  .type <- feature_type
  
  .test_data <- df %>% filter(subject_id %in% folds_out_subject_ids[[.fold]])
  .train_data <- df %>% filter(!subject_id %in% folds_out_subject_ids[[.fold]])
  
  .feat_cols <- get_features(.type)
  .params <- params_list[[.fold]]
  
  cat(sprintf("\n[%s] fold %d/5 | number of features: %d\n", .type, .fold, length(.feat_cols)))
  
  train_start <- Sys.time()
  
  .dtrain <- xgb.DMatrix(
    data = as.matrix(.train_data %>% select(all_of(.feat_cols))),
    label = .train_data[[outcome_col]]
  )
  .dtest <- xgb.DMatrix(
    data = as.matrix(.test_data %>% select(all_of(.feat_cols))),
    label = .test_data[[outcome_col]]
  )
  
  .model <- xgb.train(
    params = .params,
    data = .dtrain, 
    nrounds = 500, 
    verbose = 0, 
    early_stopping_rounds = 50,
    watchlist = list(train = .dtrain, test = .dtest),
    maximize = FALSE
  )
  
  train_end <- Sys.time()
  train_duration <- difftime(train_end, train_start, units = "secs")
  cat(sprintf("  Training time: %.1f s\n", train_duration))
  
  .preds <- predict(.model, .dtest)
  .actual <- .test_data[[outcome_col]]
  
  # discrimination
  .roc_obj <- roc(.actual, .preds, direction = "<", quiet = TRUE)
  .auc <- as.numeric(.roc_obj$auc)
  
  .auprc <- if (sum(.actual == 1) > 0 && sum(.actual == 0) > 0) {
    pr.curve(scores.class0 = .preds[.actual == 1],
             scores.class1 = .preds[.actual == 0], curve = TRUE)$auc.integral
  } else {
    NA_real_
  }
  
  .optimal_threshold <- as.numeric(
    coords(.roc_obj, x = "best", best.method = "youden", transpose = FALSE)$threshold
  )
  
  .pred_class <- as.integer(.preds >= .optimal_threshold)
  
  tn <- sum(.pred_class == 0 & .actual == 0)
  fp <- sum(.pred_class == 1 & .actual == 0)
  fn <- sum(.pred_class == 0 & .actual == 1)
  tp <- sum(.pred_class == 1 & .actual == 1)
  
  .sensitivity <- safe_div(tp, tp + fn)
  .specificity <- safe_div(tn, tn + fp)
  .ppv <- safe_div(tp, tp + fp)
  .npv <- safe_div(tn, tn + fn)
  .f1 <- safe_div(2 * tp, 2 * tp + fp + fn)
  
  cat(sprintf("  AUC: %.3f | AUPRC: %.3f | SEN: %.3f | SPE: %.3f | PPV: %.3f | NPV: %.3f | F1: %.3f\n",
              .auc, .auprc, .sensitivity, .specificity, .ppv, .npv, .f1))
  
  # calibration
  .brier <- mean((.preds - .actual)^2)
  
  # Hosmer-Lemeshow test
  .hl_test <- tryCatch({
    q_vals <- quantile(.preds, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
    q_vals[1] <- 0; q_vals[4] <- max(.preds) + 1e-10
    
    .bin_index <- cut(.preds, breaks = q_vals, include.lowest = TRUE, labels = FALSE)
    
    hl_data <- data.frame(obs = .actual, pred = .preds, bin = .bin_index) %>%
      filter(!is.na(bin)) %>% group_by(bin) %>%
      summarise(n = n(), obs_event = sum(obs), pred_prob = mean(pred), .groups = "drop") %>%
      filter(n > 0)
    
    if(nrow(hl_data) < 2) return(list(statistic = NA_real_, p.value = NA_real_, bins = 0))
    
    hl_stat <- sum(((hl_data$obs_event - hl_data$n * hl_data$pred_prob)^2) / 
                     pmax(hl_data$n * hl_data$pred_prob * (1 - hl_data$pred_prob), 1e-10))
    hl_pval <- pchisq(hl_stat, df = nrow(hl_data) - 2, lower.tail = FALSE)
    
    list(statistic = hl_stat, p.value = hl_pval, bins = nrow(hl_data))
  }, error = function(e) {
    warning("Hosmer-Lemeshow test failed: ", e$message)
    list(statistic = NA_real_, p.value = NA_real_, bins = NA_integer_)
  })
  
  # save the model and its metadata
  if(save_model) {
    train_info <- list(
      n_train = nrow(.train_data),
      n_pos_train = sum(.train_data[[outcome_col]]),
      n_neg_train = sum(1 - .train_data[[outcome_col]]),
      n_test = nrow(.test_data),
      n_pos_test = sum(.actual),
      train_duration_sec = as.numeric(train_duration),
      best_iteration = .model$niter
    )
    
    .metadata <- save_model_metadata(
      model = .model,
      feat_cols = .feat_cols,
      threshold = .optimal_threshold,
      model_type = .type,
      fold_idx = .fold,
      train_data_info = train_info,
      save_path = MODEL_SAVE_PATH
    )
  }
  
  # return the results
  return(tibble(
    fold = .fold,
    auc = .auc, auprc = .auprc, 
    sensitivity = .sensitivity, specificity = .specificity,
    ppv = .ppv, npv = .npv, f1 = .f1,
    threshold = .optimal_threshold,
    n_test = nrow(.test_data), n_pos = sum(.actual),
    actual_values = list(.actual),
    pred_values = list(.preds),
    brier = .brier,
    hl_statistic = .hl_test$statistic,
    hl_pvalue = .hl_test$p.value,
    hl_bins = .hl_test$bins,
    model_file = if(save_model) .metadata$xgb_model_file else NA_character_,
    metadata_file = if(save_model) file.path(MODEL_SAVE_PATH, sprintf("metadata_%s_fold%d.rds", .type, .fold)) else NA_character_
  ))
}

# ===== Training and saving the final model (all data) =====
train_final_model <- function(model_type, params_list) {
  cat(sprintf("\n%s", paste0(rep("=", 60), collapse = "")))
  cat(sprintf("\nTraining the final model: %s (all training data)\n", model_type))
  cat(paste0(rep("=", 60), collapse = ""), "\n")
  
  .feat_cols <- get_features(model_type)
  
  # train on all the data
  .dtrain_full <- xgb.DMatrix(
    data = as.matrix(df %>% select(all_of(.feat_cols))),
    label = df[[outcome_col]]
  )
  
  # use the fold-5 parameters (averaging them would be an option too; fold 5 keeps it simple)
  .params <- params_list[[5]]
  
  train_start <- Sys.time()
  
  .model_final <- xgb.train(
    params = .params,
    data = .dtrain_full, 
    nrounds = 500, 
    verbose = 0, 
    early_stopping_rounds = 50,
    watchlist = list(train = .dtrain_full),
    maximize = FALSE
  )
  
  train_end <- Sys.time()
  train_duration <- difftime(train_end, train_start, units = "secs")
  
  cat(sprintf("  Training time: %.1f s\n", train_duration))
  cat(sprintf("  Boosting rounds in the final model: %d\n", .model_final$niter))
  
  # save the final model
  model_file <- file.path(MODEL_SAVE_PATH, sprintf("xgb_model_%s_final.json", model_type))
  xgb.save(.model_final, fname = model_file)
  
  # compute the threshold on all data, for later use on the test set
  .preds_full <- predict(.model_final, .dtrain_full)
  .roc_full <- roc(df[[outcome_col]], .preds_full, direction = "<", quiet = TRUE)
  .optimal_threshold <- as.numeric(
    coords(.roc_full, x = "best", best.method = "youden", transpose = FALSE)$threshold
  )
  
  # save the metadata
  metadata <- list(
    model_type = model_type,
    fold_idx = "final",
    outcome_col = outcome_col,
    feature_names = .feat_cols,
    n_features = length(.feat_cols),
    optimal_threshold = .optimal_threshold,
    scale_pos_weight = scale_pos_weight,
    train_data_info = list(
      n_train = nrow(df),
      n_pos_train = sum(df[[outcome_col]]),
      n_neg_train = sum(1 - df[[outcome_col]]),
      train_duration_sec = as.numeric(train_duration),
      best_iteration = .model_final$niter
    ),
    saved_time = Sys.time(),
    xgb_model_file = model_file,
    is_final_model = TRUE
  )
  
  metadata_file <- file.path(MODEL_SAVE_PATH, sprintf("metadata_%s_final.rds", model_type))
  saveRDS(metadata, metadata_file)
  
  # write the feature list
  feat_file <- file.path(MODEL_SAVE_PATH, sprintf("features_%s_final.txt", model_type))
  writeLines(.feat_cols, feat_file)
  
  cat(sprintf("  Final model saved: %s\n", model_file))
  cat(sprintf("  Optimal threshold (Youden index): %.4f\n", .optimal_threshold))
  
  return(list(
    model = .model_final,
    metadata = metadata,
    threshold = .optimal_threshold
  ))
}

# ===== 8. Evaluating the two models =====
models_to_evaluate <- c("GENERAL", "AREA")
final_results <- list()
saved_models_info <- list()

total_start <- Sys.time()

for (model_type in models_to_evaluate) {
  cat(sprintf("\n%s", paste0(rep("=", 60), collapse = "")))
  cat(sprintf("\nEvaluating the %s feature set\n", model_type))
  cat(paste0(rep("=", 60), collapse = ""), "\n")
  
  param_file <- file.path(SAVE_PATH, sprintf("best_params_%s_%s.rds", model_type, outcome_col))
  if (!file.exists(param_file)) {
    stop(sprintf("Parameter file not found: %s\nRun the tuning step first.", param_file))
  }
  
  params_list <- load_params_with_weight(param_file, scale_pos_weight)
  
  model_results <- vector("list", 5)
  for (fold_i in 1:5) {
    cat(sprintf("\nFold %d/5...", fold_i))
    model_results[[fold_i]] <- evaluate_fold(fold_i, model_type, params_list, save_model = TRUE)
  }
  
  final_results[[model_type]] <- bind_rows(model_results)
  
  # train and save the final models
  final_model_info <- train_final_model(model_type, params_list)
  saved_models_info[[model_type]] <- final_model_info
  
  rm(model_results, params_list, final_model_info)
  gc(verbose = FALSE)
  
  cat(sprintf("\n%s done. Mean AUROC: %.3f\n",
              model_type,
              mean(final_results[[model_type]]$auc, na.rm = TRUE)))
}

# ===== Saving the model summary =====
models_summary <- list(
  outcome_col = outcome_col,
  scale_pos_weight = scale_pos_weight,
  models = models_to_evaluate,
  model_files = list(
    GENERAL = list(
      folds = file.path(MODEL_SAVE_PATH, sprintf("xgb_model_GENERAL_fold%d.json", 1:5)),
      final = file.path(MODEL_SAVE_PATH, "xgb_model_GENERAL_final.json")
    ),
    AREA = list(
      folds = file.path(MODEL_SAVE_PATH, sprintf("xgb_model_AREA_fold%d.json", 1:5)),
      final = file.path(MODEL_SAVE_PATH, "xgb_model_AREA_final.json")
    )
  ),
  metadata_files = list(
    GENERAL = list(
      folds = file.path(MODEL_SAVE_PATH, sprintf("metadata_GENERAL_fold%d.rds", 1:5)),
      final = file.path(MODEL_SAVE_PATH, "metadata_GENERAL_final.rds")
    ),
    AREA = list(
      folds = file.path(MODEL_SAVE_PATH, sprintf("metadata_AREA_fold%d.rds", 1:5)),
      final = file.path(MODEL_SAVE_PATH, "metadata_AREA_final.rds")
    )
  ),
  created_time = Sys.time(),
  R_version = R.version.string
)

saveRDS(models_summary, file.path(MODEL_SAVE_PATH, "models_summary.rds"))
cat(sprintf("\nModel summary saved: %s\n", file.path(MODEL_SAVE_PATH, "models_summary.rds")))

# ===== 9. Aggregating the results =====
all_general_actual <- do.call(c, final_results[["GENERAL"]]$actual_values)
all_general_preds <- do.call(c, final_results[["GENERAL"]]$pred_values)
all_area_actual <- do.call(c, final_results[["AREA"]]$actual_values)
all_area_preds <- do.call(c, final_results[["AREA"]]$pred_values)

# ===== 10. Bootstrap 95% confidence intervals =====
cat("\n===== Bootstrap confidence intervals (1000 replicates) =====\n")

bootstrap_ci <- function(values) {
  set.seed(42)
  quantile(replicate(1000, mean(sample(values, replace = TRUE))), c(0.025, 0.975))
}

summary_stats <- map_dfr(models_to_evaluate, function(model_type) {
  results <- final_results[[model_type]]
  metrics <- c("auc", "auprc", "sensitivity", "specificity", "ppv", "npv", "f1", "brier")
  
  map_dfr(metrics, function(m) {
    ci <- bootstrap_ci(results[[m]])
    tibble(
      model = model_type,
      metric = recode(m, "auc" = "AUROC", "auprc" = "AUPRC",
                  "sensitivity" = "Sensitivity", "specificity" = "Specificity",
                  "ppv" = "PPV", "npv" = "NPV", "f1" = "F1 score",
                  "brier" = "Brier score"),
      mean = mean(results[[m]], na.rm = TRUE),
      sd = sd(results[[m]], na.rm = TRUE),
      median = median(results[[m]], na.rm = TRUE),
      ci_lower = ci[1],
      ci_upper = ci[2],
      performance = sprintf("%.3f (%.3f-%.3f)", mean(results[[m]], na.rm = TRUE), ci[1], ci[2])
    )
  })
}) %>% arrange(model, metric)

# ===== Bootstrap test for the Brier score difference =====
cat("\n===== Brier score difference test (GENERAL vs AREA) =====\n")

brier_general_vals <- final_results[["GENERAL"]]$brier
brier_area_vals <- final_results[["AREA"]]$brier
observed_diff <- mean(brier_general_vals) - mean(brier_area_vals)

cat(sprintf("Brier score: GENERAL = %.4f +/- %.4f | AREA = %.4f +/- %.4f\n",
            mean(brier_general_vals), sd(brier_general_vals),
            mean(brier_area_vals), sd(brier_area_vals)))
cat(sprintf("Observed difference (GENERAL - AREA): %.4f\n", observed_diff))

# paired bootstrap test
set.seed(42)
n_bootstrap <- 1000
boot_diffs <- numeric(n_bootstrap)

for(i in 1:n_bootstrap) {
  idx <- sample(5, 5, replace = TRUE)
  boot_diffs[i] <- mean(brier_general_vals[idx]) - mean(brier_area_vals[idx])
}

brier_ci_lower <- quantile(boot_diffs, 0.025, names = FALSE)
brier_ci_upper <- quantile(boot_diffs, 0.975, names = FALSE)
brier_p_value <- 2 * min(mean(boot_diffs <= 0), mean(boot_diffs >= 0))
brier_p_value <- max(brier_p_value, 1/n_bootstrap)

cat(sprintf("Bootstrap 95%% CI: %.4f to %.4f\n", brier_ci_lower, brier_ci_upper))
cat(sprintf("p = %.4f | conclusion: %s\n",
            brier_p_value,
            ifelse(brier_p_value < 0.05, 
                   "Brier scores differ significantly between the two models",
                   "no significant difference in Brier score between the two models")))

rm(brier_general_vals, brier_area_vals, boot_diffs, n_bootstrap)
gc(verbose = FALSE)

# ===== Bootstrap test for the AUPRC difference =====
cat("\n===== AUPRC difference test (GENERAL vs AREA) =====\n")

auprc_general_vals <- final_results[["GENERAL"]]$auprc
auprc_area_vals <- final_results[["AREA"]]$auprc
observed_diff_auprc <- mean(auprc_general_vals) - mean(auprc_area_vals)

cat(sprintf("AUPRC: GENERAL = %.4f ± %.4f | AREA = %.4f ± %.4f\n",
            mean(auprc_general_vals), sd(auprc_general_vals),
            mean(auprc_area_vals), sd(auprc_area_vals)))
cat(sprintf("Observed difference (GENERAL - AREA): %.4f\n", observed_diff_auprc))

# paired bootstrap test
set.seed(42)
n_bootstrap <- 1000
boot_diffs_auprc <- numeric(n_bootstrap)

for(i in 1:n_bootstrap) {
  idx <- sample(5, 5, replace = TRUE)
  boot_diffs_auprc[i] <- mean(auprc_general_vals[idx]) - mean(auprc_area_vals[idx])
}

auprc_ci_lower <- quantile(boot_diffs_auprc, 0.025, names = FALSE)
auprc_ci_upper <- quantile(boot_diffs_auprc, 0.975, names = FALSE)
auprc_p_value <- 2 * min(mean(boot_diffs_auprc <= 0), mean(boot_diffs_auprc >= 0))
auprc_p_value <- max(auprc_p_value, 1/n_bootstrap)

cat(sprintf("Bootstrap 95%% CI: %.4f to %.4f\n", auprc_ci_lower, auprc_ci_upper))
cat(sprintf("p = %.4f | conclusion: %s\n",
            auprc_p_value,
            ifelse(auprc_p_value < 0.05, 
                   "AUPRC differs significantly between the two models",
                   "no significant difference in AUPRC between the two models")))

rm(auprc_general_vals, auprc_area_vals, boot_diffs_auprc)
gc(verbose = FALSE)

# ===== 11. DeLong test =====
cat("\n===== Overall DeLong test (GENERAL vs AREA) =====\n")
overall_roc_general <- roc(all_general_actual, all_general_preds, direction = "<", quiet = TRUE)
overall_roc_area <- roc(all_area_actual, all_area_preds, direction = "<", quiet = TRUE)
overall_test <- roc.test(overall_roc_general, overall_roc_area, method = "delong")

cat(sprintf("AUC: GENERAL = %.3f | AREA = %.3f\n", 
            as.numeric(overall_roc_general$auc),
            as.numeric(overall_roc_area$auc)))
cat(sprintf("AUROC difference: %.4f\n",
            as.numeric(overall_roc_general$auc - overall_roc_area$auc)))
cat(sprintf("Z = %.3f | p = %.4f | conclusion: %s\n",
            overall_test$statistic,
            overall_test$p.value,
            ifelse(overall_test$p.value < 0.05, "AUROC differs significantly between the two models", "no significant difference in AUROC between the two models")))

# ===== 12. Writing the Excel report =====
wb <- createWorkbook()

addWorksheet(wb, "Model performance")
writeData(wb, "Model performance", summary_stats, colNames = TRUE, startRow = 1)
setColWidths(wb, "Model performance", cols = 1:ncol(summary_stats), widths = "auto")
freezePane(wb, "Model performance", firstRow = TRUE)

# worksheet: Brier score difference test
addWorksheet(wb, "Brier test")
writeData(wb, "Brier test",
          data.frame(
            test_type = "Paired bootstrap test",
            brier_general = mean(final_results[["GENERAL"]]$brier),
            brier_area = mean(final_results[["AREA"]]$brier),
            brier_diff = observed_diff,
            ci_lower = brier_ci_lower,
            ci_upper = brier_ci_upper,
            p_value = brier_p_value,
            significant = ifelse(brier_p_value < 0.05, "Yes", "No"),
            conclusion = ifelse(brier_p_value < 0.05,
                        "Brier scores differ significantly between the two models",
                        "no significant difference in Brier score between the two models")
          ), startRow = 1)
setColWidths(wb, "Brier test", cols = 1:9, widths = "auto")
freezePane(wb, "Brier test", firstRow = TRUE)

# worksheet: AUPRC difference test
addWorksheet(wb, "AUPRC test")
writeData(wb, "AUPRC test",
          data.frame(
            test_type = "Paired bootstrap test",
            auprc_general = mean(final_results[["GENERAL"]]$auprc),
            auprc_area = mean(final_results[["AREA"]]$auprc),
            auprc_diff = observed_diff_auprc,
            ci_lower = auprc_ci_lower,
            ci_upper = auprc_ci_upper,
            p_value = auprc_p_value,
            significant = ifelse(auprc_p_value < 0.05, "Yes", "No"),
            conclusion = ifelse(auprc_p_value < 0.05,
                        "AUPRC differs significantly between the two models",
                        "no significant difference in AUPRC between the two models")
          ), startRow = 1)
setColWidths(wb, "AUPRC test", cols = 1:9, widths = "auto")
freezePane(wb, "AUPRC test", firstRow = TRUE)

addWorksheet(wb, "DeLong test")
writeData(wb, "DeLong test",
          data.frame(
            test_type = "Overall test",
            auc_general = as.numeric(overall_roc_general$auc),
            auc_area = as.numeric(overall_roc_area$auc),
            auc_diff = as.numeric(overall_roc_general$auc - overall_roc_area$auc),
            z_statistic = overall_test$statistic,
            p_value = overall_test$p.value,
            significant = ifelse(overall_test$p.value < 0.05, "Yes", "No"),
            conclusion = ifelse(overall_test$p.value < 0.05, "AUROC differs significantly between the two models", "no significant difference in AUROC between the two models")
          ), startRow = 1)
setColWidths(wb, "DeLong test", cols = 1:8, widths = "auto")
freezePane(wb, "DeLong test", firstRow = TRUE)

# worksheet: saved model information
model_info_data <- map_dfr(models_to_evaluate, function(model_type) {
  tibble(
    model = model_type,
    outcome = outcome_col,
    n_features = length(get_features(model_type)),
    cv_folds = 5,
    final_model = "Yes",
    model_dir = MODEL_SAVE_PATH,
    model_file = sprintf("xgb_model_%s_final.json", model_type),
    metadata_file = sprintf("metadata_%s_final.rds", model_type),
    feature_list_file = sprintf("features_%s_final.txt", model_type)
  )
})

addWorksheet(wb, "Saved model information")
writeData(wb, "Saved model information", model_info_data, startRow = 1)
setColWidths(wb, "Saved model information", cols = 1:ncol(model_info_data), widths = "auto")
freezePane(wb, "Saved model information", firstRow = TRUE)

# worksheet: notes
addWorksheet(wb, "Notes")
writeData(wb, "Notes", data.frame(
  item = c("Outcome", "Class weight", "Positive rate", "Feature set", "Tuning method", "Threshold rule",
         "DeLong test", "Brier test", "AUPRC test", "Calibration", "Model saving"),
  value = c(
    outcome_col,
    sprintf("%.1f", scale_pos_weight),
    sprintf("%.2f%%", 100*total_pos/nrow(df)),
    "GENERAL and AREA feature sets",
    "Bayesian optimization (tuned independently per fold)",
    "Youden index maximization",
    "AUROC difference between the GENERAL and AREA models (overall test)",
    "Paired bootstrap test on the five-fold cross-validation Brier scores",
    "Paired bootstrap test on the five-fold cross-validation AUPRC",
    "Three bins (forced equal-frequency binning)",
    sprintf("All models saved to: %s", MODEL_SAVE_PATH)
  )
), startRow = 1)

excel_final <- file.path(SAVE_PATH, sprintf("evaluation_results_%s.xlsx", outcome_col))
saveWorkbook(wb, excel_final, overwrite = TRUE)

cat(sprintf("\nDone. Excel report saved to:\n   %s\n", excel_final))


# ===== 13. Generating the figures =====
pdf_file1 <- file.path(SAVE_PATH, sprintf("TRIPOD_assessment_%s.pdf", outcome_col))
cairo_pdf(pdf_file1, width = 20, height = 14)
par(mfrow = c(2, 3), mar = c(3.5, 3.5, 2.5, 1.5), mgp = c(2, 0.8, 0), cex.main = 1.2)

# ROC
plot(overall_roc_general, col = "blue", lwd = 2.5, main = "ROC Curve",
     xlab = "1-Specificity", ylab = "Sensitivity")
lines(overall_roc_area, col = "red", lwd = 2.5)
legend("bottomright", legend = c("GENERAL", "AREA"), 
       col = c("blue", "red"), lwd = 2.5, cex = 0.9, title = "Model")

# PR
pr_general <- pr.curve(scores.class0 = all_general_preds[all_general_actual == 1],
                       scores.class1 = all_general_preds[all_general_actual == 0], curve = TRUE)
pr_area <- pr.curve(scores.class0 = all_area_preds[all_area_actual == 1],
                    scores.class1 = all_area_preds[all_area_actual == 0], curve = TRUE)
plot(pr_general$curve, type = "l", col = "blue", lwd = 2.5, 
     xlim = c(0, 1), ylim = c(0, 1),
     main = "Precision-Recall Curve", xlab = "Recall", ylab = "Precision")
lines(pr_area$curve, col = "red", lwd = 2.5)

dev.off()
cat(sprintf("Assessment figures saved to: %s\n", pdf_file1))

cat("\nAll analyses completed. Output files:\n")
cat("   1. Saved models:\n")
cat(sprintf("      - five-fold cross-validation models: %s\n", MODEL_SAVE_PATH))
cat(sprintf("      - final model (trained on all training data): xgb_model_*_final.json\n"))
cat(sprintf("      - model metadata: metadata_*_*.rds\n"))
cat(sprintf("      - feature lists: features_*_*.txt\n"))
cat(sprintf("      - model summary: models_summary.rds\n"))
cat("   2. AUPRC difference test:\n")
cat("      - paired bootstrap test on the five-fold cross-validation\n")
cat("      - results written to the 'AUPRC test' worksheet\n")
