# Soft-voting ensemble on the hold-out test set

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
library(ggbeeswarm)  # required for SHAP beeswarm plots
library(viridis)     # required for SHAP colours

set.seed(42)
rm(list = ls()); gc()

safe_div <- function(num, den) if(den > 0) num/den else 0

# ===== 2. Configuration =====
# outcome_col <- "clabsi_7d"
# outcome_col <- "clabsi_3d"
outcome_col <- "clabsi_1d"



DATA_PATH <- Sys.getenv("CLABSI_DATA_PATH", "data/")
SAVE_PATH <- Sys.getenv("CLABSI_OUTPUT_PATH", "output/")
MODEL_SAVE_PATH <- file.path(SAVE_PATH, "saved_models", outcome_col)

TEST_FILE <- file.path(DATA_PATH, "LM_test_t_vs_lab_clean.xlsx")

# ===== 3. Feature definition =====
BASE_FEATS <- c("lm","gender_numeric", "age", "admission_los", "los_cvc",
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
  selected <- intersect(feats, colnames(df_test))
  if (length(selected) == 0) stop(sprintf("Feature set %s matched no column", type))
  return(selected)
}

# ===== 4. Calibration functions =====
calculate_calibration_metrics <- function(preds, actual, n_bootstrap = 1000) {
  preds_adj <- pmax(pmin(preds, 0.999), 0.001)
  logit_pred <- log(preds_adj / (1 - preds_adj))
  
  fit <- glm(actual ~ logit_pred, family = binomial(link = "logit"))
  
  coef_summary <- summary(fit)$coefficients
  intercept_est <- coef_summary[1, 1]
  slope_est <- coef_summary[2, 1]
  intercept_se <- coef_summary[1, 2]
  slope_se <- coef_summary[2, 2]
  
  intercept_ci <- c(intercept_est - 1.96 * intercept_se, 
                    intercept_est + 1.96 * intercept_se)
  slope_ci <- c(slope_est - 1.96 * slope_se, 
                slope_est + 1.96 * slope_se)
  
  intercept_z <- intercept_est / intercept_se
  slope_z <- (slope_est - 1) / slope_se
  intercept_p <- 2 * (1 - pnorm(abs(intercept_z)))
  slope_p <- 2 * (1 - pnorm(abs(slope_z)))
  
  set.seed(42)
  boot_intercept <- numeric(n_bootstrap)
  boot_slope <- numeric(n_bootstrap)
  
  for(i in 1:n_bootstrap) {
    idx <- sample(length(preds), replace = TRUE)
    tryCatch({
      boot_fit <- glm(actual[idx] ~ logit_pred[idx], family = binomial(link = "logit"))
      boot_intercept[i] <- coef(boot_fit)[1]
      boot_slope[i] <- coef(boot_fit)[2]
    }, error = function(e) {
      boot_intercept[i] <- NA
      boot_slope[i] <- NA
    })
  }
  
  boot_intercept_ci <- quantile(boot_intercept, c(0.025, 0.975), na.rm = TRUE)
  boot_slope_ci <- quantile(boot_slope, c(0.025, 0.975), na.rm = TRUE)
  
  intercept_ideal <- ifelse(intercept_p >= 0.05, "Yes", "No")
  slope_ideal <- ifelse(slope_p >= 0.05, "Yes", "No")
  
  list(
    intercept = intercept_est,
    intercept_se = intercept_se,
    intercept_ci_lower = intercept_ci[1],
    intercept_ci_upper = intercept_ci[2],
    intercept_p = intercept_p,
    intercept_boot_ci_lower = boot_intercept_ci[1],
    intercept_boot_ci_upper = boot_intercept_ci[2],
    slope = slope_est,
    slope_se = slope_se,
    slope_ci_lower = slope_ci[1],
    slope_ci_upper = slope_ci[2],
    slope_p = slope_p,
    slope_boot_ci_lower = boot_slope_ci[1],
    slope_boot_ci_upper = boot_slope_ci[2],
    intercept_ideal = intercept_ideal,
    slope_ideal = slope_ideal
  )
}

# ===== Brier score confidence interval =====
calculate_brier_ci <- function(preds, actual, n_bootstrap = 1000) {
  set.seed(42)
  n <- length(preds)
  brier_boot <- numeric(n_bootstrap)
  
  for(i in 1:n_bootstrap) {
    idx <- sample(n, replace = TRUE)
    brier_boot[i] <- mean((preds[idx] - actual[idx])^2)
  }
  
  brier_ci <- quantile(brier_boot, c(0.025, 0.975), na.rm = TRUE)
  
  list(
    ci_lower = brier_ci[1],
    ci_upper = brier_ci[2]
  )
}

# ===== AUPRC confidence interval =====
calculate_auprc_ci <- function(preds, actual, n_bootstrap = 1000) {
  set.seed(42)
  n <- length(preds)
  auprc_boot <- numeric(n_bootstrap)
  
  for(i in 1:n_bootstrap) {
    idx <- sample(n, replace = TRUE)
    if (sum(actual[idx] == 1) > 0 && sum(actual[idx] == 0) > 0) {
      auprc_boot[i] <- pr.curve(scores.class0 = preds[idx][actual[idx] == 1],
                                scores.class1 = preds[idx][actual[idx] == 0], 
                                curve = TRUE)$auc.integral
    } else {
      auprc_boot[i] <- NA_real_
    }
  }
  
  auprc_ci <- quantile(auprc_boot, c(0.025, 0.975), na.rm = TRUE)
  
  list(
    ci_lower = auprc_ci[1],
    ci_upper = auprc_ci[2]
  )
}

# ===== Sensitivity and specificity 95% confidence intervals (bootstrap) =====
calculate_sens_spec_ci <- function(preds, actual, threshold, n_bootstrap = 1000) {
  set.seed(42)
  n <- length(preds)
  sens_boot <- numeric(n_bootstrap)
  spec_boot <- numeric(n_bootstrap)
  
  for(i in 1:n_bootstrap) {
    idx <- sample(n, replace = TRUE)
    pred_class <- as.integer(preds[idx] >= threshold)
    
    tp <- sum(pred_class == 1 & actual[idx] == 1)
    fn <- sum(pred_class == 0 & actual[idx] == 1)
    tn <- sum(pred_class == 0 & actual[idx] == 0)
    fp <- sum(pred_class == 1 & actual[idx] == 0)
    
    sens_boot[i] <- safe_div(tp, tp + fn)
    spec_boot[i] <- safe_div(tn, tn + fp)
  }
  
  list(
    sens_ci_lower = quantile(sens_boot, 0.025, na.rm = TRUE),
    sens_ci_upper = quantile(sens_boot, 0.975, na.rm = TRUE),
    spec_ci_lower = quantile(spec_boot, 0.025, na.rm = TRUE),
    spec_ci_upper = quantile(spec_boot, 0.975, na.rm = TRUE)
  )
}

# ===== Brier score difference test =====
test_brier_difference <- function(preds1, preds2, actual, n_bootstrap = 1000) {
  set.seed(42)
  n <- length(actual)
  
  brier1 <- mean((preds1 - actual)^2)
  brier2 <- mean((preds2 - actual)^2)
  obs_diff <- brier1 - brier2
  
  diff_boot <- numeric(n_bootstrap)
  for(i in 1:n_bootstrap) {
    idx <- sample(n, replace = TRUE)
    b1 <- mean((preds1[idx] - actual[idx])^2)
    b2 <- mean((preds2[idx] - actual[idx])^2)
    diff_boot[i] <- b1 - b2
  }
  
  ci <- quantile(diff_boot, c(0.025, 0.975), na.rm = TRUE)
  
  # normal approximation for the p-value, which avoids exact zeros
  se_boot <- sd(diff_boot, na.rm = TRUE)
  z_stat <- obs_diff / se_boot
  p_value <- 2 * (1 - pnorm(abs(z_stat)))
  
  list(
    brier1 = brier1,
    brier2 = brier2,
    difference = obs_diff,
    se_boot = se_boot,
    z_stat = z_stat,
    ci_lower = ci[1],
    ci_upper = ci[2],
    p_value = p_value
  )
}

# ===== AUPRC difference test =====
test_auprc_difference <- function(preds1, preds2, actual, n_bootstrap = 1000) {
  set.seed(42)
  n <- length(actual)
  
  # observed AUPRC
  auprc1 <- if (sum(actual == 1) > 0 && sum(actual == 0) > 0) {
    pr.curve(scores.class0 = preds1[actual == 1],
             scores.class1 = preds1[actual == 0], curve = TRUE)$auc.integral
  } else {
    NA_real_
  }
  
  auprc2 <- if (sum(actual == 1) > 0 && sum(actual == 0) > 0) {
    pr.curve(scores.class0 = preds2[actual == 1],
             scores.class1 = preds2[actual == 0], curve = TRUE)$auc.integral
  } else {
    NA_real_
  }
  
  obs_diff <- auprc1 - auprc2
  
  # paired bootstrap
  diff_boot <- numeric(n_bootstrap)
  for(i in 1:n_bootstrap) {
    idx <- sample(n, replace = TRUE)
    if (sum(actual[idx] == 1) > 0 && sum(actual[idx] == 0) > 0) {
      a1 <- pr.curve(scores.class0 = preds1[idx][actual[idx] == 1],
                     scores.class1 = preds1[idx][actual[idx] == 0], curve = TRUE)$auc.integral
      a2 <- pr.curve(scores.class0 = preds2[idx][actual[idx] == 1],
                     scores.class1 = preds2[idx][actual[idx] == 0], curve = TRUE)$auc.integral
      diff_boot[i] <- a1 - a2
    } else {
      diff_boot[i] <- NA_real_
    }
  }
  
  ci <- quantile(diff_boot, c(0.025, 0.975), na.rm = TRUE)
  se_boot <- sd(diff_boot, na.rm = TRUE)
  z_stat <- obs_diff / se_boot
  p_value <- 2 * (1 - pnorm(abs(z_stat)))
  
  list(
    auprc1 = auprc1,
    auprc2 = auprc2,
    difference = obs_diff,
    se_boot = se_boot,
    z_stat = z_stat,
    ci_lower = ci[1],
    ci_upper = ci[2],
    p_value = p_value
  )
}

# ===== Sensitivity and specificity difference test (bootstrap) =====
test_sens_spec_difference <- function(preds1, preds2, actual, threshold1, threshold2, n_bootstrap = 1000) {
  set.seed(42)
  n <- length(actual)
  
  # observed values
  pred1_class <- as.integer(preds1 >= threshold1)
  pred2_class <- as.integer(preds2 >= threshold2)
  
  tp1 <- sum(pred1_class == 1 & actual == 1)
  fn1 <- sum(pred1_class == 0 & actual == 1)
  tn1 <- sum(pred1_class == 0 & actual == 0)
  fp1 <- sum(pred1_class == 1 & actual == 0)
  
  tp2 <- sum(pred2_class == 1 & actual == 1)
  fn2 <- sum(pred2_class == 0 & actual == 1)
  tn2 <- sum(pred2_class == 0 & actual == 0)
  fp2 <- sum(pred2_class == 1 & actual == 0)
  
  sens1 <- safe_div(tp1, tp1 + fn1)
  sens2 <- safe_div(tp2, tp2 + fn2)
  spec1 <- safe_div(tn1, tn1 + fp1)
  spec2 <- safe_div(tn2, tn2 + fp2)
  
  obs_diff_sens <- sens1 - sens2
  obs_diff_spec <- spec1 - spec2
  
  # bootstrap
  diff_sens_boot <- numeric(n_bootstrap)
  diff_spec_boot <- numeric(n_bootstrap)
  
  for(i in 1:n_bootstrap) {
    idx <- sample(n, replace = TRUE)
    
    p1_class <- as.integer(preds1[idx] >= threshold1)
    p2_class <- as.integer(preds2[idx] >= threshold2)
    
    tp1_b <- sum(p1_class == 1 & actual[idx] == 1)
    fn1_b <- sum(p1_class == 0 & actual[idx] == 1)
    tn1_b <- sum(p1_class == 0 & actual[idx] == 0)
    fp1_b <- sum(p1_class == 1 & actual[idx] == 0)
    
    tp2_b <- sum(p2_class == 1 & actual[idx] == 1)
    fn2_b <- sum(p2_class == 0 & actual[idx] == 1)
    tn2_b <- sum(p2_class == 0 & actual[idx] == 0)
    fp2_b <- sum(p2_class == 1 & actual[idx] == 0)
    
    sens1_b <- safe_div(tp1_b, tp1_b + fn1_b)
    sens2_b <- safe_div(tp2_b, tp2_b + fn2_b)
    spec1_b <- safe_div(tn1_b, tn1_b + fp1_b)
    spec2_b <- safe_div(tn2_b, tn2_b + fp2_b)
    
    diff_sens_boot[i] <- sens1_b - sens2_b
    diff_spec_boot[i] <- spec1_b - spec2_b
  }
  
  ci_sens <- quantile(diff_sens_boot, c(0.025, 0.975), na.rm = TRUE)
  ci_spec <- quantile(diff_spec_boot, c(0.025, 0.975), na.rm = TRUE)
  
  se_sens <- sd(diff_sens_boot, na.rm = TRUE)
  se_spec <- sd(diff_spec_boot, na.rm = TRUE)
  
  z_sens <- obs_diff_sens / se_sens
  z_spec <- obs_diff_spec / se_spec
  
  p_sens <- 2 * (1 - pnorm(abs(z_sens)))
  p_spec <- 2 * (1 - pnorm(abs(z_spec)))
  
  list(
    sens1 = sens1,
    sens2 = sens2,
    spec1 = spec1,
    spec2 = spec2,
    diff_sens = obs_diff_sens,
    diff_spec = obs_diff_spec,
    sens_ci_lower = ci_sens[1],
    sens_ci_upper = ci_sens[2],
    spec_ci_lower = ci_spec[1],
    spec_ci_upper = ci_spec[2],
    z_sens = z_sens,
    z_spec = z_spec,
    p_sens = p_sens,
    p_spec = p_spec
  )
}

# ===== 5. Soft-voting ensemble prediction =====
predict_ensemble_soft <- function(model_type, df_test, feat_cols) {
  cat(sprintf("  Loading the five fold models for soft-voting ensemble...\n"))
  
  preds_matrix <- matrix(NA, nrow = nrow(df_test), ncol = 5)
  
  for(i in 1:5) {
    model_file <- file.path(MODEL_SAVE_PATH, 
                            sprintf("xgb_model_%s_fold%d.json", model_type, i))
    
    if(!file.exists(model_file)) {
      stop(sprintf("Model file not found: %s", model_file))
    }
    
    model <- xgb.load(model_file)
    dtest <- xgb.DMatrix(data = as.matrix(df_test[, feat_cols]))
    preds_matrix[, i] <- predict(model, dtest)
    
    cat(sprintf("    Fold %d: prediction done\n", i))
  }
  
  final_preds <- rowMeans(preds_matrix)
  
  list(
    ensemble_preds = final_preds,
    fold_preds = preds_matrix,
    fold_correlation = cor(preds_matrix, use = "pairwise.complete.obs")
  )
}

# ===== SHAP computation on the test set =====
calculate_shap_ensemble <- function(model_type, df_test, feat_cols) {
  cat(sprintf("  Computing ensemble SHAP values for the %s model...\n", model_type))
  
  # one SHAP matrix per fold
  shap_list <- list()
  feat_values_list <- list()
  
  for(i in 1:5) {
    model_file <- file.path(MODEL_SAVE_PATH, 
                            sprintf("xgb_model_%s_fold%d.json", model_type, i))
    
    if(!file.exists(model_file)) {
      stop(sprintf("Model file not found: %s", model_file))
    }
    
    model <- xgb.load(model_file)
    
    # prepare the data
    dtest <- xgb.DMatrix(data = as.matrix(df_test[, feat_cols]))
    
    # compute the SHAP values
    shap_contrib <- predict(model, dtest, predcontrib = TRUE)
    
    # drop the bias column from the SHAP matrix
    if(ncol(shap_contrib) == length(feat_cols) + 1) {
      shap_mat <- shap_contrib[, -ncol(shap_contrib), drop = FALSE]
    } else if(ncol(shap_contrib) == length(feat_cols)) {
      shap_mat <- shap_contrib
    } else {
      stop("SHAP dimension mismatch")
    }
    
    colnames(shap_mat) <- feat_cols
    shap_list[[i]] <- as.data.frame(shap_mat)
    feat_values_list[[i]] <- df_test[, feat_cols]
    
    cat(sprintf("    Fold %d: SHAP done\n", i))
  }
  
  # average the SHAP values (soft-voting ensemble)
  shap_ensemble <- Reduce("+", shap_list) / length(shap_list)
  
  # mean feature importance
  all_importance <- map_dfr(1:5, function(i) {
    colMeans(abs(shap_list[[i]]))
  })
  mean_importance <- colMeans(all_importance)
  
  list(
    shap_values = shap_ensemble,
    feature_values = df_test[, feat_cols],
    importance = sort(mean_importance, decreasing = TRUE),
    fold_shap = shap_list
  )
}

# ===== 6. Evaluation function =====
evaluate_model <- function(preds, actual, model_name) {
  cat(sprintf("\n  Computing performance metrics for %s...\n", model_name))
  
  roc_obj <- roc(actual, preds, direction = "<", quiet = TRUE)
  auc_val <- as.numeric(roc_obj$auc)
  ci_auc <- ci.auc(roc_obj, method = "bootstrap", boot.n = 1000)
  
  auprc_val <- if (sum(actual == 1) > 0 && sum(actual == 0) > 0) {
    pr.curve(scores.class0 = preds[actual == 1],
             scores.class1 = preds[actual == 0], curve = TRUE)$auc.integral
  } else {
    NA_real_
  }
  
  # AUPRC confidence interval
  auprc_ci <- calculate_auprc_ci(preds, actual)
  
  optimal_threshold <- as.numeric(
    coords(roc_obj, x = "best", best.method = "youden", transpose = FALSE)$threshold
  )
  
  pred_class <- as.integer(preds >= optimal_threshold)
  
  tn <- sum(pred_class == 0 & actual == 0)
  fp <- sum(pred_class == 1 & actual == 0)
  fn <- sum(pred_class == 0 & actual == 1)
  tp <- sum(pred_class == 1 & actual == 1)
  
  sensitivity <- safe_div(tp, tp + fn)
  specificity <- safe_div(tn, tn + fp)
  ppv <- safe_div(tp, tp + fp)
  npv <- safe_div(tn, tn + fn)
  f1 <- safe_div(2 * tp, 2 * tp + fp + fn)
  
  # sensitivity and specificity 95% confidence intervals
  sens_spec_ci <- calculate_sens_spec_ci(preds, actual, optimal_threshold)
  
  brier <- mean((preds - actual)^2)
  brier_ci <- calculate_brier_ci(preds, actual)
  
  eo_ratio <- mean(preds) / mean(actual)
  
  calib_metrics <- calculate_calibration_metrics(preds, actual)
  
  hl_test <- tryCatch({
    q_vals <- quantile(preds, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
    q_vals[1] <- 0; q_vals[4] <- max(preds) + 1e-10
    
    bin_index <- cut(preds, breaks = q_vals, include.lowest = TRUE, labels = FALSE)
    
    hl_data <- data.frame(obs = actual, pred = preds, bin = bin_index) %>%
      filter(!is.na(bin)) %>% 
      group_by(bin) %>%
      summarise(n = n(), obs_event = sum(obs), pred_prob = mean(pred), .groups = "drop") %>%
      filter(n > 0)
    
    if(nrow(hl_data) < 2) {
      list(statistic = NA_real_, p.value = NA_real_)
    } else {
      hl_stat <- sum(((hl_data$obs_event - hl_data$n * hl_data$pred_prob)^2) / 
                       pmax(hl_data$n * hl_data$pred_prob * (1 - hl_data$pred_prob), 1e-10))
      hl_pval <- pchisq(hl_stat, df = nrow(hl_data) - 2, lower.tail = FALSE)
      list(statistic = hl_stat, p.value = hl_pval)
    }
  }, error = function(e) {
    list(statistic = NA_real_, p.value = NA_real_)
  })
  
  cat(sprintf("    AUC: %.3f (%.3f-%.3f) | AUPRC: %.3f (%.3f-%.3f)\n", 
              auc_val, ci_auc[1], ci_auc[3], auprc_val, auprc_ci$ci_lower, auprc_ci$ci_upper))
  cat(sprintf("    SEN: %.3f (%.3f-%.3f) | SPE: %.3f (%.3f-%.3f)\n",
              sensitivity, sens_spec_ci$sens_ci_lower, sens_spec_ci$sens_ci_upper,
              specificity, sens_spec_ci$spec_ci_lower, sens_spec_ci$spec_ci_upper))
  cat(sprintf("    PPV: %.3f | NPV: %.3f | F1: %.3f\n", ppv, npv, f1))
  cat(sprintf("    Calibration: Int=%.3f (p=%.3f), Slope=%.3f (p=%.3f)\n",
              calib_metrics$intercept, calib_metrics$intercept_p,
              calib_metrics$slope, calib_metrics$slope_p))
  cat(sprintf("    Brier: %.4f (%.4f-%.4f)\n",
              brier, brier_ci$ci_lower, brier_ci$ci_upper))
  
  list(
    model_name = model_name,
    n_test = length(actual),
    n_pos = sum(actual),
    auc = auc_val,
    auc_ci_lower = ci_auc[1],
    auc_ci_upper = ci_auc[3],
    auprc = auprc_val,
    auprc_ci_lower = auprc_ci$ci_lower,
    auprc_ci_upper = auprc_ci$ci_upper,
    sensitivity = sensitivity,
    sens_ci_lower = sens_spec_ci$sens_ci_lower,
    sens_ci_upper = sens_spec_ci$sens_ci_upper,
    specificity = specificity,
    spec_ci_lower = sens_spec_ci$spec_ci_lower,
    spec_ci_upper = sens_spec_ci$spec_ci_upper,
    ppv = ppv,
    npv = npv,
    f1 = f1,
    threshold_used = optimal_threshold,
    brier = brier,
    brier_ci_lower = brier_ci$ci_lower,
    brier_ci_upper = brier_ci$ci_upper,
    eo_ratio = eo_ratio,
    hl_statistic = hl_test$statistic,
    hl_pvalue = hl_test$p.value,
    calib_intercept = calib_metrics$intercept,
    calib_intercept_se = calib_metrics$intercept_se,
    calib_intercept_ci_lower = calib_metrics$intercept_ci_lower,
    calib_intercept_ci_upper = calib_metrics$intercept_ci_upper,
    calib_intercept_p = calib_metrics$intercept_p,
    calib_slope = calib_metrics$slope,
    calib_slope_se = calib_metrics$slope_se,
    calib_slope_ci_lower = calib_metrics$slope_ci_lower,
    calib_slope_ci_upper = calib_metrics$slope_ci_upper,
    calib_slope_p = calib_metrics$slope_p,
    calib_intercept_ideal = calib_metrics$intercept_ideal,
    calib_slope_ideal = calib_metrics$slope_ideal,
    preds = preds,
    actual = actual,
    roc_obj = roc_obj
  )
}

# ===== 7. Three-bin calibration table =====
create_3bin_summary <- function(preds, actuals, model_name) {
  q_breaks <- quantile(preds, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  q_breaks[1] <- 0
  q_breaks[4] <- max(preds) + 1e-10
  
  bin_indices <- findInterval(preds, q_breaks, all.inside = TRUE)
  
  results <- lapply(1:3, function(i) {
    idx <- bin_indices == i
    data.frame(
      model = model_name,
      risk_group = factor(c("Low", "Intermediate", "High")[i],
                  levels = c("Low", "Intermediate", "High")),
      threshold_lower = if(i == 1) 0 else q_breaks[i],
      threshold_upper = q_breaks[i + 1],
      mean_predicted = mean(preds[idx]),
      mean_observed = mean(actuals[idx]),
      n = sum(idx),
      calibration_error = mean(preds[idx]) - mean(actuals[idx]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results)
}

# ===== 8. Helper for the detailed metric tables =====
create_detail_table <- function(results) {
  data.frame(
    metric = c("N", "Positives", "Positive rate (%)", "AUROC", "AUROC_CI_lower", "AUROC_CI_upper",
           "AUPRC", "AUPRC_CI_lower", "AUPRC_CI_upper",
           "Sensitivity", "Sensitivity_CI_lower", "Sensitivity_CI_upper",
           "Specificity", "Specificity_CI_lower", "Specificity_CI_upper",
           "PPV", "NPV", "F1", "Optimal threshold",
           "Brier", "Brier_CI_lower", "Brier_CI_upper", "E:O ratio", "HL statistic", "HL p-value",
           "Calibration intercept", "Intercept_CI_lower", "Intercept_CI_upper", "Intercept_p",
           "Calibration slope", "Slope_CI_lower", "Slope_CI_upper", "Slope_p",
           "Intercept ideal", "Slope ideal"),
    value = c(
      results$n_test,
      results$n_pos,
      round(100 * results$n_pos / results$n_test, 2),
      round(results$auc, 3),
      round(results$auc_ci_lower, 3),
      round(results$auc_ci_upper, 3),
      round(results$auprc, 3),
      round(results$auprc_ci_lower, 3),
      round(results$auprc_ci_upper, 3),
      round(results$sensitivity, 3),
      round(results$sens_ci_lower, 3),
      round(results$sens_ci_upper, 3),
      round(results$specificity, 3),
      round(results$spec_ci_lower, 3),
      round(results$spec_ci_upper, 3),
      round(results$ppv, 3),
      round(results$npv, 3),
      round(results$f1, 3),
      round(results$threshold_used, 4),
      round(results$brier, 4),
      round(results$brier_ci_lower, 4),
      round(results$brier_ci_upper, 4),
      round(results$eo_ratio, 3),
      round(results$hl_statistic, 3),
      round(results$hl_pvalue, 4),
      round(results$calib_intercept, 3),
      round(results$calib_intercept_ci_lower, 3),
      round(results$calib_intercept_ci_upper, 3),
      round(results$calib_intercept_p, 4),
      round(results$calib_slope, 3),
      round(results$calib_slope_ci_lower, 3),
      round(results$calib_slope_ci_upper, 3),
      round(results$calib_slope_p, 4),
      results$calib_intercept_ideal,
      results$calib_slope_ideal
    ),
    stringsAsFactors = FALSE
  )
}

# ===== 9. Main program =====
cat("===== Loading the hold-out test set =====\n")
df_test <- read_excel(TEST_FILE)
df_test[[outcome_col]] <- as.numeric(df_test[[outcome_col]])

cat(sprintf("Test set: %d samples | positive: %d (%.2f%%)\n",
            nrow(df_test), 
            sum(df_test[[outcome_col]], na.rm = TRUE),
            100*mean(df_test[[outcome_col]], na.rm = TRUE)))

all_results <- list()
shap_results <- list()  # stores the SHAP results

for(model_type in c("GENERAL", "AREA")) {
  cat(sprintf("\n%s", paste0(rep("=", 70), collapse = "")))
  cat(sprintf("\nSoft-voting ensemble evaluation: %s\n", model_type))
  cat(paste0(rep("=", 70), collapse = ""), "\n")
  
  feat_cols <- get_features(model_type)
  cat(sprintf("Number of features: %d\n", length(feat_cols)))
  
  missing_feats <- setdiff(feat_cols, colnames(df_test))
  if(length(missing_feats) > 0) {
    stop(sprintf("Missing features in the test set: %s", paste(missing_feats, collapse = ", ")))
  }
  
  actual <- df_test[[outcome_col]]
  
  # prediction
  ensemble_pred <- predict_ensemble_soft(model_type, df_test, feat_cols)
  
  # SHAP
  shap_data <- calculate_shap_ensemble(model_type, df_test, feat_cols)
  shap_results[[model_type]] <- shap_data
  
  # evaluation
  results_ensemble <- evaluate_model(
    ensemble_pred$ensemble_preds,
    actual,
    sprintf("%s_Ensemble", model_type)
  )
  
  all_results[[model_type]] <- list(
    ensemble = results_ensemble,
    fold_correlation = ensemble_pred$fold_correlation,
    fold_preds = ensemble_pred$fold_preds
  )
}

# ===== 10. GENERAL vs AREA comparison =====
cat("\n===== GENERAL vs AREA (soft-voting ensemble) =====\n")

comparison_data <- data.frame(
  model = c("GENERAL_ensemble", "AREA_ensemble"),
  AUC = c(all_results[["GENERAL"]]$ensemble$auc,
          all_results[["AREA"]]$ensemble$auc),
  AUC_95CI = c(sprintf("%.3f-%.3f", 
                       all_results[["GENERAL"]]$ensemble$auc_ci_lower,
                       all_results[["GENERAL"]]$ensemble$auc_ci_upper),
               sprintf("%.3f-%.3f",
                       all_results[["AREA"]]$ensemble$auc_ci_lower,
                       all_results[["AREA"]]$ensemble$auc_ci_upper)),
  AUPRC = c(all_results[["GENERAL"]]$ensemble$auprc,
            all_results[["AREA"]]$ensemble$auprc),
  AUPRC_95CI = c(sprintf("%.3f-%.3f",
                         all_results[["GENERAL"]]$ensemble$auprc_ci_lower,
                         all_results[["GENERAL"]]$ensemble$auprc_ci_upper),
                 sprintf("%.3f-%.3f",
                         all_results[["AREA"]]$ensemble$auprc_ci_lower,
                         all_results[["AREA"]]$ensemble$auprc_ci_upper)),
  sensitivity = c(all_results[["GENERAL"]]$ensemble$sensitivity,
          all_results[["AREA"]]$ensemble$sensitivity),
  sensitivity_95CI = c(sprintf("%.3f-%.3f",
                       all_results[["GENERAL"]]$ensemble$sens_ci_lower,
                       all_results[["GENERAL"]]$ensemble$sens_ci_upper),
               sprintf("%.3f-%.3f",
                       all_results[["AREA"]]$ensemble$sens_ci_lower,
                       all_results[["AREA"]]$ensemble$sens_ci_upper)),
  specificity = c(all_results[["GENERAL"]]$ensemble$specificity,
          all_results[["AREA"]]$ensemble$specificity),
  specificity_95CI = c(sprintf("%.3f-%.3f",
                       all_results[["GENERAL"]]$ensemble$spec_ci_lower,
                       all_results[["GENERAL"]]$ensemble$spec_ci_upper),
               sprintf("%.3f-%.3f",
                       all_results[["AREA"]]$ensemble$spec_ci_lower,
                       all_results[["AREA"]]$ensemble$spec_ci_upper)),
  F1 = c(all_results[["GENERAL"]]$ensemble$f1,
         all_results[["AREA"]]$ensemble$f1),
  Brier = c(all_results[["GENERAL"]]$ensemble$brier,
            all_results[["AREA"]]$ensemble$brier),
  Brier_95CI = c(sprintf("%.4f-%.4f",
                         all_results[["GENERAL"]]$ensemble$brier_ci_lower,
                         all_results[["GENERAL"]]$ensemble$brier_ci_upper),
                 sprintf("%.4f-%.4f",
                         all_results[["AREA"]]$ensemble$brier_ci_lower,
                         all_results[["AREA"]]$ensemble$brier_ci_upper)),
  Calibration_Intercept = c(all_results[["GENERAL"]]$ensemble$calib_intercept,
                            all_results[["AREA"]]$ensemble$calib_intercept),
  Calibration_Slope = c(all_results[["GENERAL"]]$ensemble$calib_slope,
                        all_results[["AREA"]]$ensemble$calib_slope),
  stringsAsFactors = FALSE
)

print(comparison_data)

# ===== 11. DeLong test =====
cat("\n===== DeLong test (GENERAL vs AREA) =====\n")

delong_test <- roc.test(
  all_results[["GENERAL"]]$ensemble$roc_obj,
  all_results[["AREA"]]$ensemble$roc_obj,
  method = "delong"
)

cat(sprintf("AUC GENERAL: %.3f (%.3f-%.3f)\n",
            all_results[["GENERAL"]]$ensemble$auc,
            all_results[["GENERAL"]]$ensemble$auc_ci_lower,
            all_results[["GENERAL"]]$ensemble$auc_ci_upper))
cat(sprintf("AUC AREA:    %.3f (%.3f-%.3f)\n",
            all_results[["AREA"]]$ensemble$auc,
            all_results[["AREA"]]$ensemble$auc_ci_lower,
            all_results[["AREA"]]$ensemble$auc_ci_upper))
cat(sprintf("AUROC difference:   %.4f\n",
            all_results[["GENERAL"]]$ensemble$auc - all_results[["AREA"]]$ensemble$auc))
cat(sprintf("Z=%.3f, P=%.4f (%s)\n",
            delong_test$statistic,
            delong_test$p.value,
            ifelse(delong_test$p.value < 0.05, "significant difference", "no significant difference")))

# ===== Brier difference test output =====
cat("\n===== Brier score difference test (GENERAL vs AREA) =====\n")

brier_test <- test_brier_difference(
  all_results[["GENERAL"]]$ensemble$preds,
  all_results[["AREA"]]$ensemble$preds,
  all_results[["GENERAL"]]$ensemble$actual
)

# format the p-value so that it does not print as 0
format_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.4f", p))

cat(sprintf("Brier GENERAL: %.4f (%.4f-%.4f)\n",
            brier_test$brier1,
            all_results[["GENERAL"]]$ensemble$brier_ci_lower,
            all_results[["GENERAL"]]$ensemble$brier_ci_upper))
cat(sprintf("Brier AREA:    %.4f (%.4f-%.4f)\n",
            brier_test$brier2,
            all_results[["AREA"]]$ensemble$brier_ci_lower,
            all_results[["AREA"]]$ensemble$brier_ci_upper))
cat(sprintf("Brier difference:   %.4f (%.4f-%.4f)\n",
            brier_test$difference, brier_test$ci_lower, brier_test$ci_upper))
cat(sprintf("Z=%.3f, P=%s (%s)\n",
            brier_test$z_stat,
            format_p(brier_test$p_value),
            ifelse(brier_test$p_value < 0.05, "significant difference", "no significant difference")))

# ===== AUPRC difference test output =====
cat("\n===== AUPRC difference test (GENERAL vs AREA) =====\n")

auprc_test <- test_auprc_difference(
  all_results[["GENERAL"]]$ensemble$preds,
  all_results[["AREA"]]$ensemble$preds,
  all_results[["GENERAL"]]$ensemble$actual
)

cat(sprintf("AUPRC GENERAL: %.3f (%.3f-%.3f)\n",
            auprc_test$auprc1,
            all_results[["GENERAL"]]$ensemble$auprc_ci_lower,
            all_results[["GENERAL"]]$ensemble$auprc_ci_upper))
cat(sprintf("AUPRC AREA:    %.3f (%.3f-%.3f)\n",
            auprc_test$auprc2,
            all_results[["AREA"]]$ensemble$auprc_ci_lower,
            all_results[["AREA"]]$ensemble$auprc_ci_upper))
cat(sprintf("AUPRC difference:   %.4f (%.3f-%.3f)\n",
            auprc_test$difference, auprc_test$ci_lower, auprc_test$ci_upper))
cat(sprintf("Z=%.3f, P=%s (%s)\n",
            auprc_test$z_stat,
            format_p(auprc_test$p_value),
            ifelse(auprc_test$p_value < 0.05, "significant difference", "no significant difference")))

# ===== Sensitivity and specificity difference test output =====
cat("\n===== Sensitivity and specificity difference test (GENERAL vs AREA) =====\n")

sens_spec_test <- test_sens_spec_difference(
  all_results[["GENERAL"]]$ensemble$preds,
  all_results[["AREA"]]$ensemble$preds,
  all_results[["GENERAL"]]$ensemble$actual,
  all_results[["GENERAL"]]$ensemble$threshold_used,
  all_results[["AREA"]]$ensemble$threshold_used
)

cat(sprintf("Sensitivity GENERAL: %.3f (%.3f-%.3f)\n",
            sens_spec_test$sens1,
            all_results[["GENERAL"]]$ensemble$sens_ci_lower,
            all_results[["GENERAL"]]$ensemble$sens_ci_upper))
cat(sprintf("Sensitivity AREA:    %.3f (%.3f-%.3f)\n",
            sens_spec_test$sens2,
            all_results[["AREA"]]$ensemble$sens_ci_lower,
            all_results[["AREA"]]$ensemble$sens_ci_upper))
cat(sprintf("Sensitivity difference: %.4f (%.4f-%.4f)\n",
            sens_spec_test$diff_sens, sens_spec_test$sens_ci_lower, sens_spec_test$sens_ci_upper))
cat(sprintf("Z=%.3f, P=%s (%s)\n",
            sens_spec_test$z_sens,
            format_p(sens_spec_test$p_sens),
            ifelse(sens_spec_test$p_sens < 0.05, "significant difference", "no significant difference")))

cat("\n")

cat(sprintf("Specificity GENERAL: %.3f (%.3f-%.3f)\n",
            sens_spec_test$spec1,
            all_results[["GENERAL"]]$ensemble$spec_ci_lower,
            all_results[["GENERAL"]]$ensemble$spec_ci_upper))
cat(sprintf("Specificity AREA:    %.3f (%.3f-%.3f)\n",
            sens_spec_test$spec2,
            all_results[["AREA"]]$ensemble$spec_ci_lower,
            all_results[["AREA"]]$ensemble$spec_ci_upper))
cat(sprintf("Specificity difference: %.4f (%.4f-%.4f)\n",
            sens_spec_test$diff_spec, sens_spec_test$spec_ci_lower, sens_spec_test$spec_ci_upper))
cat(sprintf("Z=%.3f, P=%s (%s)\n",
            sens_spec_test$z_spec,
            format_p(sens_spec_test$p_spec),
            ifelse(sens_spec_test$p_spec < 0.05, "significant difference", "no significant difference")))


# ===== 12. Writing the Excel report =====
cat("\n===== Writing the Excel report =====\n")

wb <- createWorkbook()

# sheet 1: model comparison (Brier CI, AUPRC CI and sensitivity/specificity CI)
addWorksheet(wb, "Model comparison")
writeData(wb, "Model comparison", comparison_data, startRow = 1)
setColWidths(wb, "Model comparison", cols = 1:ncol(comparison_data), widths = "auto")
freezePane(wb, "Model comparison", firstRow = TRUE)

# sheet 2: GENERAL detail
general_detail <- create_detail_table(all_results[["GENERAL"]]$ensemble)
addWorksheet(wb, "GENERAL details")
writeData(wb, "GENERAL details", general_detail, startRow = 1)
setColWidths(wb, "GENERAL details", cols = 1:2, widths = "auto")

# sheet 3: AREA detail
area_detail <- create_detail_table(all_results[["AREA"]]$ensemble)
addWorksheet(wb, "AREA details")
writeData(wb, "AREA details", area_detail, startRow = 1)
setColWidths(wb, "AREA details", cols = 1:2, widths = "auto")

# sheet 4: three-bin calibration
calib_3bin_general <- create_3bin_summary(
  all_results[["GENERAL"]]$ensemble$preds,
  all_results[["GENERAL"]]$ensemble$actual,
  "GENERAL"
)
calib_3bin_area <- create_3bin_summary(
  all_results[["AREA"]]$ensemble$preds,
  all_results[["AREA"]]$ensemble$actual,
  "AREA"
)
calib_3bin_all <- rbind(calib_3bin_general, calib_3bin_area)

addWorksheet(wb, "3-bin calibration")
writeData(wb, "3-bin calibration", calib_3bin_all, startRow = 1)
setColWidths(wb, "3-bin calibration", cols = 1:ncol(calib_3bin_all), widths = "auto")
freezePane(wb, "3-bin calibration", firstRow = TRUE)

# sheet 5: DeLong test
delong_table <- data.frame(
  test_item = c("AUROC_GENERAL", "AUROC_AREA", "AUROC difference", "Z statistic", "p-value", "Conclusion"),
  value = c(
    round(all_results[["GENERAL"]]$ensemble$auc, 3),
    round(all_results[["AREA"]]$ensemble$auc, 3),
    round(all_results[["GENERAL"]]$ensemble$auc - all_results[["AREA"]]$ensemble$auc, 4),
    round(delong_test$statistic, 3),
    round(delong_test$p.value, 4),
    ifelse(delong_test$p.value < 0.05, "significant difference", "no significant difference")
  ),
  stringsAsFactors = FALSE
)
addWorksheet(wb, "DeLong test")
writeData(wb, "DeLong test", delong_table, startRow = 1)
setColWidths(wb, "DeLong test", cols = 1:2, widths = "auto")

# sheet 6: Brier difference test
brier_test_table <- data.frame(
  test_item = c("Brier_GENERAL", "Brier_GENERAL_CI_lower", "Brier_GENERAL_CI_upper",
           "Brier_AREA", "Brier_AREA_CI_lower", "Brier_AREA_CI_upper",
           "Brier difference", "Difference_CI_lower", "Difference_CI_upper", "Z statistic", "p-value", "Conclusion"),
  value = c(
    round(brier_test$brier1, 4),
    round(all_results[["GENERAL"]]$ensemble$brier_ci_lower, 4),
    round(all_results[["GENERAL"]]$ensemble$brier_ci_upper, 4),
    round(brier_test$brier2, 4),
    round(all_results[["AREA"]]$ensemble$brier_ci_lower, 4),
    round(all_results[["AREA"]]$ensemble$brier_ci_upper, 4),
    round(brier_test$difference, 4),
    round(brier_test$ci_lower, 4),
    round(brier_test$ci_upper, 4),
    round(brier_test$z_stat, 3),
    format_p(brier_test$p_value),
    ifelse(brier_test$p_value < 0.05, "significant difference", "no significant difference")
  ),
  stringsAsFactors = FALSE
)
addWorksheet(wb, "Brier difference test")
writeData(wb, "Brier difference test", brier_test_table, startRow = 1)
setColWidths(wb, "Brier difference test", cols = 1:2, widths = "auto")

# sheet 7: AUPRC difference test
auprc_test_table <- data.frame(
  test_item = c("AUPRC_GENERAL", "AUPRC_GENERAL_CI_lower", "AUPRC_GENERAL_CI_upper",
           "AUPRC_AREA", "AUPRC_AREA_CI_lower", "AUPRC_AREA_CI_upper",
           "AUPRC difference", "Difference_CI_lower", "Difference_CI_upper", "Z statistic", "p-value", "Conclusion"),
  value = c(
    round(auprc_test$auprc1, 3),
    round(all_results[["GENERAL"]]$ensemble$auprc_ci_lower, 3),
    round(all_results[["GENERAL"]]$ensemble$auprc_ci_upper, 3),
    round(auprc_test$auprc2, 3),
    round(all_results[["AREA"]]$ensemble$auprc_ci_lower, 3),
    round(all_results[["AREA"]]$ensemble$auprc_ci_upper, 3),
    round(auprc_test$difference, 4),
    round(auprc_test$ci_lower, 3),
    round(auprc_test$ci_upper, 3),
    round(auprc_test$z_stat, 3),
    format_p(auprc_test$p_value),
    ifelse(auprc_test$p_value < 0.05, "significant difference", "no significant difference")
  ),
  stringsAsFactors = FALSE
)
addWorksheet(wb, "AUPRC difference test")
writeData(wb, "AUPRC difference test", auprc_test_table, startRow = 1)
setColWidths(wb, "AUPRC difference test", cols = 1:2, widths = "auto")

# sheet 8: sensitivity and specificity difference test
sens_spec_test_table <- data.frame(
  test_item = c("Sensitivity_GENERAL", "Sensitivity_GENERAL_CI_lower", "Sensitivity_GENERAL_CI_upper",
           "Sensitivity_AREA", "Sensitivity_AREA_CI_lower", "Sensitivity_AREA_CI_upper",
           "Sensitivity difference", "Difference_CI_lower", "Difference_CI_upper", "Z statistic", "p-value", "Conclusion",
           "",  # blank separator row
           "Specificity_GENERAL", "Specificity_GENERAL_CI_lower", "Specificity_GENERAL_CI_upper",
           "Specificity_AREA", "Specificity_AREA_CI_lower", "Specificity_AREA_CI_upper",
           "Specificity difference", "Difference_CI_lower", "Difference_CI_upper", "Z statistic", "p-value", "Conclusion"),
  value = c(
    round(sens_spec_test$sens1, 3),
    round(all_results[["GENERAL"]]$ensemble$sens_ci_lower, 3),
    round(all_results[["GENERAL"]]$ensemble$sens_ci_upper, 3),
    round(sens_spec_test$sens2, 3),
    round(all_results[["AREA"]]$ensemble$sens_ci_lower, 3),
    round(all_results[["AREA"]]$ensemble$sens_ci_upper, 3),
    round(sens_spec_test$diff_sens, 4),
    round(sens_spec_test$sens_ci_lower, 4),
    round(sens_spec_test$sens_ci_upper, 4),
    round(sens_spec_test$z_sens, 3),
    format_p(sens_spec_test$p_sens),
    ifelse(sens_spec_test$p_sens < 0.05, "significant difference", "no significant difference"),
    "",  # blank row
    round(sens_spec_test$spec1, 3),
    round(all_results[["GENERAL"]]$ensemble$spec_ci_lower, 3),
    round(all_results[["GENERAL"]]$ensemble$spec_ci_upper, 3),
    round(sens_spec_test$spec2, 3),
    round(all_results[["AREA"]]$ensemble$spec_ci_lower, 3),
    round(all_results[["AREA"]]$ensemble$spec_ci_upper, 3),
    round(sens_spec_test$diff_spec, 4),
    round(sens_spec_test$spec_ci_lower, 4),
    round(sens_spec_test$spec_ci_upper, 4),
    round(sens_spec_test$z_spec, 3),
    format_p(sens_spec_test$p_spec),
    ifelse(sens_spec_test$p_spec < 0.05, "significant difference", "no significant difference")
  ),
  stringsAsFactors = FALSE
)
addWorksheet(wb, "Sensitivity-specificity difference test")
writeData(wb, "Sensitivity-specificity difference test", sens_spec_test_table, startRow = 1)
setColWidths(wb, "Sensitivity-specificity difference test", cols = 1:2, widths = "auto")



# save the Excel report
excel_file <- file.path(SAVE_PATH, sprintf("test_set_evaluation_soft_voting_%s.xlsx", outcome_col))
saveWorkbook(wb, excel_file, overwrite = TRUE)
cat(sprintf("Excel report saved: %s\n", excel_file))



# ===== 13. Generating the figures =====
cat("\n===== Generating figures =====\n")

pdf_file <- file.path(SAVE_PATH, sprintf("test_set_evaluation_soft_voting_%s.pdf", outcome_col))
cairo_pdf(pdf_file, width = 16, height = 8)  # one row, two columns

par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))  # one row, two columns

# 1. ROC curve comparison
plot(all_results[["GENERAL"]]$ensemble$roc_obj, col = "blue", lwd = 2.5,
     main = "ROC curve: GENERAL vs AREA",
     xlab = "1-Specificity", ylab = "Sensitivity")
lines(all_results[["AREA"]]$ensemble$roc_obj, col = "red", lwd = 2.5)
legend("bottomright",
       legend = c(sprintf("GENERAL (%.3f)", all_results[["GENERAL"]]$ensemble$auc),
                  sprintf("AREA (%.3f)", all_results[["AREA"]]$ensemble$auc)),
       col = c("blue", "red"), lwd = 2.5, cex = 0.9)

# 2. PR curve comparison
pr_general <- pr.curve(scores.class0 = all_results[["GENERAL"]]$ensemble$preds[all_results[["GENERAL"]]$ensemble$actual == 1],
                       scores.class1 = all_results[["GENERAL"]]$ensemble$preds[all_results[["GENERAL"]]$ensemble$actual == 0], 
                       curve = TRUE)
pr_area <- pr.curve(scores.class0 = all_results[["AREA"]]$ensemble$preds[all_results[["AREA"]]$ensemble$actual == 1],
                    scores.class1 = all_results[["AREA"]]$ensemble$preds[all_results[["AREA"]]$ensemble$actual == 0], 
                    curve = TRUE)
plot(pr_general$curve, type = "l", col = "blue", lwd = 2.5,
     xlim = c(0, 1), ylim = c(0, 1),
     main = "PR curve: GENERAL vs AREA",
     xlab = "Recall", ylab = "Precision")
lines(pr_area$curve, col = "red", lwd = 2.5)
legend("bottomright",
       legend = c(sprintf("GENERAL (%.3f)", all_results[["GENERAL"]]$ensemble$auprc),
                  sprintf("AREA (%.3f)", all_results[["AREA"]]$ensemble$auprc)),
       col = c("blue", "red"), lwd = 2.5, cex = 0.9)

dev.off()
cat(sprintf("Figures saved: %s\n", pdf_file))






# ===== 14. Manuscript SHAP figures =====

library(patchwork)

# feature-name mapping
map_feature_names <- function(features) {
  mapping <- c(
    "los_cvc" = "Length of catheterization",
    "lm" = "Landmark",
    "admission_los" = "Length of stay",
    "u_catheter" = "Urinary Catheter",
    "ventilator" = "Mechanical Ventilation",
    "cvc_location_lm" = "Subclavian Catheterization",
    "gender_numeric" = "Sex",
    "age" = "Age",
    "OTI_lm" = "Orotracheal Intubation",
    # ---- TECA (area) features, naming aligned with Table 1 ----
    "temp_area" = "Temperature_TECA",
    "PCT_area"  = "PCT_TECA",
    "WBC_area"  = "WBC_TECA",
    "MAP_area"  = "MAP_TECA",
    "HR_area"   = "HR_TECA",
    "RR_area"   = "RR_TECA"
  )
  
  sapply(features, function(x) {
    if (x %in% names(mapping)) {
      return(mapping[x])
    } else {
      return(x)
    }
  })
}

# manuscript-style theme
theme_shap_paper <- function() {
  theme_classic(base_size = 14) +  # base font size raised from 12 to 14
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(color = "black", size = 0.5),
      axis.ticks.x = element_line(color = "black", size = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(hjust = 0, face = "bold", size = 24, margin = margin(b = 10)),  # plot title
      axis.title.x = element_text(face = "bold", size = 20, margin = margin(t = 10)),  # x-axis title
      axis.title.y = element_blank(),
      axis.text.y = element_text(face = "bold", size = 22, color = "black"),  # y-axis tick labels
      axis.text.x = element_text(color = "black", size = 12),  # x-axis tick labels
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12),  # legend title
      legend.text = element_text(size = 11),  # legend text
      plot.margin = margin(15, 5, 15, 15)
    )
}

# prepare the SHAP data, keeping the top N features
prepare_shap_for_plot <- function(shap_results, model_type, n_top = 10) {
  shap_data <- shap_results[[model_type]]
  
  if(is.null(shap_data) || is.null(shap_data$shap_values)) return(NULL)
  
  # top features by mean absolute SHAP value
  top_feats <- names(sort(shap_data$importance, decreasing = TRUE))[1:n_top]
  
  # mean absolute SHAP value, for the bar plot
  importance_df <- data.frame(
    feature = names(shap_data$importance),
    mean_abs_shap = as.numeric(shap_data$importance)
  ) %>%
    filter(feature %in% top_feats) %>%
    arrange(mean_abs_shap)  # ascending order puts the most important feature at the top
  
  # map the feature names
  importance_df$feature_original <- importance_df$feature
  importance_df$feature_label <- map_feature_names(importance_df$feature)
  
  available_feats <- importance_df$feature
  
  if(length(available_feats) == 0) return(NULL)
  
  # extract the data
  shap_mat <- shap_data$shap_values[, available_feats, drop = FALSE]
  feat_mat <- shap_data$feature_values[, available_feats, drop = FALSE]
  
  # long format
  plot_data <- data.frame(
    feature = rep(available_feats, each = nrow(shap_mat)),
    shap_value = as.vector(as.matrix(shap_mat)),
    feature_value = as.vector(as.matrix(feat_mat))
  )
  
  # factor levels, so that the order matches the bar plot
  plot_data$feature <- factor(plot_data$feature, levels = available_feats)
  
  # attach the mapped labels
  level_labels <- setNames(importance_df$feature_label, importance_df$feature)
  plot_data$feature_label <- level_labels[as.character(plot_data$feature)]
  
  # three-band colour scale: bottom 30% blue, middle 40% transition, top 30% red
  plot_data <- plot_data %>%
    group_by(feature) %>%
    mutate(
      # quantile thresholds
      q30 = quantile(feature_value, 0.3, na.rm = TRUE),  # lower colour break
      q70 = quantile(feature_value, 0.7, na.rm = TRUE),  # upper colour break
      min_val = min(feature_value, na.rm = TRUE),
      max_val = max(feature_value, na.rm = TRUE),
      
      # three-band mapping keeps the colour distribution balanced
      feature_value_scaled = case_when(
        feature_value <= q30 ~ 0,  # lowest 30%
        feature_value >= q70 ~ 1,  # highest 30%
        # middle 40% mapped linearly to 0.2-0.8, a neutral transition
        TRUE ~ 0.2 + 0.6 * (feature_value - q30) / (q70 - q30 + 1e-10)
      )
    ) %>%
    ungroup()
  
  return(list(
    plot_data = plot_data,
    importance_df = importance_df,
    top_feats = available_feats,
    level_labels = level_labels
  ))
}

# full SHAP figure: bar plot on the left, beeswarm on the right
create_shap_summary_plot <- function(shap_prep, model_name, outcome_name, 
                                     show_x_label = FALSE, show_legend = FALSE) {
  if(is.null(shap_prep) || is.null(shap_prep$plot_data)) return(NULL)
  
  plot_data <- shap_prep$plot_data
  importance_df <- shap_prep$importance_df
  level_labels <- shap_prep$level_labels
  
  # consistent factor levels: original names, mapped labels
  importance_df$feature <- factor(importance_df$feature, 
                                  levels = importance_df$feature,
                                  labels = importance_df$feature_label)
  plot_data$feature <- factor(plot_data$feature, 
                              levels = names(level_labels),
                              labels = level_labels)
  
  # ========== left panel: bar plot (feature importance) ==========
  p_bar <- ggplot(importance_df, aes(x = mean_abs_shap, y = feature)) +
    geom_bar(stat = "identity", fill = "#1E88E5", width = 0.7) +
    labs(x = NULL, y = NULL, title = model_name) +
    theme_shap_paper() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(face = "bold", size = 22, color = "black"),  # y-axis tick labels
      axis.text.x = element_text(color = "black", size = 11),  # x-axis ticks
      axis.title.x = element_blank(),
      plot.margin = margin(15, 5, 15, 15)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))
  
  # force the x-axis title of the left panel
  if(!show_x_label) {
    p_bar <- p_bar + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    p_bar <- p_bar + 
      labs(x = "mean(|SHAP value|)") +
      theme(axis.title.x = element_text(face = "bold", size = 20, margin = margin(t = 10)))  # x-axis title
  }
  
  # ========== right panel: beeswarm plot ==========
  p_beeswarm <- ggplot(plot_data, aes(x = shap_value, y = feature, color = feature_value_scaled)) +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey60", size = 0.8) +
    ggbeeswarm::geom_quasirandom(
      size = 1.8,  # points slightly larger, 1.5 -> 1.8
      alpha = 0.7, 
      width = 0.3,
      bandwidth = 0.5,
      method = "quasirandom"
    ) +
    # widen the colour transition band
    scale_color_gradient2(
      low = "#2166AC",  # low feature values
      mid = "#F7F7F7",  # neutral middle
      high = "#B2182B",  # high feature values
      midpoint = 0.5,  # midpoint at 0.5
      name = "Feature value",
      breaks = c(0, 1),
      labels = c("Low", "High"),
      limits = c(0, 1),
      guide = guide_colorbar(
        barwidth = 1.5,  # slightly wider
        barheight = unit(6, "cm"),
        ticks = FALSE,
        frame.colour = "black",
        frame.linewidth = 0.5,
        title = "Feature value",
        title.position = "top",
        label.position = "right",
        label.theme = element_text(size = 11),  # legend labels
        draw.ulim = TRUE,
        draw.llim = TRUE
      )
    ) +
    labs(x = NULL, y = NULL) +
    theme_shap_paper() +
    theme(
      axis.text.y = element_blank(),  # the right panel carries no y-axis labels, they are already on the left
      axis.ticks.y = element_blank()
    )
  
  # only the bottom panel shows the x-axis title
  if(!show_x_label) {
    p_beeswarm <- p_beeswarm + 
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
  } else {
    p_beeswarm <- p_beeswarm + 
      labs(x = "SHAP value") +
      theme(axis.title.x = element_text(face = "bold", size = 20, margin = margin(t = 10)))  # x-axis title
  }
  
  # legend placement
  if(!show_legend) {
    p_beeswarm <- p_beeswarm + theme(legend.position = "none")
  }
  
  # combine the panels
  combined_plot <- p_bar + p_beeswarm + 
    plot_layout(widths = c(1, 2))
  
  return(combined_plot)
}

# ========== main ==========

cat("\n===== Generating manuscript SHAP figures (bar plot + beeswarm plot) =====\n")

# prepare the data
shap_gen_prep <- prepare_shap_for_plot(shap_results, "GENERAL", n_top = 10)
shap_area_prep <- prepare_shap_for_plot(shap_results, "AREA", n_top = 10)

# one figure per model; the upper one hides the x-axis title, the lower one shows it
p_general <- create_shap_summary_plot(shap_gen_prep, "Baseline Model", outcome_col, 
                                      show_x_label = FALSE, show_legend = FALSE)
p_area <- create_shap_summary_plot(shap_area_prep, "TECA Model", outcome_col, 
                                   show_x_label = TRUE, show_legend = TRUE)

# stack the two models vertically
if(!is.null(p_general) && !is.null(p_area)) {
  final_plot <- p_general / p_area + 
    plot_layout(heights = c(1, 1))
  
  # PDF output; the canvas was enlarged to make room for the larger fonts
  pdf_file_shap <- file.path(SAVE_PATH, sprintf("test_set_SHAP_analysis_%s.pdf", outcome_col))
  
  pdf(pdf_file_shap, width = 16, height = 14, family = "Times")  # 14x12 -> 16x14, to leave room for the larger fonts
  print(final_plot)
  dev.off()
  
  cat("SHAP figures saved to:", pdf_file_shap, "\n")
} else {
  cat("Failed to generate SHAP figures: data preparation error\n")
}
