# Demo: TECA feature extraction on synthetic data (no patient data involved)
#
# Generates synthetic 72-hour time series for several biomarkers, computes
# TECA features, and (optionally) runs a tiny XGBoost example to show how the
# extracted features feed into model training.

source("functions/teca_feature_extraction.R")

set.seed(42)

# ----- 1. Thresholds (normal ranges fixed a priori for the study) -----
thresholds <- list(
  temp = c(36.0, 37.3),   # degC
  HR   = c(60, 100),      # beats/min
  RR   = c(12, 20),       # breaths/min
  MAP  = c(70, 105),      # mmHg
  WBC  = c(4.0, 10.0),    # x10^9/L
  CRP  = c(0, 10.0),      # mg/L
  PCT  = c(0, 0.05)       # ng/mL
)

# ----- 2. Synthetic landmark windows -----
# Simulate 200 landmarks; half with a "deteriorating" trajectory (larger
# cumulative excursion beyond thresholds) and half with a stable trajectory.
n_landmarks <- 200
times <- seq(0, 72, by = 6)  # 13 measurement points per 72-h window

make_window <- function(landmark_id, deteriorating) {
  do.call(rbind, lapply(names(thresholds), function(bm) {
    rng <- thresholds[[bm]]
    mid <- mean(rng)
    amp <- diff(rng) / 2
    if (deteriorating) {
      # drift upward beyond the upper threshold toward the end of the window
      drift <- pmax(0, (times - 36) / 36) * amp * runif(1, 1.0, 2.5)
      v <- mid + rnorm(length(times), 0, amp * 0.2) + drift
    } else {
      v <- mid + rnorm(length(times), 0, amp * 0.4)
    }
    data.frame(landmark_id = landmark_id,
               time = times,
               biomarker = bm,
               value = round(v, 3))
  }))
}

sim <- do.call(rbind, lapply(seq_len(n_landmarks), function(i) {
  make_window(i, deteriorating = (i %% 2 == 0))
}))

# ----- 3. TECA feature extraction -----
teca_wide <- teca_features(sim, thresholds, spread = TRUE)
cat("TECA feature matrix:", nrow(teca_wide), "landmarks x",
    ncol(teca_wide) - 1, "features\n")
print(head(teca_wide, 3))

# Synthetic labels: deteriorating trajectories are positive (demo only)
label <- data.frame(landmark_id = seq_len(n_landmarks),
                    y = as.integer(seq_len(n_landmarks) %% 2 == 0))

# ----- 4. (Optional) minimal XGBoost example -----
if (requireNamespace("xgboost", quietly = TRUE)) {
  library(xgboost)
  m <- merge(teca_wide, label, by = "landmark_id")
  X <- as.matrix(m[, setdiff(names(m), c("landmark_id", "y"))])
  y <- m$y
  bst <- xgboost(data = X, label = y, nrounds = 20, objective = "binary:logistic",
                 max_depth = 3, eta = 0.3, verbose = 0)
  pred <- predict(bst, X)
  cat("Demo AUROC (in-sample, synthetic data):",
      round(as.numeric(pROC::roc(y, pred, quiet = TRUE)$auc), 3), "\n")
} else {
  message("xgboost not installed - skipping the modelling demo. ",
          "TECA feature extraction completed successfully.")
}
