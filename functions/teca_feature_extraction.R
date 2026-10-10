# TECA: Threshold-Exceeding Cumulative Area feature extraction
#
# Reference implementation of the TECA method described in:
# "XGBoost-Based Dynamic Prediction of Central Line-Associated Bloodstream
#  Infection Using the Threshold-Exceeding Cumulative Area Feature Extraction
#  Method"
#
# For each dynamic biomarker, TECA quantifies the cumulative excess area
# between the observed time-series measurements and a clinically predefined
# normal range [lower, upper] over a feature window (72 h in the study),
# using trapezoidal integration with linear interpolation at threshold
# crossings.
#
# This file contains no patient data and no study-specific configuration;
# it is a general-purpose, self-contained implementation intended for
# reproducibility of the feature extraction step.

#' Instantaneous excess amplitude relative to the normal range
#'
#' @param v numeric vector of measurements
#' @param lower numeric lower threshold
#' @param upper numeric upper threshold
#' @return excess amplitude: distance beyond the nearest boundary, 0 if inside
excess_amplitude <- function(v, lower, upper) {
  ifelse(v < lower, lower - v,
         ifelse(v > upper, v - upper, 0))
}

#' Area of the excess-amplitude curve over one measurement subinterval
#'
#' Handles all three cases described in the manuscript:
#'   1) both endpoints inside the normal range  -> area = 0
#'   2) both endpoints outside on the same side -> trapezoidal rule
#'   3) segment crossing a threshold boundary   -> linear interpolation
#'      locates the crossing point; trapezoidal integration is applied to
#'      each resulting subsegment
#'
#' @param t1,t2 numeric subinterval start/end times
#' @param x1,x2 numeric measurements at t1 and t2
#' @param lower,upper numeric normal-range boundaries
#' @return numeric area (always >= 0)
segment_excess_area <- function(t1, t2, x1, x2, lower, upper) {
  if (t2 == t1) return(0)
  excess <- function(v) excess_amplitude(v, lower, upper)

  # Build piecewise-linear knots of the excess-amplitude curve over [t1, t2]:
  # start point, any threshold crossings, end point.
  kt <- t1
  ka <- excess(x1)
  for (b in c(lower, upper)) {
    d1 <- x1 - b
    d2 <- x2 - b
    if (d1 * d2 < 0) {                     # sign change -> crossing inside
      t_star <- t1 + (t2 - t1) * (-d1) / (d2 - d1)
      kt <- c(kt, t_star)
      ka <- c(ka, 0)                       # excess is 0 exactly at the boundary
    }
  }
  kt <- c(kt, t2)
  ka <- c(ka, excess(x2))

  o <- order(kt)
  kt <- kt[o]
  ka <- ka[o]
  # Trapezoidal integration over the knots
  sum(diff(kt) * (head(ka, -1) + tail(ka, -1)) / 2)
}

#' TECA for one biomarker over one feature window
#'
#' @param times numeric vector of measurement times (any consistent unit,
#'   e.g. hours); need not be sorted
#' @param values numeric vector of measurements, same length as \code{times}
#' @param lower,upper numeric normal-range boundaries
#' @return numeric scalar: cumulative threshold-exceeding area
#' @examples
#' # Body temperature measured 6-hourly over a 72-h window, normal 36.0-37.3 C
#' t  <- seq(0, 72, by = 6)
#' x  <- c(36.8, 37.0, 37.2, 37.8, 38.4, 38.9, 38.1, 37.6, 37.9, 38.2, 38.6, 38.0, 37.4)
#' teca(times = t, values = x, lower = 36.0, upper = 37.3)
teca <- function(times, values, lower, upper) {
  stopifnot(length(times) == length(values))
  stopifnot(lower <= upper)
  if (any(is.na(values))) {
    keep <- !is.na(values)
    times <- times[keep]
    values <- values[keep]
  }
  n <- length(values)
  if (n == 0) return(NA_real_)
  if (n == 1) return(excess_amplitude(values[1], lower, upper))

  o <- order(times)
  t <- times[o]
  x <- values[o]

  total <- 0
  for (i in seq_len(n - 1)) {
    total <- total + segment_excess_area(t[i], t[i + 1], x[i], x[i + 1],
                                         lower, upper)
  }
  total
}

#' TECA features for a landmark-level time-series table
#'
#' Computes one TECA value per (landmark, biomarker) pair.
#'
#' @param df data.frame with columns:
#'   \itemize{
#'     \item \code{landmark_id}: identifier of the landmark (one feature window)
#'     \item \code{time}: measurement time within the window
#'     \item \code{biomarker}: biomarker name
#'     \item \code{value}: measured value
#'   }
#' @param thresholds named list of \code{c(lower, upper)} per biomarker,
#'   e.g. \code{list(temp = c(36.0, 37.3), HR = c(60, 100))}
#' @return data.frame with columns \code{landmark_id}, \code{biomarker},
#'   \code{teca} (wide format if \code{spread = TRUE}: one column per biomarker,
#'   named \code{<biomarker>_area}, matching the training pipeline input)
#' @export
teca_features <- function(df, thresholds, spread = TRUE) {
  stopifnot(all(c("landmark_id", "time", "biomarker", "value") %in% names(df)))
  stopifnot(all(unique(df$biomarker) %in% names(thresholds)))

  out <- do.call(rbind, lapply(split(df, df$landmark_id), function(g) {
    do.call(rbind, lapply(names(thresholds), function(bm) {
      gbm <- g[g$biomarker == bm, ]
      if (nrow(gbm) == 0) return(NULL)
      rng <- thresholds[[bm]]
      data.frame(landmark_id = unique(g$landmark_id),
                 biomarker = bm,
                 teca = teca(gbm$time, gbm$value, rng[1], rng[2]))
    }))
  }))
  rownames(out) <- NULL

  if (spread) {
    out$col <- paste0(out$biomarker, "_area")
    wide <- reshape(out[, c("landmark_id", "col", "teca")],
                    idvar = "landmark_id", timevar = "col", direction = "wide")
    names(wide) <- sub("^teca\\.", "", names(wide))
    return(wide)
  }
  out
}
