# PSD calibration regularisation with a vanishing positive-eigenvalue threshold.
# Source this file and use siga_calibration_psd(M, n) at calibration calls.
# This is NOT an automatic modification of the public repository.
# For every sample size reported in the manuscript (n <= 2000), the threshold
# equals the original 1e-12 threshold, so this change does not alter that rule.

siga_calibration_psd <- function(M, n, tolerance = 1e-12) {
  M <- as.matrix(M)
  if (!is.numeric(M) || nrow(M) != ncol(M) || nrow(M) < 1L ||
      any(!is.finite(M))) {
    stop("M must be a nonempty finite numeric square matrix.", call. = FALSE)
  }
  if (length(n) != 1L || !is.finite(n) || n < 1 || n != floor(n)) {
    stop("n must be a positive integer sample size.", call. = FALSE)
  }
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be finite and nonnegative.", call. = FALSE)
  }
  M <- (M + t(M)) / 2
  ee <- eigen(M, symmetric = TRUE)
  threshold <- min(tolerance, 1 / n) * max(1, max(abs(ee$values)))
  values <- pmax(ee$values, 0)
  values[values < threshold] <- 0
  out <- ee$vectors %*% (values * t(ee$vectors))
  (out + t(out)) / 2
}

# In calibrate_pair_path_design_R(), n is already an argument. For example:
#   psi <- siga_calibration_psd((psi01 + psi02) / 2, n)
#   gamma <- siga_calibration_psd(gamma, n)
# If storing regularised diagnostic matrices, apply the same function to those
# matrices. Do not silently treat the two pair products as independent
# replicate observations. The estimator is the average of the TWO covariance
# estimators, not the covariance of the averaged product vectors.
