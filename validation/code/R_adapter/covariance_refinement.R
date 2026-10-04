# Drop-in covariance functions for the supplied base-R SIGA engine.
# No outcome model, true treatment effect, or true variance is used here.
# This file does not start a simulation or overwrite existing functions.

add_balancing_covariance <- function(calibration) {
  if (!inherits(calibration, "siga_pair_calibration")) stop("Invalid calibration class.")
  if (!(calibration$pbc > 0.5 && calibration$pbc < 1)) {
    stop("This adapter requires the stochastic range rule with 0.5 < pbc < 1.")
  }
  if (any(abs(calibration$weights - 1) > 1e-14)) stop("Equal weights are required.")
  patterns <- as.matrix(calibration$patterns)
  L <- t(cbind(1, 2 * patterns - 1))
  if (is.null(calibration$sigma_d[[1L]])) {
    stop("Need the first-copy covariance of D/sqrt(n), not Gamma alone.")
  }
  # sigma_d is estimated from the same allocation-only draws as Gamma.
  calibration$xi <- calibration$n * L %*% calibration$sigma_d[[1L]] %*% t(L)
  calibration$xi <- (calibration$xi + t(calibration$xi)) / 2
  calibration$balance_L <- L
  calibration
}

refined_imbalance_covariance <- function(counts, gamma, xi, L) {
  counts <- as.numeric(counts)
  if (anyNA(counts) || any(counts < 0) || sum(counts) <= 0) stop("Invalid counts.")
  sq <- sqrt(counts)
  old <- tcrossprod(sq) * gamma
  B <- sq * t(L)
  K <- crossprod(B)
  values <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  cutoff <- min(1e-12, 1 / sum(counts)) * max(1, max(values))
  if (min(values) <= cutoff) {
    return(list(Omega = old, original = old, fallback = TRUE, balance_error = NA_real_))
  }
  Q <- t(solve(K, t(B)))
  P <- Q %*% t(B)
  W <- diag(length(counts)) - P
  G <- W %*% gamma %*% t(W) + Q %*% xi %*% t(Q)
  O <- tcrossprod(sq) * G
  O <- (O + t(O)) / 2
  list(Omega = O, original = old, fallback = FALSE,
       balance_error = max(abs(L %*% O %*% t(L) - xi)))
}

siga_sampling_variance_refined <- function(score, X, calibration,
                                           truncate_kappa = TRUE) {
  if (is.null(calibration$xi)) calibration <- add_balancing_covariance(calibration)
  if (nrow(X) != calibration$n || ncol(X) != calibration$K) stop("Dimension mismatch.")
  # score_decomposition is the unchanged function in siga_pair_path_engine.R.
  dec <- score_decomposition(score, X, calibration$J)
  geometry <- refined_imbalance_covariance(dec$counts, calibration$gamma,
                                           calibration$xi, calibration$balance_L)
  Omega <- geometry$Omega
  active <- dec$counts > 0
  denom <- calibration$n - sum(active)
  tr <- sum(diag(Omega)[active] / dec$counts[active])
  kraw <- if (denom > 0) (calibration$n - tr) / denom else 0
  kappa <- if (truncate_kappa) max(kraw, 0) else kraw
  low <- drop(crossprod(dec$means, Omega %*% dec$means)) / 4
  within <- kappa * sum(dec$residual^2) / 4
  variance <- max(low + within, .Machine$double.eps)
  list(variance = variance, kappa = kappa, kappa_raw = kraw,
       low_component = low, within_component = within, decomposition = dec,
       rank_fallback = geometry$fallback, balance_error = geometry$balance_error)
}

# Example integration into an existing trial loop (not executed):
# cal <- add_balancing_covariance(cal)
# score <- make_boundary_scores(y, A, X, boundary)$unadjusted
# old <- siga_sampling_variance(score, X, cal)
# new <- siga_sampling_variance_refined(score, X, cal)
# dh <- estimate_stratum_effect_deviations(y, A, joint_stratum_id(X), cal$J, boundary)
# r_old <- corrected_variance(old, d = dh$d, calibration = cal)
# r_new <- corrected_variance(new, d = dh$d, calibration = cal)
# Regenerate the reference sequences ONCE, then compare all thresholds using
# that common Monte Carlo sample. Never replace dh$d by model_based_deviation()
# in the primary implementable-method comparison.
